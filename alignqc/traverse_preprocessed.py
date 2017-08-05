#!/usr/bin/env python
import argparse, sys, os, gzip, itertools, inspect, pickle, zlib, base64
from shutil import rmtree
from multiprocessing import cpu_count, Pool, Lock
from tempfile import mkdtemp, gettempdir
from subprocess import Popen, PIPE

from seqtools.range import GenomicRangeFromString
from seqtools.range.multi import ranges_to_coverage
from seqtools.format.gpd import GPD, SortedOutputFile as SortedGPDOutputFile
from seqtools.stream import GZippedOutputFile

## The purpose of this script is to read through a bam alignment and record as much information as possible from it.  ##
## The bam should be indexed ahead of time in our index format.

glock = Lock()

g_done = None
g_lines = None
best_gpd = None
chimera_gpd = None
technical_chimera_gpd = None
technical_atypical_chimera_gpd = None
gapped_gpd = None
gpd_total = None
g_lengths = None

def process_buffer_output(results):
  #print results
  global glock
  glock.acquire()
  global best_gpd
  global g_lengths
  global chimera_gpd
  global gapped_gpd
  global technical_chimera_gpd
  global technical_atypical_chimera
  global g_done
  global g_lines
  g_done += len(results['lengths'])
  g_lines += results['buffer_quantity']
  sys.stderr.write(str(g_lines)+" alignments   "+str(g_done)+" reads        \r")
  for v in results['best']:
    best_gpd.write(v+"\n")
  for v in results['lengths']:
    g_lengths.write(v+"\n")
  for v in results['chimera']:
    chimera_gpd.write(v+"\n")
  for v in results['gapped']:
    gapped_gpd.write(v+"\n")
  for v in results['technical_chimera']:
    technical_chimera_gpd.write(v+"\n")
  for v in results['technical_atypical_chimera']:
    technical_atypical_chimera_gpd.write(v+"\n")

  #'original_count':original_count,'gapped_count':gapped_count,'technical_atypical_chimera_count':technical_atypical_chimera_count,'techinical_chimera_count':technical_chimera_count,'chimera_count':chimera_count}

  glock.release()

def do_buffer(buffer,args):
  lengths = []
  best = []
  chimera = []
  technical_chimera = []
  gapped = []
  technical_atypical_chimera = []
  chimera_count = 0
  technical_chimera_count = 0
  gapped_count = 0
  technical_atypical_chimera_count = 0
  original_count = 0
  results = []
  buffer_quantity = 0
  for qname in sorted(buffer.keys()):
    buffer_quantity += len(buffer[qname])
    dat = [pickle.loads(zlib.decompress(base64.b64decode(x))) for x in buffer[qname]]
    for i in range(0,len(dat)): 
      if dat[i]['aligned_bases'] > 0:
        dat[i]['tx'] = GPD(dat[i]['tx'])
        dat[i]['qrng'] = GenomicRangeFromString(dat[i]['qrng'])
      else:
        dat[i]['tx'] = None
        dat[i]['qrng'] = None
    unaligned = [x for x in dat if x['aligned_bases'] == 0]
    aligned = [x for x in dat if x['aligned_bases'] > 0]
    #now dat should be set up with anything we need
    if len(aligned)==0:
      lengths.append(qname+"\tunaligned\t0\t0\t"+str(unaligned[0]['qlen']))
      continue
    #if len([x for x in aligned if x['flag'] & 2304==0])==0:
    #  lengths.append(qname+"\tunaligned\t0\t0\t"+str(aligned[0]['qlen']))
    #  continue
    if len(aligned)==1:
      if dat[0]['aligned_bases'] > 0:
        lengths.append(qname+"\toriginal\t"+str(dat[0]['aligned_bases'])+"\t"+str(dat[0]['aligned_bases'])+"\t"+str(dat[0]['qlen']))
        best.append(dat[0]['tx'].get_gpd_line())
        original_count += 1
      continue
    # we have multiple paths too look for
    qlen = max([x['qlen'] for x in aligned])
    #print aligned
    #print [x['tx'].get_gene_name() for x in aligned]
    best_possible = [i for i in range(0,len(aligned)) if aligned[i]['flag'] & 2304 == 0]
    best_ind = 0
    if len(best_possible) == 0:
      ## only secondary alignments to look at
      longest = 0
      for i in range(0,len(aligned)):
        if aligned[i]['aligned_bases'] > longest:
          longest = aligned[i]['aligned_bases']
          longind = i
    else:
      best_ind= best_possible[0]
    best.append(dat[best_ind]['tx'].get_gpd_line())
    v = check_paths(aligned,best_ind,args)

    o_qlen = dat[best_ind]['qrng'].length
    v_qlen = v['qlen']
    if v['type'] == 'chimera':
      chimera_count += 1
      for p in v['path']:
        chimera.append(p['tx'].get_gpd_line())
    elif v['type'] == 'self-chimera':
      technical_chimera_count += 1
      for p in v['path']:
        technical_chimera.append(p['tx'].get_gpd_line())
    elif v['type'] == 'self-chimera-atypical':
      technical_atypical_chimera_count += 1
      for p in v['path']:
        technical_atypical_chimera.append(p['tx'].get_gpd_line())
    elif v['type'] == 'gapped':
      gapped_count += 1
      for p in v['path']:
        gapped.append(p['tx'].get_gpd_line())
    elif v['type'] == 'original':
      original_count +=1
    else:
      sys.stderr.write("WARNING unaccounted for type\n")
    lengths.append(qname+"\t"+v['type']+"\t"+str(o_qlen)+"\t"+str(v_qlen)+"\t"+str(max([x['qlen'] for x in aligned])))
  return {'lengths':lengths,'gapped':gapped,'chimera':chimera,'best':best,'technical_chimera':technical_chimera,'technical_atypical_chimera':technical_atypical_chimera,'original_count':original_count,'gapped_count':gapped_count,'technical_atypical_chimera_count':technical_atypical_chimera_count,'technical_chimera_count':technical_chimera_count,'chimera_count':chimera_count,'buffer_quantity':buffer_quantity}


def main(args):
  chunksize = args.chunk_size
  inf = gzip.open(args.input)
  args.output = args.output.rstrip('/')
  if not os.path.exists(args.output):
    os.makedirs(args.output)
  buffer = {}
  prev = None
  global g_done
  global g_lines
  g_done = 0
  g_lines = 0
  global best_gpd
  best_gpd = SortedGPDOutputFile(args.output+'/best.sorted.gpd.gz',tempdir=args.tempdir)
  global g_lengths
  g_lengths = GZippedOutputFile(args.output+'/lengths.txt.gz')
  #cmd = "gzip"
  #lof = open(args.output+'/lengths.txt.gz','w')
  #plen = Popen(cmd.split(),stdout=lof,stdin=PIPE,close_fds=True)
  #g_lengths = plen.stdin
  global chimera_gpd
  chimera_gpd = GZippedOutputFile(args.output+'/chimera.gpd.gz')
  global technical_chimera_gpd
  technical_chimera_gpd = GZippedOutputFile(args.output+'/technical_chimeras.gpd.gz')
  global technical_atypical_chimera_gpd
  technical_atypical_chimera_gpd = GZippedOutputFile(args.output+'/technical_atypical_chimeras.gpd.gz')
  global gapped_gpd
  gapped_gpd = GZippedOutputFile(args.output+'/gapped.gpd.gz')
  global gpd_total
  gpd_total = {'original_count':0,'gapped_count':0,'technical_atypical_chimera_count':0,'techinical_chimera_count':0,'chimera_count':0,'unaligned':0}

  if args.threads > 1:
    p = Pool(processes=args.threads)
  z = 0
  for line in inf:
    (qname, data) = line.rstrip().split("\t")
    if qname!=prev:
      buftot = len(buffer.keys())
      if buftot >= chunksize:
        if args.threads > 1:
          p.apply_async(do_buffer,args=(buffer,args),callback=process_buffer_output)
        else:
          r = do_buffer(buffer,args)
          process_buffer_output(r)
        buffer = {}
    if qname not in buffer: buffer[qname] = []
    buffer[qname].append(data)
    prev = qname
  if len(buffer.keys()) > 0:
    if args.threads > 1:
      p.apply_async(do_buffer,args=(buffer,args),callback=process_buffer_output)
    else:
      r = do_buffer(buffer,args)
      process_buffer_output(r)
  if args.threads > 1:
    p.close()
    p.join()
  sys.stderr.write("\n")
  best_gpd.close()
  chimera_gpd.close()
  technical_chimera_gpd.close()
  technical_atypical_chimera_gpd.close()
  gapped_gpd.close()
  g_lengths.close()
  #plen.communicate()
  #lof.close()
  inf.close()

  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

# path
# aligned_bases - bases aligned not counting any deltions or insertions
# indecies - 
# type - original/chimera/self-chimera/gapped
# qlen - range spanned by query alignments
def check_paths(path_data,best_ind,args):
  #other_inds = [x for x in range(0,len(path_data)) if x != best_ind]
  possibles = get_index_sets(len(path_data))
  new_best = [path_data[best_ind]]
  new_bases = path_data[best_ind]['aligned_bases']
  new_inds = set([best_ind])
  new_type = 'original'
  new_qlen = path_data[best_ind]['qrng'].length
  for possible_path in possibles:
    if best_ind not in possible_path: continue # only consider path sets that have our best index in it
    res = evaluate_path(path_data,possible_path,best_ind,args)
    if res['any']:
      bases = sum([x['aligned_bases'] for x in res['path']])
      if bases > new_bases:
        new_best = res['path']
        new_bases = bases
        new_inds = set(possible_path)
        qrngs = [res['path'][0]['qrng']]
        for i in range(1,len(res['path'])):
          if qrngs[-1].overlaps(res['path'][i]['qrng']):
            qrngs[-1] = qrngs[-1].merge(res['path'][i]['qrng'])
          else: qrngs.append(res['path'][i]['qrng'])
        #new_qlen = sum([x.length() for x in qrngs])
        new_qlen = sum([x.length for x in ranges_to_coverage(qrngs)])
        if res['gapped']: new_type = 'gapped'
        elif res['chimera']: new_type = 'chimera'
        elif res['self-chimera']: new_type = 'self-chimera'
        elif res['self-chimera-atypical']: new_type = 'self-chimera-atypical'
        else:
          sys.stderr.write("WARNING: Unaccounted for type\n")
  return {'path':new_best, 'aligned_bases':new_bases, 'indecies':new_inds,'type':new_type,'qlen':new_qlen}
  #print path_data[best_ind]

# Create a dictionary with the follwing information
# path: a list of alignments order by query placement
# gapped: is it a gapped alignment
# chimera: is it a fusion-like 
# self-chimera: is it a + - of an overlapping target sequence
def evaluate_path(path_data,possible_path,best_ind,args):
  pord = sorted([path_data[i] for i in possible_path],key=lambda x: x['qrng'].start)
  best_bases = path_data[best_ind]['aligned_bases']
  bases = sum([x['aligned_bases'] for x in pord])
  res = {'path':pord,'gapped':False,'chimera':False,'self-chimera':False,'self-chimera-atypical':False,'any':False}
  if len(path_data) <= 1: return res
  if bases+bases*args.required_fractional_improvement < best_bases:  
    return res
  for p in pord:
    if p['aligned_bases'] < args.min_aligned_bases: return res
  # check for query overlaps ... not a useful build
  for i in range(0,len(pord)):
    for j in range(i+1,len(pord)):
      if args.max_query_gap:
        if pord[i]['qrng'].distance(pord[j]['qrng']) > args.max_query_gap: return res
      if pord[i]['qrng'].overlap_size(pord[j]['qrng']) > args.max_query_overlap:
        return res

  chrcount = len(set([x['tx'].range.chr for x in pord]))

  # check for target overlaps ... not gapped or chimera but maybe self-chimera
  for i in range(0,len(pord)):
    for j in range(i+1,len(pord)):  
      if pord[i]['tx'].overlap_size(pord[j]['tx']) > args.max_target_overlap:
        #res['gapped'] = False
        #res['chimera'] = False
        if pord[i]['tx'].strand != pord[j]['tx'].strand and chrcount == 1:
          res['self-chimera'] = True
          res['any'] = True
        else: 
          res['self-chimera-atypical'] = True
          res['any'] = True
        return res

  for i in range(0,len(pord)):
    for j in range(i+1,len(pord)):  
      if args.max_target_gap:
        dist = pord[i]['tx'].range.distance(pord[j]['tx'].range)
        if dist > args.max_target_gap or dist == -1: 
          res['chimera'] = True
          res['gapped'] = False
          res['any'] = True
  if len(pord) > 1 and not res['self-chimera'] and not res['chimera']:
    res['gapped'] = True
    res['any'] = True
  return res

def get_index_sets(indlen):
  r = []
  inds = range(0,indlen)
  for l in range(1,len(inds)+1):
    for subset in itertools.combinations(inds,l):
      r.append(subset)
  return r

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="use bam_preprocess result as input")
  parser.add_argument('-o','--output',help="OUTPUTDIR",required=True)
  parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")
  parser.add_argument('--minimum_intron_size',type=int,default=68)
  parser.add_argument('--chunk_size',type=int,default=10000,help="number of reads to do at once")

  # Arguments for finding alternate multi alignment paths
  parser.add_argument('--min_aligned_bases',type=int,default=50,help="Don't consider very short alignments")
  parser.add_argument('--max_query_overlap',type=int,default=10,help="Consider two alignments incompatible if greater overlap than this")
  parser.add_argument('--max_target_overlap',type=int,default=10,help="Consider two alignments incompatible if greater overlap than this")
  parser.add_argument('--max_target_gap',type=int,default=500000,help="Not a gapped alignment if gap is greater than this")
  parser.add_argument('--max_query_gap',type=int,default=500000,help="Consider a gapped alignment incompatible if greater thant this")
  parser.add_argument('--required_fractional_improvement',type=float,default=0.2,help="combination path should be this much better than best single alignment")  

  # Temporary working directory step 1 of 3 - Definition
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")
  args = parser.parse_args()

  # Temporary working directory step 2 of 3 - Creation
  setup_tempdir(args)
  return args

def setup_tempdir(args):
  if args.specific_tempdir:
    if not os.path.exists(args.specific_tempdir):
      os.makedirs(args.specific_tempdir.rstrip('/'))
    args.tempdir = args.specific_tempdir.rstrip('/')
    if not os.path.exists(args.specific_tempdir.rstrip('/')):
      sys.stderr.write("ERROR: Problem creating temporary directory\n")
      sys.exit()
  else:
    args.tempdir = mkdtemp(prefix="weirathe.",dir=args.tempdir.rstrip('/'))
    if not os.path.exists(args.tempdir.rstrip('/')):
      sys.stderr.write("ERROR: Problem creating temporary directory\n")
      sys.exit()
  if not os.path.exists(args.tempdir):
    sys.stderr.write("ERROR: Problem creating temporary directory\n")
    sys.exit()
  return 

def external_cmd(cmd):
  #need to save arguments
  cache_argv = sys.argv
  sys.argv = cmd
  args = do_inputs()
  main(args)
  #need to set the arguments back to what they were
  sys.argv = cache_argv
  return

if __name__=="__main__":
  #do our inputs
  args = do_inputs()
  main(args)
