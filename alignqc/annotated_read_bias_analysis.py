#!/usr/bin/env python
"""Look for bias across the 5' to 3' of annotated transcripts"""
import sys, argparse, re, gzip, os, inspect
from subprocess import PIPE, Popen
from multiprocessing import Pool, cpu_count
from tempfile import mkdtemp, gettempdir
from shutil import rmtree

from seqtools.format.gpd import GPDStream
from seqtools.range import GenomicRange, GenomicRangeFromString
from seqtools.range.multi import ranges_to_coverage
from seqtools.statistics import average
from seqtools.stream import LocusStream, MultiLocusStream

def main(args):

  sort_annot(args)
  sort_ref(args)

  inf0 = open(args.tempdir+'/annot.sorted.txt')
  anns = AnnotStream(inf0)
  ls0 = LocusStream(anns)

  inf1 = open(args.tempdir+'/ref.sorted.gpd')
  gs1 = GPDStream(inf1)
  ls1 = LocusStream(gs1)

  sys.stderr.write("reading read genepred\n")

  inf2 = None
  if is_gzip(args.read_genepred):
    inf2 = gzip.open(args.read_genepred)
  else:
    inf2 = open(args.read_genepred)
  gs2 = GPDStream(inf2)
  ls2 = LocusStream(gs2)
  mls = MultiLocusStream([ls0,ls1,ls2])
  sys.stderr.write("stream loci\n")
  z = 0
  totals = []
  for l in mls:
    z += 1
    [l0,l1,l2] = l.payload
    if args.minimum_read_count > len(l0): continue
    if args.minimum_read_count > len(l2): continue
    rvals = do_locus(l0,l1,l2,args)
    for rval in rvals:
      totals.append(rval)
    sys.stderr.write(str(z)+" "+l.get_range_string()+" "+str([len(x) for x in l.payload])+" "+str(len(totals))+" proccessed        \r")
  sys.stderr.write("\n")
  #results = {}
  #for i in range(1,101):
  #  results[str(i)] = []
  read_total = 0
  ############
  outs = []
  for v in totals:
    if not v: continue
    bins = sorted([int(x) for x in v[0].keys()])
    outs.append([0 for x in range(1,101)])
    read_total+=v[1]
    for i in range(1,101):
      if str(i) in v[0]:
        #results[str(i)].append(v[0][str(i)])
        outs[-1][i-1] = v[0][str(i)]
  of = sys.stdout
  if args.output and re.search('\.gz',args.output):
    of = gzip.open(args.output,'w')
  elif args.output:
    of = open(args.output,'w')
  tot = len(outs)
  #for i in range(1,101):
  #  ostr = str(i)
  #  tot = len(results[str(i)])
  #  for j in results[str(i)]:
  #    ostr += "\t"+str(j)
  #  of.write(ostr+"\n")
  tnum = 0
  for o in outs:
    tnum += 1
    of.write(str(tnum)+"\t"+"\t".join([str(x) for x in o])+"\n")
  of.close()
  if args.output_counts:
    of = open(args.output_counts,'w')
    of.write(str(tot)+"\t"+str(read_total)+"\n")
    of.close()
  sys.stderr.write(str(tot)+" total transcripts \t"+str(read_total)+" total reads\n")

  if not args.specific_tempdir:
    rmtree(args.tempdir)

def spawn_jobs(mls,args):
  z = 0
  for l in mls:
    z += 1
    [l0,l1,l2] = l.payload
    if args.minimum_read_count > len(l0): continue
    if args.minimum_read_count > len(l2): continue
    vals = do_locus(l0,l1,l2,args)
    for v in vals:
      yield v

def do_locus(annots,refs,reads,args):
  read_to_tx = {}
  tx_to_read = {}
  for a in annots:
    for b in [x.get_value() for x in a.payload]:
      if b['matching_exon_count'] < args.minimum_matched_exons: continue
      if b['read_length'] < args.minimum_read_length: continue
      read_to_tx[b['read_name']] = b['tx_name']
      if b['tx_name'] not in tx_to_read: tx_to_read[b['tx_name']] = {}
      tx_to_read[b['tx_name']][b['read_name']] = ''
  for tx in tx_to_read.keys():
    if len(tx_to_read[tx]) < args.minimum_read_count:
      del tx_to_read[tx]
  tx_to_ref = {}
  for ref in refs:
    for gpd in ref.payload:
      tx = gpd.entries.name
      tx_to_ref[tx] = gpd
  for read in reads:
    for b in read.payload:
      if b.entries.name not in read_to_tx: continue
      tx = read_to_tx[b.entries.name]
      if tx not in tx_to_read: continue
      if tx not in tx_to_ref: continue
      tx_to_read[tx][b.entries.name] = b
  rvals = []
  for tx in tx_to_read:
    rvals.append(do_tx_line([tx_to_ref[tx],tx_to_read[tx].values(),args]))
  return rvals

def sort_ref(args):
  sys.stderr.write("Sorting in reference genePred\n")
  if args.threads > 1:
     cmd = ['sort','-S2G','-k3,3','-k5,5n','-k6,6n',
            '--parallel='+str(args.threads)]
  else:
     cmd = ['sort','-S2G','-k3,3','-k5,5n','-k6,6n']
  of = open(args.tempdir+'/ref.sorted.gpd','w')
  p = Popen(cmd,stdin=PIPE,stdout=of)
  refgpd = {}
  if args.ref_genepred[-3:] == '.gz':
     inf = gzip.open(args.ref_genepred)
  else:
     inf = open(args.ref_genepred)
  #gs = GPDStream(inf)
  z = 0
  for line in inf:
    z += 1
    if z%1000==0: sys.stderr.write(str(z)+"       \r")
    #if z not in refcnt: continue
    #if refcnt[z] < args.minimum_read_count: continue
    p.stdin.write(line)
    #refgpd[z] = gpd
  p.communicate()
  sys.stderr.write("\n")
  inf.close()

def sort_annot(args):
  sys.stderr.write("Sorting read annotations\n")
  cmd = ['sort','-S2G','-k1,1','-k2,2n','-k3,3n']
  cmd2 = 'cut -f 4-'
  of0 = open(args.tempdir+'/annot.sorted.txt','w')
  p1 = Popen(cmd2.split(),stdin=PIPE,stdout=of0)
  p0 = Popen(cmd,stdin=PIPE,stdout=p1.stdin)
  inf = None
  if is_gzip(args.annotations):
    inf = gzip.open(args.annotations)
  else:
    inf = open(args.annotations)
  k = 0
  for line in inf:
    k+=1
    if k%1000==0: sys.stderr.write(str(k)+"       \r")
    f = line.rstrip().split("\t")
    r = GenomicRangeFromString(f[13])
    #r.set_payload(parse_annot(f))
    p0.stdin.write(r.chr+"\t"+str(r.start)+"\t"+str(r.end)+"\t"+line)
  sys.stderr.write("\n")
  of0.close()
  p0.communicate()
  p1.communicate()
  inf.close()

def do_tx_line(vals):
    (ref_gpd,reads,args) = vals
    allbits = []
    read_count = 0
    outrange = reads[-1].get_range()
    for read in reads:
      if not args.allow_overflowed_matches and read.get_range().start < ref_gpd.get_range().start: continue
      if not args.allow_overflowed_matches and read.get_range().end > ref_gpd.get_range().end: continue
      v = ref_gpd.union(read)
      for e in [x.rng for x in v.exons]: allbits.append(e)
      read_count += 1
    if len(allbits)==0: return None
    if read_count < args.minimum_read_count: return None
    cov = ranges_to_coverage(allbits)
    #print [x.get_payload() for x in cov]
    curr = 0
    bps = []
    for i in range(0,ref_gpd.length):
      bps.append(0)
    for rng1 in [x.rng for x in ref_gpd.exons]:
      overs =  [[z[0],z[1].payload] for z in [[y.union(rng1),y] for y in cov] if z[0]]
      for ov in overs:
        dist1 = ov[0].start - rng1.start+curr
        dist2 = ov[0].end - rng1.start+curr
        for i in range(dist1,dist2+1):
          bps[i]+=ov[1]
      curr+=rng1.length
    trimmedbps = bps
    if args.only_covered_ends:
      start = 0
      finish = len(bps)-1
      for i in range(0,len(bps)):
        if bps[i] != 0: 
          start = i
          break
      for i in reversed(range(0,len(bps))):
        if bps[i] != 0: 
          finish = i
          break
      trimmedbps = bps[start:finish+1]
    exp = float(sum(trimmedbps))/float(len(trimmedbps))
    if ref_gpd.get_strand()=='-': trimmedbps = list(reversed(trimmedbps))
    if len(trimmedbps) < args.minimum_read_count: return None
    #bin the results
    vals = {}
    for dat in [[str(1+int(100*float(i)/float(len(trimmedbps)))),float(trimmedbps[i])/float(read_count)] for i in range(0,len(trimmedbps))]:
      if dat[0] not in vals: vals[dat[0]] = []
      vals[dat[0]].append(dat[1])
    for num in vals:
      vals[num] = average(vals[num])
    return [vals, read_count, exp, len(trimmedbps),ref_gpd.get_exon_count(),outrange.get_range_string()]


def is_gzip(name):
  if re.search('\.gz$',name): return True
  return False

def do_inputs():
  parser = argparse.ArgumentParser(description="Generate a coverage per bin over the length of a molecule to observe bias in the 5' to 3' mapping of reads",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('read_genepred',help="Input genepred")
  parser.add_argument('ref_genepred',help="Reference genepred")
  parser.add_argument('annotations',help="Input annotations")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="Threads count")
  parser.add_argument('--full',action='store_true',help="only consider full length matched reads")
  parser.add_argument('--only_covered_ends',action='store_true',help="remove ends with zero coverage")
  parser.add_argument('--allow_overflowed_matches',action='store_true',help="by default we don't consider matches that arent fully within the bounds of their annotated transcript.")
  parser.add_argument('--minimum_read_count',type=int,default=5,help="minimum number of reads")
  parser.add_argument('--minimum_read_length',type=int,default=100,help="at least this many bp")
  parser.add_argument('--minimum_matched_exons',type=int,default=2,help="require reads matched at least this many exons")
  parser.add_argument('-o','--output',help="write output to file")
  parser.add_argument('--output_counts',help="write number of transcripts and reads used")
  # Temporary working directory step 1 of 3 - Definition
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is ")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used")
  args = parser.parse_args()
  setup_tempdir(args)
  return args

def parse_annot(f):
  res={
    'read_line':int(f[0]),\
    'read_name':f[1],\
    'gene_name':f[2],\
    'tx_name':f[3],\
    'type':f[4],\
    'matching_exon_count':int(f[5]),\
    'consecutive_exons':int(f[6]),\
    'read_exons':int(f[7]),\
    'tx_exons':int(f[8]),\
    'overlap':int(f[9]),\
    'read_length':int(f[10]),\
    'tx_length':int(f[11])}
  return res


class Annot:
  def __init__(self,line):
    f = line.rstrip().split("\t")
    self.value = parse_annot(f)
    self.range = GenomicRangeFromString(f[13])
  def get_range(self):
    return self.range
  def get_value(self):
    return self.value

class AnnotStream:
  def __init__(self,fh):
    self.fh = fh

  def read_entry(self):
    line = self.fh.readline()
    if not line: return False
    a = Annot(line)
    return a

  def next(self):
    r = self.read_entry()
    if not r: raise StopIteration
    else:
      return r

  def __iter__(self):
    return self

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
  cache_argv = sys.argv
  sys.argv = cmd
  args = do_inputs()
  main(args)
  sys.argv = cache_argv

if __name__=="__main__":
  args = do_inputs()
  main(args)
