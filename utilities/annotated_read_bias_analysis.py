#!/usr/bin/python
import sys, argparse, re, gzip, os, inspect

#bring in the folder to the path for our utilities
pythonfolder_loc = "../pylib"
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe() ))[0],pythonfolder_loc)))
if cmd_subfolder not in sys.path:
  sys.path.insert(0,cmd_subfolder)

from Bio.Format.GPD import GPDStream
from Bio.Range import GenomicRange, ranges_to_coverage
from Bio.Statistics import average

def main(args):

  sys.stderr.write("Reading in reference genePred\n")
  refgpd = {}
  inf = open(args.ref_genepred)
  gs = GPDStream(inf)
  z = 0
  for gpd in gs:
    z += 1
    refgpd[z] = gpd
  inf.close()

  sys.stderr.write("Reading in read annotations\n")
  inf = None
  if is_gzip(args.annotations):
    inf = gzip.open(args.annotations)
  else:
    inf = open(args.annotations)
  reflocs = {}
  rline = {}
  for line in inf:
    f = line.rstrip().split("\t")
    res={'read_line':int(f[0]),\
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
    'tx_length':int(f[11]),\
    'read_range':GenomicRange(range_string=f[12]),\
    'tx_range':GenomicRange(range_string=f[13]),\
    'ref_line':int(f[14])}
    if res['ref_line'] not in reflocs: reflocs[res['ref_line']] = []
    reflocs[res['ref_line']].append(res)
    if args.full and res['type'] != 'full': continue
    if args.minimum_matched_exons > res['matching_exon_count']: continue
    rline[res['read_line']] = res
  inf.close()

  sys.stderr.write("reading read genepred\n")
  inf = None
  if is_gzip(args.read_genepred):
    inf = gzip.open(args.read_genepred)
  else:
    inf = open(args.read_genepred)
  gs = GPDStream(inf)
  z = 0
  originals = {}
  for gpd in gs:
    z+=1
    if z not in rline: continue
    refline = rline[z]['ref_line']
    if refline not in originals: originals[refline] = {}
    originals[refline][z] = gpd
  inf.close()
  results = {}
  for i in range(1,101):
    results[str(i)] = []
  read_total = 0
  outs = {}
  for tx_line in originals:
    ref_gpd = refgpd[tx_line]
    annots = reflocs[tx_line]
    reads = originals[tx_line].values()
    v = do_tx_line(ref_gpd,annots,reads,args)
    if not v: continue
    tname = ref_gpd.get_transcript_name()
    bins = sorted([int(x) for x in v[0].keys()])
    outs[tname] = [0 for x in range(1,101)]
    read_total+=v[1]
    for i in range(1,101):
      if str(i) in v[0]:
        results[str(i)].append(v[0][str(i)])
        outs[tname][i-1] = v[0][str(i)]
      #else:
      #  results[str(i)].append(0)
  of = sys.stdout
  if args.output and re.search('\.gz',args.output):
    of = gzip.open(args.output,'w')
  elif args.output:
    of = open(args.output,'w')
  tot = len(outs.keys())
  #for i in range(1,101):
  #  ostr = str(i)
  #  tot = len(results[str(i)])
  #  for j in results[str(i)]:
  #    ostr += "\t"+str(j)
  #  of.write(ostr+"\n")
  for tname in outs:
    of.write(tname+"\t"+"\t".join([str(x) for x in outs[tname]])+"\n")
  of.close()
  if args.output_counts:
    of = open(args.output_counts,'w')
    of.write(str(tot)+"\t"+str(read_total)+"\n")
    of.close()
  sys.stderr.write(str(tot)+" total transcripts \t"+str(read_total)+" total reads\n")

def do_tx_line(ref_gpd,annots,reads,args):
    allbits = []
    read_count = 0
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
    for i in range(0,ref_gpd.get_length()):
      bps.append(0)
    for rng1 in [x.rng for x in ref_gpd.exons]:
      overs =  [[z[0],z[1].get_payload()] for z in [[y.union(rng1),y] for y in cov] if z[0]]
      for ov in overs:
        dist1 = ov[0].start - rng1.start+curr
        dist2 = ov[0].end - rng1.start+curr
        for i in range(dist1,dist2+1):
          bps[i]+=ov[1]
      curr+=rng1.length()
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
    return [vals, read_count, exp, len(trimmedbps),ref_gpd.get_exon_count()]


def is_gzip(name):
  if re.search('\.gz$',name): return True
  return False

def do_inputs():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('read_genepred',help="Input genepred")
  parser.add_argument('ref_genepred',help="Reference genepred")
  parser.add_argument('annotations',help="Input annotations")
  parser.add_argument('--full',action='store_true',help="only consider full length matched reads")
  parser.add_argument('--only_covered_ends',action='store_true',help="remove ends with zero coverage")
  parser.add_argument('--allow_overflowed_matches',action='store_true',help="by default we don't consider matches that arent fully within the bounds of their annotated transcript.")
  parser.add_argument('--minimum_read_count',type=int,default=5,help="minimum number of reads")
  parser.add_argument('--minimum_read_length',type=int,default=100,help="at least this many bp")
  parser.add_argument('--minimum_matched_exons',type=int,default=2,help="require reads matched at least this many exons")
  parser.add_argument('-o','--output',help="write output to file")
  parser.add_argument('--output_counts',help="write number of transcripts and reads used")
  args = parser.parse_args()
  return args

def external_cmd(cmd):
  cache_argv = sys.argv
  sys.argv = cmd.split()
  args = do_inputs()
  main(args)
  sys.argv = cache_argv

if __name__=="__main__":
  args = do_inputs()
  main(args)
