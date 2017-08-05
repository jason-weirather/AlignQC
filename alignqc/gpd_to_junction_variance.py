#!/usr/bin/env python
"""measure distances of junctions to references junctions"""
import argparse, sys, os, gzip
from shutil import rmtree
from tempfile import mkdtemp, gettempdir
from time import sleep
from multiprocessing import cpu_count, Pool

from seqtools.format.gpd import GPD, GPDStream
from seqtools.range import GenomicRange, GenomicRangeFromString
from seqtools.range.multi import sort_ranges, BedArrayStream, BedStream
from seqtools.stream import MultiLocusStream

rcnt = 0

def main(args):

  # Start by reading in our reference GPD
  exon_start_strings = set()
  exon_end_strings = set()
  inf = None
  if args.reference[-3:] == '.gz':
    inf = gzip.open(args.reference)
  else:
    inf = open(args.reference)
  z = 0
  for line in inf:
    z += 1
    if z%1000 == 0: sys.stderr.write("ref transcripts: "+str(z)+" starts: "+str(len(exon_start_strings))+" ends: "+str(len(exon_end_strings))+"           \r")
    gpd = GPD(line)
    if gpd.get_exon_count() < 2: continue
    for j in gpd.junctions:
      exon_start_strings.add(j.right.get_range_string())
      exon_end_strings.add(j.left.get_range_string())
  inf.close()
  sys.stderr.write("\n")
  # Now convert each start or end to list of ranges to overlap
  sys.stderr.write("finding start windows\n")
  starts = get_search_ranges_from_strings(exon_start_strings,args)
  sh = BedArrayStream(starts)
  ofst = open(args.tempdir+'/starts.bed','w')
  for v in sh:
    ofst.write("\t".join([str(x) for x in v.get_bed_array()]+[v.payload.get_range_string()])+"\n")
  ofst.close()
  starts = None
  sh = None
  sys.stderr.write("finding end windows\n")
  ends = get_search_ranges_from_strings(exon_end_strings,args)
  eh = BedArrayStream(ends)
  ofen = open(args.tempdir+'/ends.bed','w')
  for v in eh:
    ofen.write("\t".join([str(x) for x in v.get_bed_array()]+[v.payload.get_range_string()])+"\n")
  ofen.close()
  # switch to file based operations to save some memory
  ends = None
  eh = None
  sleep(10) # give garbage collection time to flush these big memory objects
  sinf = open(args.tempdir+'/starts.bed')
  einf = open(args.tempdir+'/ends.bed')
  sh = BedStream(sinf)
  eh = BedStream(einf)
  # now stream in our reads
  sys.stderr.write("working through reads\n")
  inf = sys.stdin
  if args.input != '-':
    if args.input[-3:] == '.gz': inf = gzip.open(args.input)
    else: inf = open(args.input)
  gh = GPDStream(inf)
  mls = MultiLocusStream([gh,sh,eh])
  z = 0
  global rcnt
  rcnt = 0
  #start_distances = []
  #end_distances = []
  buffer = []
  max_buffer = 100
  tos = open(args.tempdir+'/starts.txt','w')
  toe = open(args.tempdir+'/ends.txt','w')
  p = Pool(processes=args.threads)
  csize=10
  results = p.imap_unordered(process_locus,gen_locus(mls,args),chunksize=csize)
  for r in results:
    if len(r[0]) > 0:
      tos.write("\n".join([str(x) for x in r[0]])+"\n")
    if len(r[1]) > 0:
      toe.write("\n".join([str(x) for x in r[1]])+"\n")
  tos.close()
  toe.close()
  inf.close()
  sys.stderr.write("\n")
  # now we have the distance, we don't actually know if a start is  a start or an end from what have
  distances = {}
  sys.stderr.write("Reading start distances\n")
  inf = open(args.tempdir+'/starts.txt')
  for line in inf:
    d = int(line.rstrip())
    if d not in distances: distances[d] = 0
    distances[d] += 1
  inf.close()
  sys.stderr.write("Reading end distances\n")
  inf = open(args.tempdir+'/ends.txt')
  for line in inf:
    d = int(line.rstrip())
    if d not in distances: distances[d] = 0
    distances[d] += 1
  inf.close()
  # now output results
  of = sys.stdout
  if args.output: of = open(args.output,'w')
  for d in sorted(distances.keys()):
    of.write(str(d)+"\t"+str(distances[d])+"\n")
  of.close()

  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

def gen_locus(mls,args):
  global rcnt
  z = 0
  for es in mls:
    z += 1
    if z%1000 == 0: sys.stderr.write(es.get_range_string()+" locus: "+str(z)+" reads: "+str(rcnt)+"        \r")
    if len(es.payload[0]) == 0: continue
    rcnt += len(es.payload[0])
    yield [es,args]

class Queue:
  def __init__(self,val):
    self.val = [val]
  def get(self):
    return self.pop(0)

def process_locus(vals):
    (es,args) = vals
    out_start_distances = []
    out_end_distances = []
    streams = es.payload
    reads = streams[0]
    starts = streams[1]
    for i in range(0,len(starts)):
      v = starts[i].payload
      starts[i].set_payload(GenomicRangeFromString(v))
    ends = streams[2]
    for i in range(0,len(ends)):
      v = ends[i].payload
      ends[i].set_payload(GenomicRangeFromString(v))
    #if len(starts) == 0 and len(ends) == 0: continue
    for read in reads:
      if read.get_exon_count() < 2: continue
      # multi exon
      for j in read.junctions:
        ex_start = j.right
        ex_end = j.left
        # we can look for evidence for each
        #sys.exit()
        evstart = [x.payload.start-ex_start.start for x in starts if x.overlaps(ex_start) and x.payload.distance(ex_start) <= args.window]
        if len(evstart) > 0:
          # get the best distance
          min_dist = args.window+10
          num = None
          for v in evstart:
            if abs(v) < min_dist:
              min_dist = abs(v)
              num = v
          if num is not None:
            out_start_distances.append(num)
        evend = [ex_end.end-x.payload.end for x in ends if x.overlaps(ex_end) and x.payload.distance(ex_end) <= args.window]
        if len(evend) > 0:
          # get the best distance
          min_dist = args.window+10
          num = None
          for v in evend:
            if abs(v) < min_dist:
              min_dist = abs(v)
              num = v
          if num is not None:
            out_end_distances.append(num)
    return [out_start_distances,out_end_distances]

def get_search_ranges_from_strings(position_strings,args):
  results = []
  for pstr in position_strings:
    rng = GenomicRangeFromString(pstr)
    rngw = rng.copy()
    rngw.start = max(rng.start-args.window,1)
    rngw.end =rng.start+args.window
    rngw.set_payload(rng)
    results.append(rngw)
  return sort_ranges(results)

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Given your SORTED read mappings in GPD format, how do the junction sites differ from the nearest known reference?",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="INPUT reads GPD or '-' for STDIN")
  parser.add_argument('-r','--reference',required=True,help="Reference GPD")
  parser.add_argument('-o','--output',help="OUTPUTFILE or STDOUT if not set")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")
  parser.add_argument('-w','--window',type=int,default=30,help="Window for how far to search for a nearby reference")

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
  cache_argv = sys.argv
  sys.argv = cmd
  args = do_inputs()
  main(args)
  sys.argv = cache_argv

if __name__=="__main__":
  args = do_inputs()
  main(args)
