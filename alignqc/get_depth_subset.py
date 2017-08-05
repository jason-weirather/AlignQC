#!/usr/bin/env python
"""Get the bed depths for each feature exon intron intergenic"""
import sys, argparse, re, gzip

from seqtools.range.multi import BedStream, intersect_range_array
from seqtools.stream import MultiLocusStream

def main(args):

  inf1 = None
  if re.search('\.gz$',args.depth_bed):
    inf1 = gzip.open(args.depth_bed)
  else:
    inf1 = open(args.depth_bed)
  inf2 = None
  if re.search('\.gz$',args.feature_bed):
    inf2 = gzip.open(args.feature_bed)
  else:
    inf2 = open(args.feature_bed)


  of = sys.stdout
  if args.output:
    if re.search('\.gz$',args.output):
      of = gzip.open(args.output,'w')
    else:
      of = open(args.output,'w')
  bs1 = BedStream(inf1)
  bs2 = BedStream(inf2)
  mls = MultiLocusStream([bs1,bs2])
  for overlapped in mls:
    [b1s,b2s] = overlapped.payload
    if len(b1s)==0 or len(b2s)==0: continue
    for b1 in b1s:
      m = intersect_range_array(b1,b2s,is_sorted=True)
      for rng in m:
        of.write("\t".join([str(x) for x in rng.get_bed_array()])+"\t"+b1.payload+"\n")
  of.close()

def external_cmd(cmd):
  cache_argv = sys.argv
  sys.argv = cmd
  args = do_inputs()
  main(args)
  sys.argv = cache_argv

def do_inputs():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('depth_bed',help="depth bed (sorted)")
  parser.add_argument('feature_bed',help="feature bed (sorted and merged)")
  parser.add_argument('-o','--output',help="write to output")
  args = parser.parse_args()
  return args

if __name__=="__main__":
  args = do_inputs()
  main(args)
