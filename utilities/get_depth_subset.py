#!/usr/bin/python
import sys, argparse, re, gzip, inspect, os

#bring in the folder to the path for our utilities
pythonfolder_loc = "../pylib"
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe() ))[0],pythonfolder_loc)))
if cmd_subfolder not in sys.path:
  sys.path.insert(0,cmd_subfolder)

from Bio.Range import BedStream, union_range_array
from Bio.Stream import MultiLocusStream

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
    [b1s,b2s] = overlapped.get_payload()
    if len(b1s)==0 or len(b2s)==0: continue
    for b1 in b1s:
      m = union_range_array(b1,b2s,is_sorted=True)
      for rng in m:
        of.write("\t".join([str(x) for x in rng.get_bed_array()])+"\t"+b1.get_payload()+"\n")
  of.close()

def external_cmd(cmd):
  cache_argv = sys.argv
  sys.argv = cmd.split()
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
