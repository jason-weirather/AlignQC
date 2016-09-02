#!/usr/bin/python
import sys, argparse, gzip, re
from Bio.Format.GPD import GPDStream
from Bio.Stream import LocusStream
from Bio.Range import ranges_to_coverage

def main(args):
  
  inf = sys.stdin
  if args.input != '-':
    if re.search('\.gz$',args.input):
      inf = gzip.open(args.input)
    else:
      inf = open(args.input)
  of = sys.stdout
  if args.output:
    if re.search('\.gz$',args.output):
      of = gzip.open(args.output,'w')
    else:
      of = open(args.output,'w')
  loci = LocusStream(GPDStream(inf))
  for locus in loci:
    exranges = []
    for entry in locus.get_payload():
      for exon in entry.exons:
        exranges.append(exon.get_range())
    covs = ranges_to_coverage(exranges)    
    for cov in covs:    
      of.write("\t".join([str(x) for x in cov.get_bed_coordinates()])+"\t"+str(+cov.get_payload())+"\n")
  of.close()
  inf.close()

def do_inputs():
  parser = argparse.ArgumentParser(description="Convert sorted gpd file to bed depth",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN")
  parser.add_argument('-o','--output',help="OUTPUT or dont set for STDOUT")
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
