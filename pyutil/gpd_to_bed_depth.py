#!/usr/bin/python
import sys, argparse, gzip, re
from Bio.Format.GPD import GPDStream
from Bio.Stream import LocusStream
from Bio.Range import ranges_to_coverage

from multiprocessing import cpu_count, Pool

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
  p = Pool(processes=args.threads)
  loci = LocusStream(GPDStream(inf))
  csize=100
  results = p.imap(func=do_locus,iterable=generate_gpd(loci),chunksize=csize)
  for covs in results:
    for cov in covs:
      of.write(cov)
  of.close()
  inf.close()

def do_locus(locus):
  exranges = []
  for entry in locus.get_payload():
    for exon in entry.exons:
      exranges.append(exon.get_range())
  covs = ranges_to_coverage(exranges)
  output = []
  for cov in covs:    
    output.append("\t".join([str(x) for x in cov.get_bed_coordinates()])+"\t"+str(+cov.get_payload())+"\n")
  return output

def generate_gpd (loci):
  for locus in loci:
    yield locus

def do_inputs():
  parser = argparse.ArgumentParser(description="Convert sorted gpd file to bed depth",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN")
  parser.add_argument('--threads',default=cpu_count(),type=int,help="number of threads to use (default cpu_count)")
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
