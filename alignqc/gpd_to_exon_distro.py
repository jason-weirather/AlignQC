#!/usr/bin/env python
import sys, argparse, re, gzip, inspect, os

from multiprocessing import Pool, cpu_count

from seqtools.format.gpd import GPDStream

def main(args):


  inf = None
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

  gs = GPDStream(inf)
  p = Pool(processes=args.threads)
  csize = 50000
  results = p.imap(func=do_gpd,iterable=generate_gpd(gs),chunksize=csize)
  for res in results:
    of.write(res)
  of.close()  

def do_gpd(gpd):
    return str(gpd.length)+"\t"+str(gpd.get_exon_count())+"\n"

def generate_gpd(gs):
  for gpd in gs:
    yield gpd

def external_cmd(cmd):
  cache_argv = sys.argv
  sys.argv = cmd
  args = do_inputs()
  main(args)
  sys.argv = cache_argv

def do_inputs():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="best gpd")
  parser.add_argument('-o','--output',help="write to output")
  parser.add_argument('--threads',default=cpu_count(),type=int,help="number of threads")
  args = parser.parse_args()
  return args  

if __name__=="__main__":
  args = do_inputs()
  main(args)
