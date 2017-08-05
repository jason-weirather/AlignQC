#!/usr/bin/env python
"""Get data about the lengths of the annotated sequences"""
import sys, argparse, gzip, re, os, inspect

from seqtools.format.gpd import GPDStream

def main(args):

  inf = None
  if re.search('\.gz',args.best_gpd):
    inf = gzip.open(args.best_gpd)
  else:
    inf = open(args.best_gpd)
  gs = GPDStream(inf)
  z = 0
  data = {}
  for gpd in gs:
    z += 1
    data[z] = [gpd.length,gpd.get_exon_count()]
    #gpd.length
  inf.close()
  inf = None
  if re.search('\.gz',args.best_annotation):
    inf = gzip.open(args.best_annotation)
  else:
    inf = open(args.best_annotation)
  done_reads = set()
  of = sys.stdout
  if args.output:
    if re.search('\.gz$',args.output):
      of = gzip.open(args.output,'w')
    else:
      of = open(args.output,'w')
  for line in inf:
    f = line.rstrip().split("\t")
    read_id = int(f[0])
    type = f[4]
    done_reads.add(read_id)
    of.write(type+"\t"+str(data[read_id][0])+"\t"+str(data[read_id][1])+"\n")
  for i in [x for x in range(1,z+1) if x not in done_reads]:
    of.write('unannotated'+"\t"+str(data[i][0])+"\t"+str(data[i][1])+"\n")
  of.close()

def do_inputs():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('best_gpd',help="Best alignments")
  parser.add_argument('best_annotation',help="Best annotations")
  parser.add_argument('-o','--output',help="output file")
  args = parser.parse_args()
  return args

def external_cmd(cmd):
  cache_argv = sys.argv
  sys.argv = cmd
  args = do_inputs()
  main(args)
  sys.argv = cache_argv
  
if __name__=="__main__":
  args = do_inputs()
  main(args)
