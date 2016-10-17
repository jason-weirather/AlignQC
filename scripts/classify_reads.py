#!/usr/bin/env python
import sys, argparse, gzip

#from Bio.Format.GPD import GPD

def main(args):
  rinf = None
  if args.read_annotations[-3:] == '.gz':
    rinf = gzip.open(args.read_annotations)
  else:
    rinf = open(args.read_annotations)
  
  ginf = None
  if args.best_gpd[-3:] == '.gz':
    ginf = gzip.open(args.best_gpd)
  else:
    ginf = open(args.best_gpd)

  of = sys.stdout
  if args.output:
    if args.output[-3:] == '.gz':
      of = gzip.open(args.output,'w')
    else:
      of = open(args.output,'w')

  seen_reads = set()
  sys.stderr.write("traversing annotations\n")
  for line in rinf:
    f = line.rstrip().split("\t")
    read_name = f[1]
    gene_name = f[2]
    match_type = f[4]
    matching_exons = int(f[5])
    consecutive_exons = int(f[6])
    read_exons = int(f[7])
    if read_exons < 2:  continue
    seen_reads.add(read_name)
    if match_type == 'full':
      of.write(read_name+"\t"+gene_name+"\tfull"+"\n")
    elif consecutive_exons == read_exons:
      of.write(read_name+"\t"+gene_name+"\tsubset"+"\n")
    else:
      of.write(read_name+"\t"+gene_name+"\tnovel-isoform"+"\n")
  sys.stderr.write("traversing reads\n")
  z = 0
  for line in ginf:
    z+=1
    if z%1000 == 0: sys.stderr.write(str(z)+" reads processed\r")
    #gpd = GPD(line)
    f = line.rstrip().split("\t")
    if f[1] in seen_reads: continue
    if int(f[8]) < 2: continue
    of.write(f[1]+"\t"+"\tnovel-locus"+"\n")
  sys.stderr.write("\n")
def do_inputs():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('read_annotations',help="annotbest from alignqc is sufficient")
  parser.add_argument('best_gpd',help="best gpd from alignqc is sufficient")
  parser.add_argument('-o','--output',help="output the refined read classifications")
  args = parser.parse_args()
  return args

def external_cmd(cmd):
  cache_argv = sys.argv
  sys.argv = cmd.split()
  args = do_inputs()
  main(args)
  sys.argv = cache_argv
  return

if __name__=="__main__":
  args = do_inputs()
  main(args)

