#!/usr/bin/python
import sys, argparse, gzip

from Bio.Format.GPD import GPD

def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('read_annotations',help="annotbest from alignqc is sufficient")
  parser.add_argument('best_gpd',help="best gpd from alignqc is sufficient")
  args = parser.parse_args()
  
  rinf = None
  if args.read_annotations[-3:] == '.gz':
    rinf = gzip.open(args.read_annotations)
  else:
    rinf = open(args.read_annotations)
  
  ginf = None
  if args.read_annotations[-3:] == '.gz':
    ginf = gzip.open(args.best_gpd)
  else:
    ginf = open(args.best_gpd)
  seen_reads = set()
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
      print read_name+"\t"+gene_name+"\tfull"
    elif consecutive_exons == read_exons:
      print read_name+"\t"+gene_name+"\tsubset"
    else:
      print read_name+"\t"+gene_name+"\tnovel-isoform"
  for line in ginf:
    gpd = GPD(line)
    if gpd.get_transcript_name() in seen_reads: continue
    if gpd.get_exon_count() < 2: continue
    print gpd.get_transcript_name()+"\t"+"\tnovel-locus"
if __name__=="__main__":
  main()
