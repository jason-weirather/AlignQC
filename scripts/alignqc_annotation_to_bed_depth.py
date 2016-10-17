#!/usr/bin/env python
import sys, argparse, gzip
from Bio.Format.GPD import GPD
from Bio.Range import ranges_to_coverage

def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN")
  parser.add_argument('genepred',help="the genepred used for this alignqc")
  parser.add_argument('--min_exons',type=int,default=1,help="At least this number of exons")
  parser.add_argument('--full',action='store_true',help="only use full matches")
  parser.add_argument('-o','--output',help="OUTPUT file or nothing for STDOUT")
  args = parser.parse_args()
  
  inf = sys.stdin
  if args.input != '-':
    if args.input[-3:]=='.gz':
      inf = gzip.open(args.input)
    else: inf = open(args.input)
  genes = {}
  sys.stderr.write("Reading annotation file\n")
  for line in inf:
    f = line.rstrip().split("\t")  
    gene = f[2]
    tx = f[3]
    type = f[4]
    if args.full and type != 'full': continue
    if gene not in genes:
      genes[gene] = {}
      genes[gene]['transcripts'] = {}
      genes[gene]['cnt'] = 0
    if tx not in genes[gene]['transcripts']:
      genes[gene]['transcripts'][tx] = 0
    genes[gene]['cnt'] += 1
    genes[gene]['transcripts'][tx] += 1
  inf.close()

  txs = {}
  sys.stderr.write("Reading genepred file\n")
  z = 0
  with open(args.genepred) as inf:
    for line in inf:
      z +=1
      if z%1000==0: sys.stderr.write(str(z)+"   \r")
      gpd = GPD(line)
      exs = []
      for ex in gpd.exons:
        exs.append(ex.get_range())
      txs[gpd.get_transcript_name()] = exs
  sys.stderr.write("\n")
  vals = []
  sys.stderr.write("Traversing annotation file\n")
  for gene in genes:
    for tx in genes[gene]['transcripts']:
      v = genes[gene]['transcripts'][tx]
      exons = txs[tx]
      if len(exons) < args.min_exons: continue
      for i in range(0,v):
        vals += exons[:]
  sys.stderr.write("Generating coverage file "+str(len(vals))+"\n")
  of = sys.stdout
  if args.output:
    if args.output[-3:]=='.gz':
      of = gzip.open(args.output,'w')
    else:
      of = open(args.output,'w')
  covs = ranges_to_coverage(vals)
  for v in covs:
    of.write(v.chr+"\t"+str(v.start-1)+"\t"+str(v.end)+"\t"+str(v.get_payload())+"\n")
  #    of.write(tx+"\t"+gene+"\t"+str(genes[gene]['transcripts'][tx])+"\t"+str(genes[gene]['cnt'])+"\n")
  of.close()
if __name__=="__main__":
  main()
