#!/usr/bin/env python
import argparse, sys, os, gzip
from shutil import rmtree
from multiprocessing import cpu_count
from tempfile import mkdtemp, gettempdir

def main(args):
  chrcovs = {}
  total = {}
  inf = gzip.open(args.input)
  for line in inf:
      f = line.rstrip().split("\t")
      chr = f[0]
      start = int(f[1])
      finish = int(f[2])
      depth = int(f[3])
      if chr not in chrcovs: chrcovs[chr] = {}
      if depth not in chrcovs[chr]: chrcovs[chr][depth] = 0
      chrcovs[chr][depth] += finish-start
      if depth not in total:  total[depth] = 0
      total[depth] += finish-start
  inf.close()
  chrlens = {}
  with open(args.reflens) as inf:
    for line in inf:
      f = line.rstrip().split("\t")
      chrlens[f[0]]=int(f[1])
  total_len = sum(chrlens.values())
  cov_len = sum(total.values())
  print total_len
  print cov_len
  depths = sorted(total.keys())
  #bases = total_len-cov_len
  prev = total_len-cov_len
  oflpt = gzip.open(args.output+"/line_plot_table.txt.gz",'w')
  for d in depths:
    oflpt.write(str(d)+"\t"+str(prev+1)+"\t"+str(total_len)+"\n")
    oflpt.write(str(d)+"\t"+str(prev+total[d])+"\t"+str(total_len)+"\n")
    prev = prev+total[d]
  oflpt.close()
  oftdt = gzip.open(args.output+"/total_distro_table.txt.gz",'w')
  for d in depths:
    oftdt.write(str(d)+"\t"+str(total[d])+"\t"+str(cov_len)+"\t"+str(total_len)+"\n")
  oftdt.close()
  ofcdt = gzip.open(args.output+"/chr_distro_table.txt.gz",'w')
  for chr in sorted(chrcovs.keys()):
    covered_bases = sum(chrcovs[chr].values())
    for depth in sorted(chrcovs[chr].keys()):
      ofcdt.write(chr + "\t" + str(depth)+"\t"+str(chrcovs[chr][depth])+"\t"+str(covered_bases)+"\t"+str(chrlens[chr])+"\n")
  ofcdt.close()

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="depth INPUT FILE")
  parser.add_argument('reflens',help="reflens INPUT FILE")
  parser.add_argument('-o','--output',help="OUTPUT directory")
  args = parser.parse_args()
  return args

def external_cmd(cmd):
  cache_argv = sys.argv
  sys.argv = cmd
  args = do_inputs()
  main(args)
  sys.argv = cache_argv

if __name__=="__main__":
  #do our inputs
  args = do_inputs()
  main(args)
