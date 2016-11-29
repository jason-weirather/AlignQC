#!/usr/bin/python
import sys, argparse, gzip
from Bio.Format.Fasta import FastaData
import random

g_version = None

def main(args):
  random.seed(args.seed)
  sum = 0
  if args.reference_genome:
    ref = FastaData(open(args.reference_genome).read())
    for name in ref.keys():
      sum += len(ref[name])
  else:
    with open(args.reference_lengths) as inf:
      for line in inf:
        f = line.rstrip().split("\t")
        sum += int(f[1])
  c = args.minimum_coverage
  z = 0
  values = {}
  while c < sum:
    z += 1
    values[c] = z
    c = c*5
    if c >= sum: break
    z += 1
    values[c] = z
    c = c*2
  z +=1
  values[sum] = z
  for c in sorted(values.keys()):
    values[c] = z-values[c]+1
  ### Now values contains the stratified coverage values
  if args.output_key:
    of = open(args.output_key,'w')
    of.write("bp_size\tstrata_label\n")
    for c in sorted(values.keys()):
      of.write(str(c)+"\t"+str(values[c])+"\n")
    of.close()
  inf = sys.stdin
  if args.input != '-': 
    if args.input[-3:]=='.gz': inf = gzip.open(args.input)
    else: inf = open(args.input)

  of = sys.stdout
  if args.output:
    if args.output[-3:]=='.gz': of = gzip.open(args.output,'w')
    else: of = open(args.output,'w')
  depths = {}
  vals = []
  z = 0
  for line in inf:
    z += 1
    if z % 100000 == 0: sys.stderr.write(str(z)+"    bed entries read   \r")
    f = line.rstrip().split("\t")
    addition = 0
    if not args.dont_make_unique: addition = +args.unique_scale*random.random()
    vals.append([f[0],int(f[1]),int(f[2]),float(f[3])+addition])
  z = 0
  sys.stderr.write("\n")
  for f in vals:
    z += 1
    if z % 100000 == 0: sys.stderr.write(str(z)+"    bed entries read   \r")
    #keep track of the number of bases at each depth
    depth = f[3]
    cov = f[2]-f[1]
    if depth not in depths:  depths[depth] = 0
    depths[depth] += cov
    #vals.append([f[0],int(f[1]),int(f[2]),depth])
  sys.stderr.write("\n")
  #total_bases = sum(depths.values())
  #thresh = {}
  #for strata in stratas:
  #  pos = 0
  #  cur = float(i)*float(total_bases)/float(args.strata)
  stratas = sorted(values.keys())
  pos = 0
  depth_strata = {}
  for d in reversed(sorted(depths.keys())):
    pos += depths[d]
    while stratas[0] < pos:
      stratas.pop(0)
    depth_strata[d] = values[stratas[0]]
    #print str(d)+"\t"+str(values[stratas[0]])
    #if float(pos) > cur:
    #  thresh[d] = [pos,i]
    #  break
  vals[0][3] = depth_strata[vals[0][3]]
  buffer = vals[0]
  for val in vals[1:]:
    val[3] = depth_strata[val[3]]
    if val[1]==buffer[2] and val[3]==buffer[3] and val[0]==buffer[0]:
      #print 'hello'
      buffer[2] = val[2]
      continue
    else:
      of.write(buffer[0]+"\t"+str(buffer[1])+"\t"+str(buffer[2])+"\t"+str(buffer[3])+"\n")
      buffer = val
  of.write(buffer[0]+"\t"+str(buffer[1])+"\t"+str(buffer[2])+"\t"+str(buffer[3])+"\n")
  of.close()


def do_inputs():
  parser = argparse.ArgumentParser(description="convert bed depth from depth to strata",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="bed or - for STDIN")
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument('-r','--reference_genome',help="fasta reference genome")
  group.add_argument('-l','--reference_lengths',help="lenths of reference chromosomes TSV <chr> <length>")
  #parser.add_argument('reference_genome',help="fasta reference genome")
  #parser.add_argument('strata',type=int,help="number of strata to group reads into")
  parser.add_argument('--minimum_coverage',default=10000,type=int,choices=[1,10,100,1000,10000,100000,1000000,10000000,100000000],help="at least this many bases.")
  parser.add_argument('--output_key',help="the key file")
  parser.add_argument('-o','--output',help="output file")
  parser.add_argument('--dont_make_unique',action='store_true',help="do not add small number to make expression unique")
  parser.add_argument('--seed',type=int,default=1,help="int to seed rng")
  parser.add_argument('--unique_scale',type=float,default=0.001,help="add float scaled by this*rng")
  args = parser.parse_args()
  return args

def external_cmd(cmd,version=None):
  #set version by input
  global g_version
  g_version = version

  cache_argv = sys.argv
  sys.argv = cmd.split()
  args = do_inputs()
  main(args)
  sys.argv = cache_argv
  
if __name__=="__main__":
  args = do_inputs()
  main(args)
