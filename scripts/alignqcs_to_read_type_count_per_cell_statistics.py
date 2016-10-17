#!/usr/bin/env python
import sys, argparse, re, os
from subprocess import Popen, PIPE
from Bio.Statistics import average, median, standard_deviation

def main():
  parser = argparse.ArgumentParser(description="Assume standard pacbio ccs and subread read name formats, and ONT name formats where _pass_2D _pass_tem _pass_com _fail_2D _fail_tem or _fail_com have been appended to the name.",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('inputs',nargs='+',help="Specify xhtml files")
  parser.add_argument('--aligned',action='store_true',help="restrict output to aligned reads only")
  args = parser.parse_args()
  
  p_pbccs = re.compile('^m[^\/]+\/\d+\/ccs$')
  p_pbsub = re.compile('^m[^\/]+\/\d+\/\d+_\d+$')
  p_ontpass2D = re.compile('^\S+_pass_2D$')
  p_ontpasstem = re.compile('^\S+_pass_tem$')
  p_ontpasscom = re.compile('^\S+_pass_com$')
  p_ontfail2D = re.compile('^\S+_fail_2D$')
  p_ontfailtem = re.compile('^\S+_fail_tem$')
  p_ontfailcom = re.compile('^\S+_fail_com$')
  results = []
  for fname in args.inputs:
    sys.stderr.write("processing "+fname+"\n")
    c = {
    'pacbio':0,
    'ont':0,
    'ccs':0,
    'sub':0,
    'ontpass2D':0,
    'ontpasstem':0,
    'ontpasscom':0,
    'ontfail2D':0,
    'ontfailtem':0,
    'ontfailcom':0,
    'ontpass':0,
    'ontfail':0,
    'ont2D':0,
    'ont1D':0,
    'other':0 }
    cmd = 'alignqc dump '+fname+' -e lengths.txt'
    p = Popen(cmd.split(),stdout=PIPE)
    for line in p.stdout:
      f = line.rstrip().split("\t")
      bp = int(f[3])
      if args.aligned and bp == 0: continue
      rname = f[0]
      if p_pbccs.match(rname):
        c['ccs']+=1
        c['pacbio']+=1
      elif p_pbsub.match(rname):
        c['sub'] += 1
        c['pacbio']+=1
      elif p_ontpass2D.match(rname):
        c['ontpass2D'] += 1
        c['ont2D'] += 1
        c['ont']+=1
        c['ontpass']+=1
      elif p_ontpasstem.match(rname):
        c['ontpasstem'] += 1
        c['ont1D'] += 1
        c['ont']+=1
        c['ontpass']+=1
      elif p_ontpasscom.match(rname):
        c['ontpasscom'] += 1
        c['ont1D'] += 1
        c['ont']+=1
        c['ontpass']+=1
      elif p_ontfail2D.match(rname):
        c['ontfail2D'] += 1
        c['ont2D'] += 1
        c['ont']+=1
        c['ontfail']+=1
      elif p_ontfailtem.match(rname):
        c['ontfailtem'] += 1
        c['ont1D'] += 1
        c['ont']+=1
        c['ontfail']+=1
      elif p_ontfailcom.match(rname):
        c['ontfailcom'] += 1
        c['ont1D'] += 1
        c['ont']+=1
        c['ontfail']+=1
      else:
        c['other']+=1
    p.communicate()
    results.append(c)
  k = results[0].keys()
  print "feature\tsum\tmedian\taverage\tstandard deviation"
  for feature in sorted(k):
    arr = [x[feature] for x in results]
    print feature+"\t"+str(sum(arr))+"\t"+str(median(arr))+"\t"+str(average(arr))+"\t"+str(standard_deviation(arr))

if __name__=="__main__":
  main()
