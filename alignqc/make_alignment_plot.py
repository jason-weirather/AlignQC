#!/usr/bin/env python
import argparse, sys, os, gzip, re
from shutil import rmtree, copy
from multiprocessing import cpu_count
from tempfile import mkdtemp, gettempdir
from subprocess import call

def main(args):
  udir = os.path.dirname(os.path.realpath(__file__))
  #sys.stderr.write("Making text report\n")

  sys.stderr.write("making plot\n")
  for ofile in args.output:
    cmd = [args.rscript_path,udir+'/plot_gapped_alignment_statistics.r',
           args.input,ofile]
    sys.stderr.write(" ".join(cmd)+"\n")
    call(cmd)

  if args.output_stats:
    do_stats(args)
  sys.stderr.write("Finished.\n")
  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

def do_stats(args):
  total_reads = 0
  unaligned_reads = 0
  aligned_reads = 0
  single_align_reads = 0
  gapped_align_reads = 0
  chimera_align_reads = 0
  selfchimera_align_reads = 0
  transchimera_align_reads = 0
  total_bases = 0
  unaligned_bases = 0
  aligned_bases = 0
  single_align_bases = 0
  gapped_align_bases = 0
  chimera_align_bases = 0
  selfchimera_align_bases = 0
  transchimera_align_bases = 0

  inf = None
  if re.search('\.gz',args.input):
    inf = gzip.open(args.input)
  else:
    inf = open(args.input)

  for line in inf:
      (name, type, single, both, rlen) = line.rstrip().split("\t")
      single = int(single)
      both = int(both)
      rlen = int(rlen)
      total_reads += 1
      if type=="unaligned": unaligned_reads += 1
      else: aligned_reads += 1
      if type=="original": single_align_reads +=1
      if type=="gapped": 
        gapped_align_reads += 1
        gapped_align_bases += both-single
      if type=="chimera": 
        transchimera_align_reads += 1
        transchimera_align_bases += both-single
      if type=="self-chimera" or type=="self-chimera-atypical": 
        selfchimera_align_reads += 1
        selfchimera_align_bases += both-single
      if re.search('chimera',type): 
        chimera_align_reads +=1
        chimera_align_bases += both-single
      if type!="unaligned":
        total_bases += rlen
        unaligned_bases += (rlen-both)
        aligned_bases += both
        single_align_bases += single
  of = open(args.output_stats,'w')
  of.write("TOTAL_READS\t"+str(total_reads)+"\n")
  of.write("UNALIGNED_READS\t"+str(unaligned_reads)+"\n")
  of.write("ALIGNED_READS\t"+str(aligned_reads)+"\n")
  of.write("SINGLE_ALIGN_READS\t"+str(single_align_reads)+"\n")
  of.write("GAPPED_ALIGN_READS\t"+str(gapped_align_reads)+"\n")
  of.write("CHIMERA_ALIGN_READS\t"+str(chimera_align_reads)+"\n")
  of.write("TRANSCHIMERA_ALIGN_READS\t"+str(transchimera_align_reads)+"\n")
  of.write("SELFCHIMERA_ALIGN_READS\t"+str(selfchimera_align_reads)+"\n")
  of.write("TOTAL_BASES\t"+str(total_bases)+"\n")
  of.write("UNALIGNED_BASES\t"+str(unaligned_bases)+"\n")
  of.write("ALIGNED_BASES\t"+str(aligned_bases)+"\n")
  of.write("SINGLE_ALIGN_BASES\t"+str(single_align_bases)+"\n")
  of.write("GAPPED_ALIGN_BASES\t"+str(gapped_align_bases)+"\n")
  of.write("CHIMERA_ALIGN_BASES\t"+str(chimera_align_bases)+"\n")
  of.write("TRANSCHIMERA_ALIGN_BASES\t"+str(transchimera_align_bases)+"\n")
  of.write("SELFCHIMERA_ALIGN_BASES\t"+str(selfchimera_align_bases)+"\n")

  of.close()
  inf.close()
def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="INPUT lengths.txt file")
  parser.add_argument('-o','--output',nargs='+',help="OUTPUT FILE can put multiple")
  parser.add_argument('--output_stats',help="Save some summary statistics")
  parser.add_argument('--rscript_path',default='Rscript',help="Path of Rscript")

  # Temporary working directory step 1 of 3 - Definition
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")
  args = parser.parse_args()
  # Temporary working directory step 2 of 3 - Creation
  setup_tempdir(args)
  return args

def setup_tempdir(args):
  if args.specific_tempdir:
    if not os.path.exists(args.specific_tempdir):
      os.makedirs(args.specific_tempdir.rstrip('/'))
    args.tempdir = args.specific_tempdir.rstrip('/')
    if not os.path.exists(args.specific_tempdir.rstrip('/')):
      sys.stderr.write("ERROR: Problem creating temporary directory\n")
      sys.exit()
  else:
    args.tempdir = mkdtemp(prefix="weirathe.",dir=args.tempdir.rstrip('/'))
    if not os.path.exists(args.tempdir.rstrip('/')):
      sys.stderr.write("ERROR: Problem creating temporary directory\n")
      sys.exit()
  if not os.path.exists(args.tempdir):
    sys.stderr.write("ERROR: Problem creating temporary directory\n")
    sys.exit()
  return 

def external_cmd(cmd):
  cache_argv = sys.argv
  sys.argv = cmd
  #do our inputs
  args = do_inputs()
  main(args)
  sys.argv = cache_argv

if __name__=="__main__":
  #do our inputs
  args = do_inputs()
  main(args)
