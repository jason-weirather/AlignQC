#!/usr/bin/python
import argparse, sys, os, inspect
from shutil import rmtree
from multiprocessing import cpu_count
from tempfile import mkdtemp, gettempdir

from Bio.Format.Bed import Bed12

#bring in the folder to the path for our utilities
pythonfolder_loc = "../pyutil"
#pythonfolder_loc = "../../Au-public/iron/utilities"
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe() ))[0],pythonfolder_loc)))
if cmd_subfolder not in sys.path:
  sys.path.insert(0,cmd_subfolder)

import bam_bgzf_index

g_version = None

def main(args):
  ofname = args.input[0]+'.bgi'
  if args.output:
    ofname = args.output+'.bgi'
  if len(args.input) > 1 and args.output: 
    sys.stderr.write("ERROR cannot specify specific output if more than one input.  Must use default output (not set) for multiple inputs.  will append .bgi")
  if os.path.exists(ofname):
    sys.stderr.write("ERROR: "+ofname+"\nalready exists\n")
    sys.exit()
  for iname in args.input:
    if len(args.input) > 0:
      ofname = iname+'.bgi'
    cmd = 'bam_bgzf_index.py --threads '+str(args.threads)+' -o '+ofname+' --specific_tempdir '+args.tempdir+' '+iname
    sys.stderr.write('Indexing '+iname+"\n")
    bam_bgzf_index.external_cmd(cmd)
  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

def external_cmd(cmd,version=None):
  #set version by input
  global g_version
  g_version = version

  cache_argv = sys.argv
  sys.argv = cmd.split()
  args = do_inputs()
  main(args)
  sys.argv = cache_argv

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Index BAM File(s).  Creates an index used by alignqc.  Its gzip format, and contains some useful information in tsv\n<qname> <target range> <bgzf file block start> <bgzf inner block start> <aligned base count> <flag>\nFlag may be different than BAM line because it now only contains one primary alignment for each read name.",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('-o','--output',help="Specify the output location (only can do one input)")
  parser.add_argument('input',nargs='+',help="Input BAM(s)")
  parser.add_argument('--threads',default=cpu_count(),help="number of processes")
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

if __name__=="__main__":
  main()
