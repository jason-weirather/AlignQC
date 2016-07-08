#!/usr/bin/python
import argparse, os, inspect, sys
from subprocess import Popen, PIPE
from tempfile import mkdtemp, gettempdir
from multiprocessing import cpu_count
from shutil import rmtree


import prepare_all_data
import create_html

g_version = None

def main(args):

  if not args.output and not args.portable_output and not args.output_folder:
    sys.stderr.write("ERROR: must specify some kind of output\n")
    sys.exit()

  ## Check and see if directory for outputs exists
  if args.output_folder:
    if os.path.isdir(args.output_folder):
      sys.stderr.write("ERROR: output directory already exists.  Remove it to to use this location\n")
      sys.exit()


  global g_version
  #Make sure rscript is installed
  try:
    cmd = 'Rscript --version'
    prscript = Popen(cmd.split(),stdout=PIPE,stderr=PIPE)
    rline = prscript.communicate()
    sys.stderr.write("Using Rscript version:\n")
    sys.stderr.write(rline[1].rstrip()+"\n")
  except:
    sys.stderr.write("ERROR: Rscript not installed\n")
    sys.exit()

  if args.no_reference:
    sys.stderr.write("WARNING: No reference specified.  Will be unable to report error profile\n")
  if args.no_annotation:
    sys.stderr.write("WARNING: No annotation specified.  Will be unable to report feature specific outputs\n")

  prepare_all_data.external(args)
  create_html.external(args,version=g_version)

  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

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
  parser=argparse.ArgumentParser(description="Create an output report",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  label1 = parser.add_argument_group(title="Input parameters",description="Required BAM file.  If reference or annotation is not set, use --no_reference or --no_annotation respectively to continue.")
  label1.add_argument('input',help="INPUT BAM file")

  group1 = label1.add_mutually_exclusive_group(required=True)
  group1.add_argument('-r','--reference',help="Reference Fasta")
  group1.add_argument('--no_reference',action='store_true',help="No Reference Fasta")
  group2 = label1.add_mutually_exclusive_group(required=True)
  group2.add_argument('-a','--annotation',help="Reference annotation genePred")
  group2.add_argument('--no_annotation',action='store_true',help="No annotation is available")

  # output options
  label2 = parser.add_argument_group(title="Output parameters",description="At least one output parameter must be set")
  label2.add_argument('-o','--output',help="OUTPUT xhtml with data")
  label2.add_argument('--portable_output',help="OUTPUT file in a small xhtml format")
  label2.add_argument('--output_folder',help="OUTPUT folder of all data")

  label3 = parser.add_argument_group(title="Performance parameters")
  label3.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")

  # Temporary working directory step 1 of 3 - Definition
  label4 = parser.add_argument_group(title="Temporary folder parameters")
  group = label4.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")

  ### Parameters for alignment plots
  label5 = parser.add_argument_group(title="Alignment plot parameters")
  label5.add_argument('--min_intron_size',type=int,default=68,help="minimum intron size when smoothing indels")
  label5.add_argument('--min_aligned_bases',type=int,default=50,help="for analysizing alignment, minimum bases to consider")
  label5.add_argument('--max_query_overlap',type=int,default=10,help="for testing gapped alignment advantage")
  label5.add_argument('--max_target_overlap',type=int,default=10,help="for testing gapped alignment advantage")
  label5.add_argument('--max_query_gap',type=int,help="for testing gapped alignment advantge")
  label5.add_argument('--max_target_gap',type=int,default=500000,help="for testing gapped alignment advantage")
  label5.add_argument('--required_fractional_improvement',type=float,default=0.2,help="require gapped alignment to be this much better (in alignment length) than single alignment to consider it.")
  
  ### Parameters for locus analysis
  label6 = parser.add_argument_group(title="Locus parameters",description="Optionally produce plots and data regarding clusters of sequences")
  label6.add_argument('--do_loci',action='store_true',help="this analysis is time consuming at the moment\n")
  label6.add_argument('--min_depth',type=float,default=1.5,help="require this or more read depth to consider locus")
  label6.add_argument('--min_coverage_at_depth',type=float,default=0.8,help="require at leas this much of the read be covered at min_depth")
  label6.add_argument('--min_exon_count',type=int,default=2,help="Require at least this many exons in a read to consider assignment to a locus")
  label6.add_argument('--locus_downsample',type=int,default=100,help="Limit how deep to search loci\n")

  ### Params for alignment error plot
  label7 = parser.add_argument_group(title="Alignment error parameters")
  label7.add_argument('--alignment_error_scale',nargs=6,type=float,help="<ins_min> <ins_max> <mismatch_min> <mismatch_max> <del_min> <del_max>")
  label7.add_argument('--alignment_error_max_length',type=int,default=1000000,help="The maximum number of alignment bases to calculate error from")
  
  ### Params for context error plot
  label8 = parser.add_argument_group(title="Context error parameters")
  label8.add_argument('--context_error_scale',nargs=6,type=float,help="<ins_min> <ins_max> <mismatch_min> <mismatch_max> <del_min> <del_max>")
  label8.add_argument('--context_error_stopping_point',type=int,default=10000,help="Sample at least this number of each context")

  ## Params for rarefraction plots
  label9 = parser.add_argument_group(title="Rarefraction plot parameters")
  label9.add_argument('--samples_per_xval',type=int,default=500)

  args = parser.parse_args()
  setup_tempdir(args)
  return args

if __name__=='__main__':
  do_inputs()
  main(args)
