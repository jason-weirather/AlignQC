#!/usr/bin/env python
import argparse, sys, os, random, inspect, os, time, gc, gzip
from shutil import rmtree
from multiprocessing import cpu_count
from tempfile import mkdtemp, gettempdir
from subprocess import call

from seqtools.format.sam.bam.files import BAMFile
from seqtools.format.sam.bam.bamindex import BAMIndexRandomAccessPrimary as BIRAP
from seqtools.errors import ErrorProfileFactory
from seqtools.format.fasta import FASTAData

# Take the bam file as an input and produce plots and data file for context errors.

def main(args):
  # make our error profile report
  sys.stderr.write("Reading reference fasta\n")
  if args.reference[-3:] == '.gz':
     ref = FASTAData(gzip.open(args.reference).read())
  else:
     ref = FASTAData(open(args.reference).read())
  sys.stderr.write("Reading index\n")
  epf = ErrorProfileFactory()
  bf = BAMFile(args.input,BAMFile.Options(reference=ref))
  bind = None
  if args.random:
    if args.input_index:
      bind = BIRAP(index_file=args.input_index,alignment_file=args.input)
    else:
      bind = BIRAP(index_file=args.input+'.bgi',alignment_file=args.input)
    z = 0
    strand = 'target'
    if args.query: strand = 'query'
    con = 0
    while True:
      #rname = random.choice(bf.index.get_names())
      #print rname
      #coord = bf.index.get_longest_target_alignment_coords_by_name(rname)
      #print coord
      coord = bind.get_random_coord()
      if not coord: continue
      e = bf.fetch_by_coord(coord)
      if not e.is_aligned(): continue
      epf.add_alignment(e)
      z+=1
      if z%100==1:
        con = epf.get_min_context_count(strand)
      sys.stderr.write(str(z)+" alignments, "+str(con)+" min context coverage\r")
      if args.max_alignments <= z: break
      if args.stopping_point <= con: break

  else:
    z = 0
    strand = 'target'
    if args.query: strand = 'query'
    con = 0
    for e in bf:
      if e.is_aligned():
        epf.add_alignment(e)
        z+=1
        if z%100==1:
          con = epf.get_min_context_count(strand)
        sys.stderr.write(str(z)+" alignments, "+str(con)+" min context coverage\r")
        if args.max_alignments <= z: break
        if args.stopping_point <= con: break
  sys.stderr.write("\n")
  #if bf.index:
  #  bf.index.destroy()
  bf = None
  if bind:
    bind.destroy()
  sys.stderr.write('working with:'+"\n")
  sys.stderr.write(str(z)+" alignments, "+str(con)+" min context coverage"+"\n")
  epf.write_context_error_report(args.tempdir+'/err.txt',strand)

  for ofile in args.output:
    cmd = [args.rscript_path,
           os.path.dirname(os.path.realpath(__file__))+'/plot_base_error_context.r',
           args.tempdir+'/err.txt',ofile]
    if args.scale:
      cmd += [str(x) for x in args.scale]
    sys.stderr.write(" ".join(cmd)+"\n")
    call(cmd)
  sys.stderr.write("finished\n")
  if args.output_raw:
    of = open(args.output_raw,'w')
    with open(args.tempdir+"/err.txt") as inf:
      for line in inf:
        of.write(line)
  epf.close()
  time.sleep(5)
  gc.collect()
  time.sleep(5)
  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="INPUT FILE or '-' for STDIN")
  parser.add_argument('--input_index',help="Index file for bam if other than default location")
  parser.add_argument('-r','--reference',required=True,help="Reference Genome")
  parser.add_argument('-o','--output',nargs='+',required=True,help="OUTPUTFILE(s)")
  parser.add_argument('--output_raw',help="Save the raw data")
  parser.add_argument('--scale',type=float,nargs=6,help="<insertion_min> <insertion_max> <mismatch_min> <mismatch_max> <deletion_min> <deletion_max>")
  parser.add_argument('--max_alignments',type=int,default=10000000000,help="The maximum number of alignments to scan")
  parser.add_argument('--stopping_point',type=int,default=1000,help="Stop after you see this many of each context")
  #parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")
  # Temporary working directory step 1 of 3 - Definition
  group1 = parser.add_mutually_exclusive_group(required=True)
  group1.add_argument('--target',action='store_true',help="Context on the target sequence")
  group1.add_argument('--query',action='store_true',help="Context on the query sequence")
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")
  parser.add_argument('--random',action='store_true',help="Randomly select alignments, requires an indexed bam")
  parser.add_argument('--rscript_path',default='Rscript',help="Path to Rscript")
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
  args = do_inputs()
  main(args)
  sys.argv = cache_argv

if __name__=="__main__":
  #do our inputs
  args = do_inputs()
  main(args)
