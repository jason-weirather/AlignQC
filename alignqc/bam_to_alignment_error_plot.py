#!/usr/bin/env python
"""Calculate non-context based errors"""

import argparse, sys, os, time, gc, gzip
from shutil import rmtree
from tempfile import mkdtemp, gettempdir
from subprocess import call

from seqtools.format.sam.bam.files import BAMFile
from seqtools.errors import ErrorProfileFactory
from seqtools.format.fasta import FASTAData
from seqtools.format.sam.bam.bamindex import BAMIndexRandomAccessPrimary as BIRAP

def main(args):

  sys.stderr.write("Reading our reference Fasta\n")
  if args.reference[-3:] == '.gz':
     ref = FASTAData(gzip.open(args.reference).read())
  else:
     ref = FASTAData(open(args.reference).read())
  sys.stderr.write("Finished reading our reference Fasta\n")
  bf = BAMFile(args.input,BAMFile.Options(reference=ref))
  bind = None
  epf = ErrorProfileFactory()
  if args.random:
    sys.stderr.write("Reading index\n")
    if args.input_index:
      bind = BIRAP(index_file=args.input_index,alignment_file=args.input)
    else:
      bind = BIRAP(index_file=args.input+'.bgi',alignment_file=args.input)
    z = 0
    while True:
      coord = bind.get_random_coord()
      if not coord: continue
      e = bf.fetch_by_coord(coord)
      if not e.is_aligned(): continue
      epf.add_alignment(e)
      z+=1
      #print z
      if z %100==1:
        con = epf.get_alignment_errors().alignment_length
        if args.max_length <= con: break
        sys.stderr.write(str(con)+"/"+str(args.max_length)+" bases from "+str(z)+" alignments\r")
    sys.stderr.write("\n")
  else:
    z = 0
    for e in bf:
      if e.is_aligned():
        epf.add_alignment(e)
        z+=1
        #print z
        if z %100==1:
          con = epf.get_alignment_errors().alignment_length
          if args.max_length <= con: break
          sys.stderr.write(str(con)+"/"+str(args.max_length)+" bases from "+str(z)+" alignments\r")
    sys.stderr.write("\n")
  of = open(args.tempdir+'/report.txt','w')
  of.write(epf.get_alignment_errors().get_report())
  of.close()

  for ofile in args.output:
    cmd = [args.rscript_path,
           os.path.dirname(os.path.realpath(__file__))+'/plot_alignment_errors.r',
           args.tempdir+'/report.txt',ofile]
    if args.scale:
      cmd += [str(x) for x in args.scale]
    sys.stderr.write(" ".join(cmd)+"\n")
    call(cmd)


  if args.output_raw:
    of = open(args.output_raw,'w')
    with open(args.tempdir+"/report.txt") as inf:
      for line in inf:
        of.write(line)
    of.close()
  if args.output_stats:
    of = open(args.output_stats,'w')
    of.write(epf.get_alignment_errors().get_stats())
    of.close()
  sys.stderr.write("finished\n")
  if bind:
    bind.destroy()
  bf = None
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
  parser.add_argument('input',help="BAMFILE input")
  parser.add_argument('--input_index',help="BAMFILE index")
  parser.add_argument('-r','--reference',help="Fasta reference file",required=True)
  parser.add_argument('--scale',type=float,nargs=6,help="<ins_min> <ins_max> <mismatch_min> <mismatch_max> <del_min> <del_max>")
  parser.add_argument('-o','--output',nargs='+',help="OUTPUTFILE for pdf plot",required=True)
  parser.add_argument('--output_stats',help="OUTPUTFILE for a stats report")
  parser.add_argument('--output_raw',help="OUTPUTFILE for the raw data")
  parser.add_argument('--random',action='store_true',help="randomly select alignments, requires indexed file")
  parser.add_argument('--max_length',type=int,default=100000,help="maximum number of alignment bases to use")
  # Temporary working directory step 1 of 3 - Definition
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")
  parser.add_argument('--rscript_path',default='Rscript',help="Rscript path")
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
