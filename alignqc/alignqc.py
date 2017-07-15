#!/usr/bin/env python
"""The front end of alignqc"""
import argparse, sys

# Import our main launchers for each mode
import analyze
import dump
import compare
from seqtools.cli.utilities.bam_bgzf_index import external_cmd as index_bam

version = '1.3'

def main():
  operable_argv = [sys.argv[0]]+sys.argv[2:]
  sys.argv = sys.argv[:2]
  #do our inputs
  args = do_inputs()
  if args.mode == 'analyze':
    analyze.external_cmd(" ".join(operable_argv),version=version)
  elif args.mode == 'dump':
    dump.external_cmd(" ".join(operable_argv),version=version)
  elif args.mode == 'compare':
    compare.external_cmd(" ".join(operable_argv),version=version)
  else:
    sys.stderr.write("Run mode not yet implemented\n")

def do_inputs():
  # Setup command line inputs
  global version
  parser=argparse.ArgumentParser(description="Version "+str(version)+"\nReview reports about alignments.",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('mode',choices=['analyze','dump','compare'],help="MODE of program to run")
  args = parser.parse_args()
  return args


if __name__=="__main__":
  main()
