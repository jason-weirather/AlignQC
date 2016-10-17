#!/usr/bin/env python
import sys, argparse, os, inspect

#bring in the folder to the path for our utilities
pythonfolder_loc = "../pylib"
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe() ))[0],pythonfolder_loc)))
if cmd_subfolder not in sys.path:
  sys.path.insert(0,cmd_subfolder)

from Bio.Format.Sam import BAMFile

def main(args):
  
  bf = BAMFile(args.input)
  chrlens = bf.get_header().get_sequence_lengths()
  of_chrlens = open(args.output,'w')
  for qname in sorted(chrlens.keys()):
    of_chrlens.write(qname+"\t"+str(chrlens[qname])+"\n")
  of_chrlens.close()

def do_inputs():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use BAM file")
  parser.add_argument('-o','--output',help="output file")
  args = parser.parse_args()
  return args

def external_cmd(cmd):
  cache_argv=sys.argv
  sys.argv = cmd.split()
  args = do_inputs()
  main(args)
  sys.argv = cache_argv

if __name__=="__main__":
  args = do_inputs()
  main(args)
