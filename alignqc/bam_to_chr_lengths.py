"""Pretty simply get the reference sequence lengths from a bam header"""

import sys, argparse, os, inspect

from seqtools.format.sam.bam.files import BAMFile

def main(args):
  
  bf = BAMFile(args.input)
  chrlens = bf.header.sequence_lengths
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
  sys.argv = cmd
  args = do_inputs()
  main(args)
  sys.argv = cache_argv

if __name__=="__main__":
  args = do_inputs()
  main(args)
