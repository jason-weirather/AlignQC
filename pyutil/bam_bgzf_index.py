#!/usr/bin/python
import argparse, sys, os, gzip
from shutil import rmtree, copy
from multiprocessing import cpu_count, Pool
from tempfile import mkdtemp, gettempdir
from Bio.Format.Sam import BAMFile, SAM, check_flag
from subprocess import Popen, PIPE

def main(args):
  #do our inputs
  ind_path = args.input+'.bgi'
  if args.output: ind_path = args.output
  of = gzip.open(args.tempdir+'/myfile.bgi','w')
  if args.output: ind_path = args.output
  if os.path.isfile(ind_path) and not args.output:
    sys.stderr.write("ERROR index file already there.  Delete it if you want to rebuild a new one.\n")
    sys.exit()
  bf = BAMFile(args.input)

  # Reada the location of each line
  entries = []
  z = 0
  sys.stderr.write("read basics\n")
  for e in bf:
    flag = e.value('flag')
    entries.append([e.get_block_start(),e.get_inner_start(),e.value('qname'),flag])
    z+=1
    if z%1000==0:
      sys.stderr.write(str(z)+"       \r")
  sys.stderr.write("\n")
  # now analyze entries to look for best
  sys.stderr.write("check for best set\n")
  fail_primary = False
  reads = {}
  best_flag_set = False
  for i in range(0,len(entries)):
    name = entries[i][2]
    if name not in reads:
      reads[name] = {}
    type = 'u'
    if check_flag(entries[i][3],64):
      type = 'l'
    elif check_flag(entries[i][3],128):
      type = 'r'
    if not check_flag(entries[i][3],2304):
      if type not in reads[name]: reads[name][type] = 0
      reads[name][type] += 1
      if reads[name][type] > 1:
        fail_primary = True
        break
  one_count = True #one alignment per read?
  for name in reads:
    if fail_primary: break
    for type in reads[name]:
      if reads[name][type] != 1:
        fail_primary = True
        one_count = False
        break
  if fail_primary:
    sys.stderr.write("Failed to find a single primary alignment for each read\n")
  # Reading more detailed information
  results = []
  if args.threads > 1:
    p = Pool(processes=args.threads)
  final = []
  chunksize = 10000
  for i in range(0,len(entries),chunksize):
    subset = [x[0:2] for x in entries[i:i+chunksize]]
    sys.stderr.write(str(i)+'/'+str(len(entries))+"       \r")
    if args.threads > 1:
      results.append(p.apply_async(do_chunk,args=(subset[0],len(subset),args,)))
    else:
      results.append(Queue(do_chunk(subset[0],len(subset),args)))
  sys.stderr.write("\n")
  if args.threads > 1:
    p.close()
    p.join()
  sys.stderr.write("combining results\n")
  for r in results:
    for v in r.get():
      final.append(v)
  best = set()
  if fail_primary:
    sys.stderr.write("Now find best entry for each read\n")
    reads = {}
    for i in range(0,len(entries)):
      name = entries[i][2]
      if name not in reads:
        reads[name] = {}
      type = 'u'
      if check_flag(entries[i][3],64):
        type = 'l'
      elif check_flag(entries[i][3],128):
        type = 'r'
      if type not in reads[name]:  reads[name][type] = {'ind':0,'len':-1}
      if final[i][1] > reads[name][type]['len']:
        reads[name][type]['len'] = final[i][1]
        reads[name][type]['ind'] = i
    for name in reads:
      for type in reads[name]:
        best.add(reads[name][type]['ind'])
  elif one_count: #only one read per alignment
    for i in range(0,len(entries)):
      if final[i][1] > 0: best.add(i)
  # now we can print out results
  for i in range(0,len(entries)):
    flag = entries[i][3]
    if i not in best:
      flag = flag | 2304
    of.write(entries[i][2]+"\t"+final[i][0]+"\t"+str(entries[i][0])+"\t"+str(entries[i][1])+"\t"+str(final[i][1])+"\t"+str(flag)+"\n")
  of.close()
  #bf.write_index(args.tempdir+'/myfile.bgi',verbose=True)
  copy(args.tempdir+'/myfile.bgi',ind_path)  

  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

class Queue:
  def __init__(self,val):
    self.val = val
  def get(self):
    return self.val

def do_chunk(coords,ecount,args):
  bf = BAMFile(args.input,blockStart=coords[0],innerStart=coords[1])
  results = []
  for i in range(0,ecount):
    e = bf.read_entry()
    if e.is_aligned():
      rng = e.get_target_range()
      results.append([rng.get_range_string(),e.get_aligned_bases_count()])
    else:
      results.append(['',0])
  return results

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Generate our .bgi index for a bam file",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="INPUT BAM FILE")
  parser.add_argument('--output','-o',help="Specifiy path to write index")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")
  # Temporary working directory step 1 of 3 - Definition
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")
  group.add_argument('--no_primary_search',action='store_true',help="Dont try to assign a primary flag")
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
  sys.argv = cmd.split()
  args = do_inputs()
  main(args)
  sys.argv = cache_argv

if __name__=="__main__":
  args = do_inputs()
  main(args)
