"""It is necessary to traverse the bam file sort some data by read name"""

import argparse, sys, os, gzip, pickle, zlib, base64
from shutil import rmtree, copy
from multiprocessing import cpu_count, Pool, Lock
from tempfile import mkdtemp, gettempdir
from subprocess import Popen, PIPE

from seqtools.format.sam.bam.files import BAMFile
from seqtools.range import GenomicRange

import seqtools.cli.utilities.bam_bgzf_index as bam_bgzf_index


## The purpose of this script is to read through a bam alignment and record as much information as possible from it.  ##
## The bam should be indexed ahead of time in our index format.

gfinished = None
gtotal = None
glock = Lock()

g_count = 0
g_sortpipe = None

def do_chunk(ilines,infile,args):
  """Takes in a the lines from the index file to work on in array form,
     and the bam file name, and the arguments

     returns a list of the necessary data for chimera detection ready for sorting
  """
  ilines = [x.rstrip().split("\t") for x in ilines]
  coord = [int(x) for x in ilines[0][2:4]]
  bf = BAMFile(infile,BAMFile.Options(blockStart=coord[0],innerStart=coord[1]))
  results = []
  for i in range(0,len(ilines)):
    flag = int(ilines[i][5])
    e = bf.read_entry()
    #if not e: break
    value = None
    if e.is_aligned():
      tx = e.get_target_transcript(args.minimum_intron_size)
      value =  {'qrng':e.actual_original_query_range.get_range_string(),'tx':tx.get_gpd_line(),'flag':flag,'qlen':e.original_query_sequence_length,'aligned_bases':e.get_aligned_bases_count()}
      results.append(e.entries.qname+"\t"+base64.b64encode(
                                      zlib.compress(
                                       pickle.dumps(value))))
      #results.append([e.value('qname'),zlib.compress(pickle.dumps(value))])
    else:
      value =  {'qrng':'','tx':'','flag':flag,'qlen':e.original_query_sequence_length,'aligned_bases':0}
      results.append(e.entries.qname+"\t"+base64.b64encode(
                                      zlib.compress(
                                       pickle.dumps(value))))
      #results.append([e.value('qname'),zlib.compress(pickle.dumps(value))])
  return results

def process_chunk(res):
  global glock
  glock.acquire()
  #global g_preordered
  global g_sortpipe
  global g_count
  g_count += len(res)
  for val in res:
    g_sortpipe.stdin.write(val+"\n")
  sys.stderr.write(str(g_count)+"     \r")
  glock.release()

def main(args):
  bind_path = args.input+'.bgi'
  if not os.path.isfile(bind_path):
    bind_path = args.tempdir+'/myindex.bgi'
    cmd = ["bam_bgzf_index.py",args.input,"-o",bind_path,"--threads",str(args.threads)]
    bam_bgzf_index.external_cmd(cmd)
    #call(cmd.split())

  #parallel_thread = ''
  #if args.threads > 1: parallel_thread = ' --parallel='+str(args.threads)+' '
  #cmd1 = 'sort '+parallel_thread+' -k1,1 -T '+args.tempdir+'/'
  if args.threads > 1:
    cmd1 = ['sort','-k1,1','-T',args.tempdir+'/',
            '--parallel='+str(args.threads)]
  else:
    cmd1 = ['sort','-k1,1','-T',args.tempdir+'/']
  cmd2 = 'gzip'
  global g_sortpipe
  global g_count
  g_count = 0
  of = open(args.output,'wb')
  if os.name != 'nt':
     gzippipe = Popen(cmd2.split(),stdout=of,stdin=PIPE,close_fds=True)
     g_sortpipe = Popen(cmd1,stdout=gzippipe.stdin,stdin=PIPE,close_fds=True)
  else:
     sys.stderr.write("WARNING: Windows OS detected. operating in single thread mode.\n")
     if args.threads > 1: raise ValueError('Error. --threads must be 1 for windows operation')
     gzippipe = Popen(cmd2,stdout=of,stdin=PIPE, shell=True)
     g_sortpipe = Popen(cmd1,stdout=gzippipe.stdin,stdin=PIPE, shell=True)
  inf = gzip.open(bind_path)
  chunksize = args.chunk_size
  buffer = []
  if args.threads > 1:
    p = Pool(processes=args.threads)
  for line in inf:
    buffer.append(line)
    if len(buffer)>=chunksize:
      if args.threads > 1:
        p.apply_async(do_chunk,args=(buffer[:],args.input,args),callback=process_chunk)
      else:
        r = do_chunk(buffer[:],args.input,args)
        process_chunk(r)
      buffer = []
  if len(buffer) > 0:
    if args.threads > 1:
        p.apply_async(do_chunk,args=(buffer[:],args.input,args),callback=process_chunk)
    else:
      r= do_chunk(buffer[:],args.input,args)
      process_chunk(r)
  if args.threads > 1:
    p.close()
    p.join()
  inf.close()
  sys.stderr.write("\n")
  g_sortpipe.communicate()
  gzippipe.communicate()
  of.close()


def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="BAMFILE input")
  parser.add_argument('-o','--output',help="gzipped output",required=True)
  parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")
  parser.add_argument('--minimum_intron_size',type=int,default=68)
  parser.add_argument('--chunk_size',type=int,default=10000,help="number of alignments to process at a time")

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
  #need to save arguments
  cache_argv = sys.argv
  sys.argv = cmd
  args = do_inputs()
  main(args)
  #need to set the arguments back to what they were
  sys.argv = cache_argv
  return

if __name__=="__main__":
  #do our inputs
  args = do_inputs()
  main(args)
