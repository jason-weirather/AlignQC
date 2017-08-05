#!/usr/bin/env python
import argparse, sys, os, re, gzip
from shutil import rmtree, copyfile
from multiprocessing import cpu_count
from tempfile import mkdtemp, gettempdir

import dump
from seqtools.statistics import average, N50, median, standard_deviation


g_version = None

def main(args):
  #do our inputs
  args.output = args.output.rstrip('/')
  if not os.path.exists(args.output):
    os.makedirs(args.output) 
  
  all_files = []
  z = 0
  sys.stderr.write("Collecting information regarding available inputs\n")
  for infile in args.xhtml_inputs:
    z += 1
    outfile = args.tempdir+'/'+str(z)+'.list'
    cmd = ['dump.py',infile,'-l','-o',outfile]
    dump.external_cmd(cmd)
    res = {'fname':infile,'xlist':set(),'index':z}
    with open(outfile) as inf:
      for line in inf:
        line= line.rstrip()
        primary = re.match('(\S+)',line).group(1)
        res['xlist'].add(primary)
        m = re.match('\S+\s+\[(.+)\]',line)
        if m:
          alts = m.group(1).split(' ')
          for alt in alts: res['xlist'].add(alt)
    all_files.append(res)
  sys.stderr.write("Extracting data for comparison\n")
  for file in all_files:
    sys.stderr.write("  Extracting for "+str(file['fname'])+" "+str(file['index'])+"\n")
    if 'alignment_stats.txt' in file['xlist']:
      ofile = args.tempdir+'/'+str(file['index'])+'.alignment_stats.txt'
      cmd = ['dump.py',file['fname'],'-e','alignment_stats.txt','-o',ofile]
      dump.external_cmd(cmd)
      file['alignment_stats'] = {}
      with open(ofile) as inf:
        for line in inf:
          f = line.rstrip().split("\t")
          file['alignment_stats'][f[0]] = int(f[1])
    if 'error_stats.txt' in file['xlist']:
      ofile = args.tempdir+'/'+str(file['index'])+'.error_stats.txt'
      cmd = ['dump.py',file['fname'],'-e','error_stats.txt','-o',ofile]
      dump.external_cmd(cmd)
      file['error_stats'] = {}
      with open(ofile) as inf:
        for line in inf:
          f = line.rstrip().split("\t")
          file['error_stats'][f[0]] = f[1]
    if 'lengths.txt.gz' in file['xlist']:
      ofile = args.tempdir+'/'+str(file['index'])+'.lengths.txt.gz'
      cmd = ['dump.py',file['fname'],'-e','lengths.txt.gz','-o',ofile]
      dump.external_cmd(cmd)
      file['lengths'] = {'average':'','median':'','N50':'','stddev':'','average_aligned':'','median_aligned':'','N50_aligned':'','stddev_aligned':''}
      inf = gzip.open(ofile)
      lengths = []
      for line in inf:
        f = line.rstrip().split("\t")
        lengths.append([int(f[3]),int(f[4])])
      if len(lengths) > 0:
        file['lengths']['average']=average([x[1] for x in lengths])
        file['lengths']['median']=median([x[1] for x in lengths])
        file['lengths']['N50']=N50([x[1] for x in lengths])
      if len([x[0] for x in lengths if x[0] != 0]) > 0:
        file['lengths']['average_aligned']=average([x[0] for x in lengths if x[0] != 0])
        file['lengths']['median_aligned']=median([x[0] for x in lengths if x[0] != 0])
        file['lengths']['N50_aligned']=N50([x[0] for x in lengths if x[0] != 0])
      if len(lengths) > 2:
        file['lengths']['stddev']=standard_deviation([x[1] for x in lengths])
      if len([x[0] for x in lengths if x[0] != 0]) > 2:
        file['lengths']['stddev_aligned']=standard_deviation([x[0] for x in lengths if x[0] != 0])
      
  # Now we can output table
  ofname = args.tempdir+'/stats_table.txt'
  of = open(ofname,'w')
  header =  'File'+"\t"
  #Basic
  header += 'Reads'+"\t"
  header += 'Avg_Read_Length'+"\t"
  header += 'Median_Read_Length'+"\t"
  header += 'N50_Read_Length'+"\t"
  header += 'Stddev_Read_Length'+"\t"
  header += 'Aligned_Read_Count'+"\t"
  header += 'Aligned_Reads'+"\t"
  header += 'Avg_Aligned_Length'+"\t"
  header += 'Median_Aligned_Length'+"\t"
  header += 'N50_Aligned_Length'+"\t"
  header += 'Stddev_Aligned_Length'+"\t"
  header += 'Chimeric_Total_Reads'+"\t"
  header += 'Chimeric_Trans_Reads'+"\t"
  header += 'Chimeric_Self_Reads'+"\t"
  header += 'Bases'+"\t"
  header += 'Bases_Aligned'+"\t"
  # Error
  header += 'Error_Rate'+"\t"
  header += 'Mismatches'+"\t"
  header += 'Deletions_Total'+"\t"
  header += 'Deletions_Homopolymer'+"\t"
  header += 'Insertions_Total'+"\t"
  header += 'Insertions_Homopolymer'
  of.write(header+"\n")
  for file in all_files:
    of.write(file['fname']+"\t")
    basic = get_basic(file)
    of.write("\t".join([str(x) for x in basic]))    
    of.write("\t")
    error = get_error(file)
    of.write("\t".join([str(x) for x in error]))    
    of.write("\n")
  of.close()

  copyfile(args.tempdir+'/stats_table.txt',args.output+'/stats_table.txt')
  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

# File is the dict that has the 'fname' and 'alignment_stats' and 'error_stats'
# if there is error stats return the array to output
def get_error(file):
  error = ['' for x in range(0,7)]
  if 'error_stats' not in file: return error
  dat = file['error_stats']
  if dat['ALIGNMENT_BASES'] > 0:
    error[0] = float(dat['ANY_ERROR'])/float(dat['ALIGNMENT_BASES'])
  if dat['ALIGNMENT_BASES'] > 0:
    error[1] = float(dat['MISMATCHES'])/float(dat['ALIGNMENT_BASES'])
  if dat['ALIGNMENT_BASES'] > 0:
    error[2] = float(dat['ANY_DELETION'])/float(dat['ALIGNMENT_BASES'])
  if dat['ANY_DELETION'] > 0:
    error[3] = float(dat['HOMOPOLYMER_DELETION'])/float(dat['ANY_DELETION'])
  if dat['ALIGNMENT_BASES'] > 0:
    error[4] = float(dat['ANY_INSERTION'])/float(dat['ALIGNMENT_BASES'])
  if dat['ANY_INSERTION'] > 0:
    error[5] = float(dat['HOMOPOLYMER_INSERTION'])/float(dat['ANY_INSERTION'])
  return error
# File is the dict that has the 'fname' and 'alignment_stats' and 'error_stats'
# if there is alignment stats return the array to output
def get_basic(file):
  basic = ['' for x in range(0,16)]
  if 'alignment_stats' not in file: return basic
  dat = file['alignment_stats']
  basic[0] = dat['TOTAL_READS']

  basic[1] = file['lengths']['average']
  basic[2] = file['lengths']['median']
  basic[3] = file['lengths']['N50']
  basic[4] = file['lengths']['stddev']

  basic[5] = dat['ALIGNED_READS']
  if dat['TOTAL_READS'] > 0:
    basic[6] = float(dat['ALIGNED_READS'])/float(dat['TOTAL_READS'])

  basic[7] = file['lengths']['average_aligned']
  basic[8] = file['lengths']['median_aligned']
  basic[9] = file['lengths']['N50_aligned']
  basic[10] = file['lengths']['stddev_aligned']  

  if dat['ALIGNED_READS'] > 0:
     basic[11] = float(dat['CHIMERA_ALIGN_READS'])/float(dat['ALIGNED_READS'])
  if dat['CHIMERA_ALIGN_READS'] > 0:
    basic[12] = float(dat['TRANSCHIMERA_ALIGN_READS'])/float(dat['CHIMERA_ALIGN_READS'])
  if dat['CHIMERA_ALIGN_READS'] > 0:
    basic[13] = float(dat['SELFCHIMERA_ALIGN_READS'])/float(dat['CHIMERA_ALIGN_READS'])
  basic[14] = dat['TOTAL_BASES']
  if dat['TOTAL_BASES'] > 0:
    basic[15] = float(dat['ALIGNED_BASES'])/float(dat['TOTAL_BASES'])
  return basic

def external_cmd(cmd,version=None):
  #set version by input
  global g_version
  g_version = version

  cache_argv = sys.argv
  sys.argv = cmd
  args = do_inputs()
  main(args)
  sys.argv = cache_argv


def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('xhtml_inputs',nargs='+',help="xhtml analysis files")
  parser.add_argument('-o','--output',required=True,help="OUTPUT directory")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")
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

if __name__=="__main__":
  args = do_inputs()
  main(args)
