#!/usr/bin/env python
import argparse, sys, os, re, base64, zlib, gzip, StringIO
from shutil import rmtree
from multiprocessing import cpu_count
from tempfile import mkdtemp, gettempdir
from xml.etree import ElementTree

from seqtools.format.bed import Bed12

g_version = None

def fixtag(ns, tag, nsmap):
  return '{'+nsmap[ns]+'}'+tag

def main(args):
  of = sys.stdout
  if args.output != '-':
    of = open(args.output,'w')
  names = {}
  tree = ElementTree.parse(args.xhtml_input)
  if args.verbose:
    sys.stderr.write("Traversing xhtml\n")
  for v in [x for x in tree.iter() if x.tag=='{http://www.w3.org/1999/xhtml}a']:
    data = v.attrib
    name = None
    if 'download' in data:
      name = data['download']
    elif 'id' in data:
      name = data['id']
    if not name:
      if args.verbose:
        sys.stderr.write("warning no name for linked data\n")
      continue
    info = data['href']
    m = re.match('data:[^,]+base64,',info)
    if not m: continue
    if args.list:
      names[name] = ''
    m = re.match('data:[^,]+base64,(.*)$',info)
    if not m:
      sys.stderr.write("warning unable to get base64 string")
      continue
    v = base64.b64decode(m.group(1))
    names[name] = v
  if args.verbose:
    sys.stderr.write("Finished traversing xhtml\n")
  if args.list:
    for name in sorted(names.keys()):
      newname = name
      if name[-3:]=='.gz':
        newname = name[:-3]
      if is_ucsc_bed(newname):
        of.write(name+" ["+newname+" "+newname[:-4]+".gpd]\n")        
      else: 
        if name[-3:]=='.gz':
          of.write(name +" ["+newname+"]\n")
        else:
          of.write(newname+"\n")
    return
  shortnames = {}
  for name in names:
    if name[-3:] == '.gz': shortnames[name[:-3]] = name
    else: shortnames[name] = name
  ### if we are still here we must be doing an extract
  exname = args.extract
  out_gzipped = False
  if args.extract[-3:]=='.gz': 
    out_gzipped = True
    exname = args.extract[:-3]
  #handle case of a gpd conversion
  if is_ucsc_gpd(exname):
    oname = exname[:-4]+'.bed.gz'
    if oname not in names:
      sys.stderr.write("ERROR '"+args.extract+"' is not found. Use --list option to see what is available\n")
      sys.exit()
    sio = StringIO.StringIO(zlib.decompress(names[oname],15+32))
    header = sio.readline()
    for v in sio:
      b = Bed12(v)
      print b.get_gpd_line()
    return
  # figure out what the stored name is
  oname = shortnames[exname]
  in_gzipped = False
  if oname[-3:]=='.gz': in_gzipped = True
  if exname not in shortnames:
    sys.stderr.write("ERROR '"+args.extract+"' is not found. Use --list option to see what is available\n")
    sys.exit()
  if in_gzipped and not out_gzipped:
    of.write(zlib.decompress(names[oname],15+32)) #special for gzip format
  else:
    of.write(names[oname])
  #of.write(names[args.extract])
  of.close()
  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

def is_ucsc_bed(newname):
      if newname == 'best.sorted.bed' or \
         newname== 'chimera.bed' or \
         newname== 'gapped.bed' or \
         newname== 'techinical_chimeras.bed' or \
         newname == 'techinical_atypical_chimeras.bed':
        return True
      return False
def is_ucsc_gpd(newname):
      if newname == 'best.sorted.gpd' or \
         newname== 'chimera.gpd' or \
         newname== 'gapped.gpd' or \
         newname== 'technical_chimeras.gpd' or \
         newname == 'technical_atypical_chimeras.gpd':
        return True
      return False
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
  parser=argparse.ArgumentParser(description="Extract data from xhtml output of alignqc through the command line",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('xhtml_input',help="INPUT XHTML FILE")
  parser.add_argument('-o','--output',default='-',help="OUTPUTFILE or STDOUT if not set")
  parser.add_argument('-v','--verbose',action='store_true',help="Show all stderr messages\n")
  group1 = parser.add_mutually_exclusive_group(required=True)
  group1.add_argument('-l','--list',action='store_true',help='show available data')
  group1.add_argument('-e','--extract',help='dump this data')
  #parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")
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
  main()
