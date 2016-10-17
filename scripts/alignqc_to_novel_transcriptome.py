#!/usr/bin/env python
import argparse, sys, os, gzip, inspect
from shutil import rmtree, copy
from multiprocessing import cpu_count
from tempfile import mkdtemp, gettempdir
from subprocess import PIPE, Popen

from Bio.Format.GPD import GPDStream
from Bio.Stream import LocusStream

import classify_reads

#bring in the folder to the path for our utilities
#pythonfolder_loc = "../pyutil"
pythonfolder_loc = "../../Au-public/iron/utilities"
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe() ))[0],pythonfolder_loc)))
if cmd_subfolder not in sys.path:
  sys.path.insert(0,cmd_subfolder)

import gpd_to_nr
import gpd_annotate

def main():
  #do our inputs
  args = do_inputs()


  # first we need to run the classify
  classify_reads.external_cmd('classify_reads.py '+args.input_annot+' '+args.input_gpd+' -o '+args.tempdir+'/classify.txt.gz')

  get_novel_sets(args.tempdir+'/classify.txt.gz',args.input_gpd,args.tempdir+'/novel_isoform_reads.gpd.gz',args.tempdir+'/novel_locus_reads.gpd.gz',args)
  # Now we can make a new non-redundant set of genpreds from the novel isoforms
  sys.stderr.write("making NR novel isoforms\n")
  cmd = 'gpd_to_nr.py '+args.tempdir+'/novel_isoform_reads.gpd.gz '+\
        ' -j '+str(args.junction_tolerance)+' --threads '+str(args.threads)+\
        ' --minimum_junction_end_support '+str(args.minimum_junction_end_support)+\
        ' --minimum_support '+str(args.minimum_support)+\
        ' --gene_names '+\
        ' -o '+args.tempdir+'/novel_isoforms_nr.gpd.gz'
  gpd_to_nr.external_cmd(cmd)

  sys.stderr.write("reannotating novel based on our new gpd\n")
  # Now we reannotate the novel based on the these newly annotated isoforms
  cmd = 'gpd_anntotate.py '+args.tempdir+'/novel_locus_reads.gpd.gz '+\
        ' --threads '+str(1)+' '+\
        ' -r '+args.tempdir+'/novel_isoforms_nr.gpd.gz '+\
        ' -o '+args.tempdir+'/novel_locus_reads.annot.txt.gz'
  gpd_annotate.external_cmd(cmd)
  
  # now this new annotation should be classified
  # the new isoform will be in novel_isoform_reads.gpd.gz
  cmd = 'classify_reads.py '+args.tempdir+'/novel_locus_reads.annot.txt.gz '+args.tempdir+'/novel_locus_reads.gpd.gz -o '+args.tempdir+'/classify_novel.txt.gz'
  sys.stderr.write(cmd+"\n")
  classify_reads.external_cmd(cmd)
  get_novel_sets(args.tempdir+'/classify_novel.txt.gz',args.tempdir+'/novel_locus_reads.gpd.gz',args.tempdir+'/novel_isoform_reads2.gpd.gz',args.tempdir+'/novel_locus_reads2.gpd.gz',args)

  # now lets combine our novel isoform reads making sure to sort them
  of = open(args.tempdir+'/new_novel_isoform_reads.gpd.gz','w')
  cmd2 = 'gzip'
  p2 = Popen(cmd2.split(),stdout=of,stdin=PIPE)
  cmd1 = 'sort -k3,3 -k5,5n -k6,6n'
  p1 = Popen(cmd1.split(),stdout=p2.stdin,stdin=PIPE)
  inf = gzip.open(args.tempdir+'/novel_isoform_reads.gpd.gz')
  for line in inf:  p1.stdin.write(line)
  inf.close()
  inf = gzip.open(args.tempdir+'/novel_isoform_reads2.gpd.gz')
  for line in inf: p1.stdin.write(line)
  inf.close()
  p1.communicate()
  p2.communicate()
  of.close()
  
  # Now we can make a new non-redundant set of genpreds from the novel isoforms
  sys.stderr.write("making NR novel isoforms\n")
  cmd = 'gpd_to_nr.py '+args.tempdir+'/new_novel_isoform_reads.gpd.gz '+\
        ' -j '+str(args.junction_tolerance)+' --threads '+str(args.threads)+\
        ' --minimum_junction_end_support '+str(args.minimum_junction_end_support)+\
        ' --minimum_support '+str(args.minimum_support)+\
        ' --gene_names '+\
        ' -o '+args.tempdir+'/novel_isoforms_nr2.gpd.gz'
  gpd_to_nr.external_cmd(cmd)

  #Only need to reannotate if we are interested in whats left over
  #sys.stderr.write("reannotating novel based on our new gpd\n")
  ## Now we reannotate the novel based on the these newly annotated isoforms
  #cmd = 'gpd_anntotate.py '+args.tempdir+'/novel_locus_reads.gpd.gz '+\
  #      ' --threads '+str(args.threads)+' '+\
  #      ' -r '+args.tempdir+'/novel_isoforms_nr2.gpd.gz '+\
  #      ' -o '+args.tempdir+'/novel_locus_reads.annot.txt.gz'
  #gpd_annotate.external_cmd(cmd)  

  sys.stderr.write("now work on the novel loci\n")
  # Now lets work on the novel locus
  of = open(args.tempdir+'/sorted_novel_locus_reads.gpd.gz','w')
  cmd2 = 'gzip'
  p2 = Popen(cmd2.split(),stdout=of,stdin=PIPE)
  cmd1 = 'sort -k3,3 -k5,5n -k6,6n'
  p1 = Popen(cmd1.split(),stdout=p2.stdin,stdin=PIPE)
  inf = gzip.open(args.tempdir+'/novel_locus_reads2.gpd.gz')
  for line in inf:  p1.stdin.write(line)
  inf.close()
  p1.communicate()
  p2.communicate()
  of.close()

  sys.stderr.write("making NR novel loci\n")
  cmd = 'gpd_to_nr.py '+args.tempdir+'/sorted_novel_locus_reads.gpd.gz '+\
        ' -j '+str(args.junction_tolerance)+' --threads '+str(args.threads)+\
        ' --minimum_junction_end_support '+str(args.minimum_junction_end_support)+\
        ' --minimum_support '+str(args.minimum_support)+\
        ' -o '+args.tempdir+'/novel_locus_nr.gpd.gz'
  gpd_to_nr.external_cmd(cmd)

  sys.stderr.write("sort the novel isoforms\n")
  of = open(args.tempdir+'/novel_isoforms_nr.sorted.gpd.gz','w')
  cmd2 = 'gzip'
  p2 = Popen(cmd2.split(),stdout=of,stdin=PIPE)
  cmd1 = 'sort -k3,3 -k5,5n -k6,6n'
  p1 = Popen(cmd1.split(),stdout=p2.stdin,stdin=PIPE)
  inf = gzip.open(args.tempdir+'/novel_isoforms_nr2.gpd.gz')
  for line in inf:  p1.stdin.write(line)
  inf.close()
  p1.communicate()
  p2.communicate()
  of.close()

  sys.stderr.write("sort the novel loci\n")
  of = open(args.tempdir+'/novel_loci_nr.sorted.gpd.gz','w')
  cmd2 = 'gzip'
  p2 = Popen(cmd2.split(),stdout=of,stdin=PIPE)
  cmd1 = 'sort -k3,3 -k5,5n -k6,6n'
  p1 = Popen(cmd1.split(),stdout=p2.stdin,stdin=PIPE)
  inf = gzip.open(args.tempdir+'/novel_locus_nr.gpd.gz')
  for line in inf:  p1.stdin.write(line)
  inf.close()
  p1.communicate()
  p2.communicate()
  of.close()

  # Now we can rename totally novel genes based on locus overlap
  of = open(args.tempdir+'/novel_loci_nr_named.sorted.gpd.gz','w')
  cmd2 = 'gzip'
  p2 = Popen(cmd2.split(),stdout=of,stdin=PIPE)
  cmd1 = 'sort -k3,3 -k5,5n -k6,6n'
  p1 = Popen(cmd1.split(),stdout=p2.stdin,stdin=PIPE)

  inf = gzip.open(args.tempdir+'/novel_loci_nr.sorted.gpd.gz')
  gs = GPDStream(inf)
  ls = LocusStream(gs)
  z = 0
  for rng in ls:
    z+=1
    rng_string = rng.get_range_string()
    gpds = rng.get_payload()
    for gpd in gpds:
      gene_name = 'LOC'+str(z)+'|'+str(len(gpds))+'|'+rng_string
      f = gpd.get_gpd_line().rstrip().split("\t")
      f[0] = gene_name
      gpd_line = "\t".join(f)
      p1.stdin.write(gpd_line+"\n")
  p1.communicate()
  p2.communicate()
  of.close()

  # we are almost done but we need to make sure these genepreds aren't subsets of known genes
  sys.stderr.write("reannotating novel-isoform by reference\n")
  cmd = 'gpd_anntotate.py '+args.tempdir+'/novel_isoforms_nr.sorted.gpd.gz '+\
        ' --threads '+str(1)+' '+\
        ' -r '+args.reference_annotation_gpd+\
        ' -o '+args.tempdir+'/novel_isoforms_nr.annot.txt.gz'
  gpd_annotate.external_cmd(cmd)  
  cmd = 'classify_reads.py '+args.tempdir+'/novel_isoforms_nr.annot.txt.gz '+args.tempdir+'/novel_isoforms_nr.sorted.gpd.gz -o '+args.tempdir+'/classify_novel_isoform_ref.txt.gz'
  sys.stderr.write(cmd+"\n")
  classify_reads.external_cmd(cmd)

  # now we can screen to make sure things in the novel isoform file really are novel isoforms
  blacklist = set()
  finf = gzip.open(args.tempdir+'/classify_novel_isoform_ref.txt.gz')
  for line in finf:
    f = line.rstrip().split("\t")
    if f[2]=='subset' or f[2]=='full': blacklist.add(f[0])
  finf.close()
  fof = gzip.open(args.tempdir+'/novel_isoforms_nr.filtered.sorted.gpd.gz','w')
  finf = gzip.open(args.tempdir+'/novel_isoforms_nr.sorted.gpd.gz')
  for line in finf:
    f = line.rstrip().split("\t")
    if f[1] in blacklist: continue
    fof.write(line)
  finf.close()
  fof.close()


  sys.stderr.write("reannotating novel-locus by reference\n")
  cmd = 'gpd_anntotate.py '+args.tempdir+'/novel_loci_nr_named.sorted.gpd.gz '+\
        ' --threads '+str(1)+' '+\
        ' -r '+args.reference_annotation_gpd+\
        ' -o '+args.tempdir+'/novel_loci_nr_named.annot.txt.gz'
  gpd_annotate.external_cmd(cmd)  
  cmd = 'classify_reads.py '+args.tempdir+'/novel_loci_nr_named.annot.txt.gz '+args.tempdir+'/novel_loci_nr_named.sorted.gpd.gz -o '+args.tempdir+'/classify_novel_loci.txt.gz'
  sys.stderr.write(cmd+"\n")
  classify_reads.external_cmd(cmd)

  # now we can screen to make sure things in the novel isoform file really are novel isoforms
  blacklist = set()
  finf = gzip.open(args.tempdir+'/classify_novel_loci.txt.gz')
  for line in finf:
    f = line.rstrip().split("\t")
    if f[2]=='subset' or f[2]=='full': blacklist.add(f[0])
  finf.close()
  fof = gzip.open(args.tempdir+'/novel_loci_nr_named.filtered.sorted.gpd.gz','w')
  finf = gzip.open(args.tempdir+'/novel_loci_nr_named.sorted.gpd.gz')
  for line in finf:
    f = line.rstrip().split("\t")
    if f[1] in blacklist: continue
    fof.write(line)
  finf.close()
  fof.close()




  if not os.path.exists(args.output):
    os.makedirs(args.output)

  copy(args.tempdir+'/novel_loci_nr_named.filtered.sorted.gpd.gz',args.output+'/novel_loci_nr_named.sorted.gpd.gz')
  copy(args.tempdir+'/novel_isoforms_nr.filtered.sorted.gpd.gz',args.output+'/novel_isoforms_nr.sorted.gpd.gz')
  
  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

def get_novel_sets(classification,input_gpd,out_iso,out_locus,args):
  # now we want to create a non redundant version of the novel isoforms
  novel_isoforms = set()
  novel_isoform_genes = {}
  novel_loci = set()
  inf = gzip.open(classification)
  for line in inf:
    f = line.rstrip().split("\t")
    if f[2] == 'novel-isoform':
      novel_isoforms.add(f[0])
      novel_isoform_genes[f[0]]=f[1] # save the gene name
    elif f[2] == 'novel-locus':
      novel_loci.add(f[0])
  inf.close()
  sys.stderr.write("outputing novel isoforms to a file\n")
  tof = gzip.open(out_iso,'w')
  lof = gzip.open(out_locus,'w')
  inf_gpd = None;
  if input_gpd[-3:]=='.gz':
    inf_gpd = gzip.open(input_gpd)
  else:
    inf_gpd = open(input_gpd)
  z = 0
  for line in inf_gpd:
    z += 1
    if z % 1000 == 0: sys.stderr.write(str(z)+" reads processed\r")
    f = line.rstrip().split("\t")
    if f[1] in novel_isoforms:
      f[0] = novel_isoform_genes[f[0]]
      newline = "\t".join(f)
      tof.write(newline+"\n")
    elif f[1] in novel_loci:
      lof.write(line)  
  inf_gpd.close()
  tof.close()
  lof.close()
  sys.stderr.write("\n")

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input_annot',help="<input annotbest.txt>")
  parser.add_argument('input_gpd',help="<input best.sorted.gpd>")
  parser.add_argument('-a','--reference_annotation_gpd',required=True,help="Reference annotation GPD")
  parser.add_argument('-o','--output',required=True,help="OUTPUT DIRECTORY")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")

  # Run parameters
  parser.add_argument('-j','--junction_tolerance',type=int,default=10,help="number of bp to tolerate junction mismatch on either side")
  parser.add_argument('--minimum_junction_end_support',type=int,default=2,help="minimum coverage of end exons")
  parser.add_argument('--minimum_support',type=int,default=2,help="minimum supporting reads")
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
