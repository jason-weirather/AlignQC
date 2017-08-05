#!/usr/bin/env python
"""This script creates reference bed files with locations of mutually exclusive
   genomic features
   1) exons 2) introns 3) intergenic regions
 
   It also annotates reads according to these features.
 
   PRE:  Requires the read mappings in gpd format
        Requires a transcript reference annotation in genepred format
        Requires the chromosome lengths in TSV format
   POST: Output a folder of beds describing the locations of genomic features
        Output a file describing the membership of each read
"""
import sys, argparse, gzip, re, os, inspect, itertools
from multiprocessing import Pool, cpu_count

from seqtools.format.gpd import GPDStream
from seqtools.range.multi import merge_ranges, subtract_ranges, BedArrayStream, sort_ranges
from seqtools.range import GenomicRange
from seqtools.stream import MultiLocusStream

def main(args):

  inf = None
  chrlens = {}
  chrbed = []
  if re.search('\.gz$',args.chromosome_lengths): 
    inf = gzip.open(args.chromosome_lengths)
  else:
    inf = open(args.chromosome_lengths)
  for line in inf:
    f = line.rstrip().split("\t")
    chrlens[f[0]] = int(f[1])
    chrbed.append(GenomicRange(f[0],1,int(f[1])))
  inf.close()  

  inf = None
  exonbed = []
  txbed = []
  sys.stderr.write("Reading Exons\n")
  if re.search('\.gz$',args.annotation_gpd): 
    inf = gzip.open(args.annotation_gpd)
  else:
    inf = open(args.annotation_gpd)
  gs = GPDStream(inf)
  for gpd in gs:
    exonbed += [x.range for x in gpd.exons]
    txbed.append(gpd.range)
  inf.close()  
  sys.stderr.write("Merging "+str(len(txbed))+" transcripts\n")
  txbed = merge_ranges(txbed)
  sys.stderr.write(str(len(txbed))+" transcripts after merging\n")
  sys.stderr.write("Finding intergenic\n")
  intergenicbed = subtract_ranges(chrbed,txbed)
  sys.stderr.write("Found "+str(len(intergenicbed))+" intergenic regions\n")
  intergenicbp = sum([x.length for x in intergenicbed])
  sys.stderr.write("Intergenic size: "+str(intergenicbp)+"\n")
  sys.stderr.write("Merging "+str(len(exonbed))+" exons\n")
  exonbed = merge_ranges(exonbed)
  sys.stderr.write(str(len(exonbed))+" exons after merging\n")
  sys.stderr.write("Finding introns\n")
  intronbed = subtract_ranges(txbed,exonbed)
  sys.stderr.write("Found "+str(len(intronbed))+" introns\n")
  
  chrbp = sum([x.length for x in chrbed])
  sys.stderr.write("Genome size: "+str(chrbp)+"\n")
  txbp = sum([x.length for x in txbed])
  sys.stderr.write("Tx size: "+str(txbp)+"\n")
  exonbp = sum([x.length for x in exonbed])
  sys.stderr.write("Exon size: "+str(exonbp)+"\n")
  intronbp = sum([x.length for x in intronbed])
  sys.stderr.write("Intron size: "+str(intronbp)+"\n")
  #sys.stderr.write(str(txbp+intergenicbp)+"\n")

  if args.output_beds:
    if not os.path.exists(args.output_beds): os.makedirs(args.output_beds)
    with open(args.output_beds+'/chrs.bed','w') as of1:
      for rng in chrbed: of1.write("\t".join([str(x) for x in rng.get_bed_array()])+"\n")
    with open(args.output_beds+'/exon.bed','w') as of1:
      for rng in exonbed: of1.write("\t".join([str(x) for x in rng.get_bed_array()])+"\n")
    with open(args.output_beds+'/intron.bed','w') as of1:
      for rng in intronbed: of1.write("\t".join([str(x) for x in rng.get_bed_array()])+"\n")
    with open(args.output_beds+'/intergenic.bed','w') as of1:
      for rng in intergenicbed: of1.write("\t".join([str(x) for x in rng.get_bed_array()])+"\n")
    with open(args.output_beds+'/tx.bed','w') as of1:
      for rng in txbed: of1.write("\t".join([str(x) for x in rng.get_bed_array()])+"\n")

  inf = None
  if re.search('\.gz$',args.reads_gpd): 
    inf = gzip.open(args.reads_gpd)
  else:
    inf = open(args.reads_gpd)
  reads = {}

  gs = GPDStream(inf)
  for gpd in gs:
    reads[gpd.gene_name] = {}
  sys.stderr.write("Checking "+str(len(reads.keys()))+" Aligned Reads\n")

  #now we know all features we can annotate reads

  sys.stderr.write("Read through our reads and bed entries\n")
  sys.stderr.write("Annotate intron\n")
  intron = annotate_gpds(args,intronbed)
  intronnames = set(intron.keys())
  sys.stderr.write("Annotate intergenic\n")
  intergenic = annotate_gpds(args,intergenicbed)
  intergenicnames = set(intergenic.keys())
  sys.stderr.write("Annotate exons\n")
  exons = annotate_gpds(args,exonbed)
  exonnames = set(exons.keys())
  allnames = exonnames|intronnames|intergenicnames
  sys.stderr.write(str(len(allnames))+" reads attributed to a feature\n")
  vals = set(reads.keys())-allnames
  if len(vals) > 0:
    sys.stderr.write("WARNING unable to ascribe annotation to "+str(len(vals))+" reads\n")
  donenames = set()
  of = sys.stdout
  if args.output:
    if re.search('\.gz$',args.output):
      of = gzip.open(args.output,'w')
    else:
      of = open(args.output,'w')
  for name in allnames:
    exonfrac = 0
    intronfrac = 0
    intergenicfrac = 0
    readlen = 0
    exoncount = 0
    if name in exons:  
      exonfrac = float(exons[name][1])/float(exons[name][0])
      readlen = exons[name][0]
      exoncount = exons[name][2]
    if name in intron:  
      intronfrac = float(intron[name][1])/float(intron[name][0])
      readlen = intron[name][0]
      exoncount = intron[name][2]
    if name in intergenic:  
      intergenicfrac = float(intergenic[name][1])/float(intergenic[name][0])
      readlen = intergenic[name][0]
      exoncount = intergenic[name][2]
    vals = {'exon':exonfrac,'intron':intronfrac,'intergenic':intergenicfrac}
    type = None
    if exonfrac >= 0.5: 
      type = 'exon'
    elif intronfrac >= 0.5:
      type = 'intron'
    elif intergenicfrac >= 0.5:
      type = 'intergenic'
    else:
      type = sorted(vals.keys(),key=lambda x: vals[x])[-1]
      if vals[type] == 0: 
        sys.stderr.write("WARNING trouble setting type\n")
    if not type: continue
    of.write(name+"\t"+type+"\t"+str(exoncount)+"\t"+str(readlen)+"\n")
  of.close()

def generate_locus(mls):
  for es in mls:
    [gpds,inbeds] = es.payload
    if len(gpds) == 0 or len(inbeds) == 0:
      continue
    yield es

def annotate_gpds(args,inputbed):
  if args.threads > 1:
    p = Pool(processes=args.threads)
  bas = BedArrayStream(sort_ranges(inputbed))
  inf = None
  if re.search('\.gz$',args.reads_gpd): 
    inf = gzip.open(args.reads_gpd)
  else:
    inf = open(args.args.reads_gpd)
  gs = GPDStream(inf)
  mls = MultiLocusStream([gs,bas])
  results = {}
  # try and implement as a multiprocessing map function
  csize = 100 #control how many jobs to send to one thread at a time
  if args.threads > 1:
    results2 = p.imap_unordered(func=annotate_inner,iterable=generate_locus(mls),chunksize=csize)
  else:
    results2 = itertools.imap(annotate_inner,generate_locus(mls))
  for chunk in results2:
    for res in chunk:
      results[res[0]] = res[1:]
  inf.close()
  return results

def annotate_inner(es):
  results = []
  [gpds,inbeds] = es.payload
  for gpd in gpds:
    orig = gpd.length
    tot = 0
    for rng1 in [x.range for x in gpd.exons]:
      tot += sum([y.overlap_size(rng1) for y in inbeds])
    if tot > 0:
      results.append([gpd.gene_name,orig,tot,gpd.get_exon_count()])
  return results

def do_inputs():
  parser = argparse.ArgumentParser(description="Assign genomic features to reads based on where they majority of the read lies.  In the event of a tie prioritize exon over intron and intron over intergenic.",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('reads_gpd',help="reads gpd")
  parser.add_argument('annotation_gpd',help="reference annotations gpd")
  parser.add_argument('chromosome_lengths',help="reference lengths table")
  parser.add_argument('--output_beds',help="save features")
  parser.add_argument('-o','--output',help="output results")
  parser.add_argument('--threads',default=cpu_count(),type=int,help="number of threads default cpu_count()")
  args = parser.parse_args()
  return args
  

def external_cmd(cmd):
  cache_argv = sys.argv
  sys.argv = cmd
  args = do_inputs()
  main(args)
  sys.argv = cache_argv

if __name__=="__main__":
  args = do_inputs()
  main(args)
