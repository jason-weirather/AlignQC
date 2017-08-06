#!/usr/bin/env python
"""This script will call most of the individual modules analyzing the data"""

import argparse, sys, os, time, re, gzip, locale, inspect, time
from subprocess import Popen, PIPE


# BAM imports
import bam_preprocess
import traverse_preprocessed
import bam_to_chr_lengths
import get_platform_report
import gpd_loci_analysis
import gpd_to_exon_distro
import make_alignment_plot
import depth_to_coverage_report
import locus_bed_to_rarefraction

# BAM + reference imports
import bam_to_context_error_plot
import bam_to_alignment_error_plot

# BAM + annotation
import  annotate_from_genomic_features
import  get_depth_subset
import  annotated_length_analysis
import  gpd_annotation_to_rarefraction
import  annotated_read_bias_analysis
import  gpd_to_junction_variance

# BAM
from seqtools.cli.utilities.gpd_to_bed_depth import external_cmd as gpd_to_bed_depth

from seqtools.cli.utilities.bed_depth_to_stratified_coverage import external_cmd as bed_depth_to_stratified_coverage

from seqtools.cli.utilities.gpd_to_UCSC_bed12 import external_cmd as gpd_to_UCSC_bed12

# BAM + annotation
from seqtools.cli.utilities.gpd_annotate import external_cmd as gpd_annotate

# read count
rcnt = -1
tlog = None

def main(args):

  if not os.path.exists(args.tempdir+'/plots'):
    os.makedirs(args.tempdir+'/plots')
  if not os.path.exists(args.tempdir+'/data'):
    os.makedirs(args.tempdir+'/data')
  if not os.path.exists(args.tempdir+'/logs'):
    os.makedirs(args.tempdir+'/logs')
  if not os.path.exists(args.tempdir+'/temp'):
    os.makedirs(args.tempdir+'/temp')

  global tlog
  tlog = TimeLog(args.tempdir+'/logs/time.log')

  ## Extract data that can be realized from the bam
  make_data_bam(args)


  ## Extract data that can be realized from the bam and reference
  if args.genome:
    make_data_bam_reference(args)

  ## Extract data that can be realized from bam and reference annotation
  if args.gpd:
    make_data_bam_annotation(args)

  # Write params file
  of = open(args.tempdir+'/data/params.txt','w')
  for arg in vars(args):
    of.write(arg+"\t"+str(getattr(args,arg))+"\n")
  of.close()

class TimeLog:
  def __init__(self,fname):
    self.fh = open(fname,'w')
    self.recording = False
    self.st = time.time()
  def start(self,msg):
    self.st = time.time()
    self.fh.write(msg+"\n")
    self.fh.flush()
  def write(self,msg):
    self.fh.write('$ '+msg+"\n")
  def stop(self):
    self.fh.write("--- "+str(time.time()-self.st)+ " seconds ---\n")
    self.fh.flush()

def make_data_bam(args):
  global tlog
  # Get the data necessary for making tables and reports

  tlog.start("traverse bam and preprocess")
  # 1. Traverse bam and store alignment mappings ordered by query name
  udir = os.path.dirname(os.path.realpath(__file__))
  cmd =  [udir+'/bam_preprocess.py',args.input,'--minimum_intron_size',
          str(args.min_intron_size),'-o',args.tempdir+'/temp/alndata.txt.gz',
          '--threads',str(args.threads),'--specific_tempdir',
          args.tempdir+'/temp/']
  sys.stderr.write("Creating initial alignment mapping data\n")
  sys.stderr.write(" ".join(cmd)+"\n")
  bam_preprocess.external_cmd(cmd)
  tlog.write(" ".join(cmd))
  tlog.stop()
  tlog.start("traverse preprocessed")
  # 2. Describe the alignments by traversing the previously made file
  cmd = [udir+'/traverse_preprocessed.py',args.tempdir+'/temp/alndata.txt.gz',
         '-o',args.tempdir+'/data/','--specific_tempdir',args.tempdir+'/temp/',
         '--threads',str(args.threads)]
  if args.min_aligned_bases:
    cmd += ['--min_aligned_bases',str(args.min_aligned_bases)]
  if args.max_query_overlap:
    cmd += ['--max_query_overlap',str(args.max_query_overlap)]
  if args.max_target_overlap:
    cmd += ['--max_target_overlap',str(args.max_target_overlap)]
  if args.max_query_gap:
    cmd += ['--max_query_gap',str(args.max_query_gap)]
  if args.max_target_gap:
    cmd += ['--max_target_gap',str(args.max_target_gap)]
  if args.required_fractional_improvement:
    cmd += ['--required_fractional_improvement',
            str(args.required_fractional_improvement)]
  sys.stderr.write("Traverse bam for alignment analysis\n")
  sys.stderr.write(" ".join(cmd)+"\n")
  traverse_preprocessed.external_cmd(cmd)
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("get chromosome lengths from bam")
  # 3. get chromosome lengths from bam
  cmd = [udir+'/bam_to_chr_lengths.py',args.input,
         '-o',args.tempdir+'/data/chrlens.txt']
  sys.stderr.write("Writing chromosome lengths from header\n")
  sys.stderr.write(" ".join(cmd)+"\n")
  bam_to_chr_lengths.external_cmd(cmd)
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("look for platform specific read names")
  # Now we can find any known reads
  # 4. Go through read names to find if there are platform-specific read names present
  sys.stderr.write("Can we find any known read types\n")
  cmd = [udir+'/get_platform_report.py',args.tempdir+'/data/lengths.txt.gz',
         args.tempdir+'/data/special_report']
  sys.stderr.write(" ".join(cmd)+"\n")
  get_platform_report.external_cmd(cmd)
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("check for pacbio in case we make a special graph")
  # Check for pacbio to see if we need to make a graph for it
  do_pb = False
  with open(args.tempdir+'/data/special_report') as inf:
    for line in inf:
      f = line.rstrip().split("\t")
      if f[0]=='PB':
        do_pb = True
        break
  if do_pb:
    cmd = [args.rscript_path,udir+'/plot_pacbio.r',
           args.tempdir+'/data/special_report.pacbio',
           args.tempdir+'/plots/pacbio.png']
    sys.stderr.write(" ".join(cmd)+"\n")
    mycall(cmd,args.tempdir+'/logs/special_report_pacbio_png')
    cmd = [args.rscript_path,udir+'/plot_pacbio.r',
           args.tempdir+'/data/special_report.pacbio',
           args.tempdir+'/plots/pacbio.pdf']
    sys.stderr.write(" ".join(cmd)+"\n")
    mycall(cmd,args.tempdir+'/logs/special_report_pacbio_pdf')
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("get a depth for our best alignments")
  # 5. Go through the genepred file and get a depth bed for our best alignments
  sys.stderr.write("Go through genepred best alignments and make a bed depth file\n")
  cmd = ["gpd_to_bed_depth.py",args.tempdir+'/data/best.sorted.gpd.gz',
         '-o',args.tempdir+'/data/depth.sorted.bed.gz',"--threads",
         str(args.threads)]
  sys.stderr.write("Generate the depth bed for the mapped reads\n")
  sys.stderr.write(" ".join(cmd)+"\n")
  gpd_to_bed_depth(cmd)

  sys.stderr.write("Stratify the depth to make it plot quicker and cleaner\n")
  cmd = ["bed_depth_to_stratified_coverage.py",
         args.tempdir+'/data/depth.sorted.bed.gz','-l',
         args.tempdir+"/data/chrlens.txt",
         '-o',args.tempdir+'/temp/depth.coverage-strata.sorted.bed.gz',
         '--output_key',args.tempdir+'/temp/coverage-strata.key',
         '--minimum_coverage','100000']
  bed_depth_to_stratified_coverage(cmd)

  global rcnt #read count
  rcnt = 0
  tinf = gzip.open(args.tempdir+'/data/lengths.txt.gz')
  for line in tinf:  rcnt += 1
  tinf.close()
  tlog.write(" ".join(cmd))
  tlog.stop()

  # For now reporting loci will be optional until it can be tested and optimized.
  if args.do_loci:
    tlog.start("do locus search")
    # 6. Go through the best alignments and look for loci
    sys.stderr.write("Approximate loci and mapped read distributions among them.\n")
    cmd = [udir+"/gpd_loci_analysis.py",
           args.tempdir+'/data/best.sorted.gpd.gz','-o',
           args.tempdir+'/data/loci-all.bed.gz','--output_loci',
           args.tempdir+'/data/loci.bed.gz','--downsample',
           str(args.locus_downsample),'--threads',str(args.threads)]
    if args.min_depth:
      cmd += ['--min_depth',str(args.min_depth)]
    if args.min_depth:
      cmd += ['--min_coverage_at_depth',str(args.min_coverage_at_depth)]
    if args.min_exon_count:
      cmd += ['--min_exon_count',str(args.min_exon_count)]
    sys.stderr.write(" ".join(cmd)+"\n")
    gpd_loci_analysis.external_cmd(cmd)

    cmd = [udir+"/locus_bed_to_rarefraction.py",
           args.tempdir+'/data/loci.bed.gz','-o',
           args.tempdir+'/data/locus_rarefraction.txt','--threads',
           str(args.threads),'--original_read_count',str(rcnt)]
    sys.stderr.write("Make rarefraction curve\n")
    sys.stderr.write(" ".join(cmd)+"\n")
    locus_bed_to_rarefraction.external_cmd(cmd)

    sys.stderr.write("Make locus rarefraction plot\n")
    for ext in ['png','pdf']:
      cmd = [args.rscript_path,udir+'/plot_annotation_rarefractions.r',
             args.tempdir+'/plots/locus_rarefraction.'+ext,'locus',
             args.tempdir+'/data/locus_rarefraction.txt','#FF000088']
      sys.stderr.write(" ".join(cmd)+"\n")
      mycall(cmd,args.tempdir+'/logs/plot_locus_rarefraction_'+ext)
    tlog.write(" ".join(cmd))
    tlog.stop()

  tlog.start("get ready for alignment plot")
  # 7. Alignment plot preparation
  sys.stderr.write("Get ready for alignment plot\n")
  cmd = [udir+'/make_alignment_plot.py',args.tempdir+'/data/lengths.txt.gz',
         '--rscript_path',args.rscript_path,'--output_stats',
         args.tempdir+'/data/alignment_stats.txt','--output',
         args.tempdir+'/plots/alignments.png',
         args.tempdir+'/plots/alignments.pdf']
  sys.stderr.write("Make alignment plots\n")
  sys.stderr.write(" ".join(cmd)+"\n")
  make_alignment_plot.external_cmd(cmd)
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("make depth reports")
  # 8. Make depth reports
  sys.stderr.write("Making depth reports\n")
  cmd = [udir+'/depth_to_coverage_report.py',
         args.tempdir+'/data/depth.sorted.bed.gz',
         args.tempdir+'/data/chrlens.txt','-o',args.tempdir+'/data']
  sys.stderr.write(" ".join(cmd)+"\n")
  depth_to_coverage_report.external_cmd(cmd)
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("make coverage plots")
  # do the depth graphs
  sys.stderr.write("Making coverage plots\n")
  cmd = [args.rscript_path,udir+'/plot_chr_depth.r',
         args.tempdir+'/data/line_plot_table.txt.gz',
         args.tempdir+'/data/total_distro_table.txt.gz',
         args.tempdir+'/data/chr_distro_table.txt.gz',
         args.tempdir+'/plots/covgraph.png']
  sys.stderr.write(" ".join(cmd)+"\n")
  mycall(cmd,args.tempdir+'/logs/covgraph_png')
  cmd = [args.rscript_path,udir+'/plot_chr_depth.r',
         args.tempdir+'/data/line_plot_table.txt.gz',
         args.tempdir+'/data/total_distro_table.txt.gz',
         args.tempdir+'/data/chr_distro_table.txt.gz',
         args.tempdir+'/plots/covgraph.pdf']
  sys.stderr.write(" ".join(cmd)+"\n")
  mycall(cmd,args.tempdir+'/logs/covgraph_pdf')
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("make chr depth plots")
  # do depth plots
  sys.stderr.write("Making chr depth plots\n")
  cmd = [args.rscript_path,udir+'/plot_depthmap.r',
         args.tempdir+'/temp/depth.coverage-strata.sorted.bed.gz',
         args.tempdir+'/data/chrlens.txt',
         args.tempdir+'/temp/coverage-strata.key',
         args.tempdir+'/plots/perchrdepth.png']
  sys.stderr.write(" ".join(cmd)+"\n")
  mycall(cmd,args.tempdir+'/logs/perchr_depth_png')
  cmd = [args.rscript_path,udir+'/plot_depthmap.r',
         args.tempdir+'/temp/depth.coverage-strata.sorted.bed.gz',
         args.tempdir+'/data/chrlens.txt',
         args.tempdir+'/temp/coverage-strata.key',
         args.tempdir+'/plots/perchrdepth.pdf']
  sys.stderr.write(" ".join(cmd)+"\n")
  mycall(cmd,args.tempdir+'/logs/perchr_depth_pdf')
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("get the exon size distribution")
  #Get the exon distribution
  sys.stderr.write("Get the exon distributions\n")
  cmd = [udir+'/gpd_to_exon_distro.py',args.tempdir+'/data/best.sorted.gpd.gz',
         '-o',args.tempdir+'/data/exon_size_distro.txt.gz','--threads',
         str(args.threads)]
  sys.stderr.write(" ".join(cmd)+"\n")
  gpd_to_exon_distro.external_cmd(cmd)
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("plot exon distro")
  cmd = [args.rscript_path,udir+'/plot_exon_distro.r',
         args.tempdir+'/data/exon_size_distro.txt.gz',
         args.tempdir+'/plots/exon_size_distro.png']
  sys.stderr.write(" ".join(cmd)+"\n")
  mycall(cmd,args.tempdir+'/logs/exon_size_distro_png')
  cmd = [args.rscript_path,udir+'/plot_exon_distro.r',
         args.tempdir+'/data/exon_size_distro.txt.gz',
         args.tempdir+'/plots/exon_size_distro.pdf']
  sys.stderr.write(" ".join(cmd)+"\n")
  mycall(cmd,args.tempdir+'/logs/exon_size_distro_pdf')
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("make bed file")
  # Make a UCSC compatible bed file
  sys.stderr.write("Make a UCSC genome browser compatible bed file\n")
  cmd = ['gpd_to_UCSC_bed12.py','--headername',args.input+':best',
         args.tempdir+'/data/best.sorted.gpd.gz','-o',
         args.tempdir+'/data/best.sorted.bed.gz','--color','red']
  sys.stderr.write(" ".join(cmd)+"\n")
  gpd_to_UCSC_bed12(cmd)
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("make bed file")
  cmd = ['gpd_to_UCSC_bed12.py','--headername',args.input+':trans-chimera',
        args.tempdir+'/data/chimera.gpd.gz','-o',
        args.tempdir+'/data/chimera.bed.gz','--color','blue']
  sys.stderr.write(" ".join(cmd)+"\n")
  gpd_to_UCSC_bed12(cmd)
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("make bed file")
  cmd = ['gpd_to_UCSC_bed12.py','--headername',args.input+':gapped',
         args.tempdir+'/data/gapped.gpd.gz','-o',
         args.tempdir+'/data/gapped.bed.gz','--color','orange']
  sys.stderr.write(" ".join(cmd)+"\n")
  gpd_to_UCSC_bed12(cmd)
  tlog.write(" ".join(cmd))
  tlog.stop()

  cmd = ['gpd_to_UCSC_bed12.py','--headername',args.input+':self-chimera',
         args.tempdir+'/data/technical_chimeras.gpd.gz','-o',
         args.tempdir+'/data/technical_chimeras.bed.gz','--color','green']
  sys.stderr.write(" ".join(cmd)+"\n")
  gpd_to_UCSC_bed12(cmd)
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("make bed file")
  cmd = ['gpd_to_UCSC_bed12.py','--headername',args.input+':self-atypical',
         args.tempdir+'/data/technical_atypical_chimeras.gpd.gz','-o',
         args.tempdir+'/data/technical_atypical_chimeras.bed.gz','--color',
         'purple']
  sys.stderr.write(" ".join(cmd)+"\n")
  gpd_to_UCSC_bed12(cmd)
  tlog.write(" ".join(cmd))
  tlog.stop()

def make_data_bam_reference(args):
  global tlog

  # make the context error plots
  udir = os.path.dirname(os.path.realpath(__file__))

  # Find the index file that was generated earlier.
  indfile = None
  if os.path.exists(args.tempdir+'/temp/myindex.bgi'):
    indfile = args.tempdir+'/temp/myindex.bgi'

  tlog.start("Get context error")
  # 1. Context error
  cmd = [udir+'/bam_to_context_error_plot.py',args.input,'-r',args.genome,
         '--target','--output_raw',args.tempdir+'/data/context_error_data.txt',
         '-o',args.tempdir+'/plots/context_plot.png',
         args.tempdir+'/plots/context_plot.pdf','--rscript_path',
         args.rscript_path,'--random','--specific_tempdir',
         args.tempdir+'/temp']
  if args.context_error_scale:
    cmd += ['--scale']+[str(x) for x in args.context_error_scale]
  if args.context_error_stopping_point:
    cmd += ['--stopping_point',str(args.context_error_stopping_point)]
  if indfile:
    cmd += ['--input_index',indfile]
  sys.stderr.write("Making context plot\n")
  sys.stderr.write(" ".join(cmd)+"\n")
  bam_to_context_error_plot.external_cmd(cmd)
  tlog.write(" ".join(cmd))
  tlog.stop()
  time.sleep(3)
  #gc.collect()

  tlog.start("alignment based error")
  # 2. Alignment overall error
  cmd = [udir+'/bam_to_alignment_error_plot.py',args.input,'-r',
         args.genome,'--output_stats',args.tempdir+'/data/error_stats.txt',
         '--output_raw',args.tempdir+'/data/error_data.txt','-o',
         args.tempdir+'/plots/alignment_error_plot.png',
         args.tempdir+'/plots/alignment_error_plot.pdf','--rscript_path',
         args.rscript_path]
  if args.alignment_error_scale:
    cmd += ['--scale']+[str(x) for x in args.alignment_error_scale]
  if args.alignment_error_max_length:
    cmd += ['--max_length',str(args.alignment_error_max_length)]
  if indfile:
    cmd += ['--input_index',indfile]
  cmd += ['--random']
  cmd += ['--specific_tempdir',args.tempdir+'/temp']
  sys.stderr.write("Making alignment error plot\n")
  sys.stderr.write(" ".join(cmd)+"\n")
  bam_to_alignment_error_plot.external_cmd(cmd)
  tlog.write(" ".join(cmd))
  tlog.stop()
  time.sleep(3)
  #gc.collect()
  return

def make_data_bam_annotation(args):
  global tlog
  udir = os.path.dirname(os.path.realpath(__file__))

  tlog.start("identify genomic features exon intron integenic")
  # 1. Use annotations to identify genomic features (Exon, Intron, Intergenic)
  # And assign membership to reads
  # Stores the feature bed files in a beds folder
  cmd = [udir+'/annotate_from_genomic_features.py','--output_beds',
         args.tempdir+'/data/beds',args.tempdir+'/data/best.sorted.gpd.gz',
         args.gpd,args.tempdir+'/data/chrlens.txt','-o',
         args.tempdir+'/data/read_genomic_features.txt.gz','--threads',
         str(args.threads)]
  sys.stderr.write("Finding genomic features and assigning reads membership\n")
  sys.stderr.write(" ".join(cmd)+"\n")
  tlog.write(" ".join(cmd))
  annotate_from_genomic_features.external_cmd(cmd)
  tlog.stop()

  time.sleep(3)
  tlog.start("get per-exon depth")
  # 2. Get depth distributions for each feature subset
  # now get depth subsets
  sys.stderr.write("get depths of features\n")
  cmd = [udir+'/get_depth_subset.py',args.tempdir+'/data/depth.sorted.bed.gz',
         args.tempdir+'/data/beds/exon.bed','-o',
         args.tempdir+'/data/exondepth.bed.gz']
  sys.stderr.write(" ".join(cmd)+"\n")
  get_depth_subset.external_cmd(cmd)
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("get per-intron subset")
  cmd = [udir+'/get_depth_subset.py',args.tempdir+'/data/depth.sorted.bed.gz',
         args.tempdir+'/data/beds/intron.bed','-o',
         args.tempdir+'/data/introndepth.bed.gz']
  sys.stderr.write(" ".join(cmd)+"\n")
  get_depth_subset.external_cmd(cmd)
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("get per intergenic depth")
  cmd = [udir+'/get_depth_subset.py',args.tempdir+'/data/depth.sorted.bed.gz',
         args.tempdir+'/data/beds/intergenic.bed','-o',
         args.tempdir+'/data/intergenicdepth.bed.gz']
  sys.stderr.write(" ".join(cmd)+"\n")
  get_depth_subset.external_cmd(cmd)
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("plot feature depth png")
  # 3. Plot the feature depth
  cmd = [args.rscript_path,udir+'/plot_feature_depth.r',
         args.tempdir+'/data/depth.sorted.bed.gz',
         args.tempdir+'/data/exondepth.bed.gz',
         args.tempdir+'/data/introndepth.bed.gz',
         args.tempdir+'/data/intergenicdepth.bed.gz',
         args.tempdir+'/plots/feature_depth.png']
  sys.stderr.write(" ".join(cmd)+"\n")
  mycall(cmd,args.tempdir+'/logs/featuredepth_png')
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("plot feature depth pdf")
  cmd = [args.rscript_path,udir+'/plot_feature_depth.r',
         args.tempdir+'/data/depth.sorted.bed.gz',
         args.tempdir+'/data/exondepth.bed.gz',
         args.tempdir+'/data/introndepth.bed.gz',
         args.tempdir+'/data/intergenicdepth.bed.gz',
         args.tempdir+'/plots/feature_depth.pdf']
  sys.stderr.write(" ".join(cmd)+"\n")
  mycall(cmd,args.tempdir+'/logs/featuredepth_pdf')
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("generate plots of which reads correspont to which features png")
  # 4. Generate plots from reads assigend to features
  sys.stderr.write("Plot read assignment to genomic features\n")
  cmd = [args.rscript_path,udir+'/plot_annotated_features.r',
         args.tempdir+'/data/read_genomic_features.txt.gz',
         args.tempdir+'/plots/read_genomic_features.png']
  sys.stderr.write(" ".join(cmd)+"\n")
  mycall(cmd,args.tempdir+'/logs/read_genomic_features_png')
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("generate plots of which reads correspont to which features pdf")
  cmd = [args.rscript_path,udir+'/plot_annotated_features.r',
         args.tempdir+'/data/read_genomic_features.txt.gz',
         args.tempdir+'/plots/read_genomic_features.pdf']
  sys.stderr.write(" ".join(cmd)+"\n")
  mycall(cmd,args.tempdir+'/logs/read_genomic_features_pdf')
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("annotate the reads")
  # 5. annotated the best mappend read mappings
  cmd = ['gpd_annotate.py',args.tempdir+'/data/best.sorted.gpd.gz','-r',
         args.gpd,'-o',args.tempdir+'/data/annotbest.txt.gz']
  if args.threads:
    cmd += ['--threads',str(args.threads)]
  sys.stderr.write("Annotating reads\n")
  sys.stderr.write(" ".join(cmd)+"\n")
  gpd_annotate(cmd)
  tlog.write(" ".join(cmd))
  tlog.stop()
  time.sleep(3)

  tlog.start("plot by transcript length png")
  # 6. Make plots of the transcript lengths
  sys.stderr.write("Make plots from transcript lengths\n")
  cmd = [args.rscript_path,udir+'/plot_transcript_lengths.r',
         args.tempdir+'/data/annotbest.txt.gz',
         args.tempdir+'/plots/transcript_distro.png']
  sys.stderr.write(" ".join(cmd)+"\n")
  mycall(cmd,args.tempdir+'/logs/transcript_distro_png')
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("plot by transcript length png")
  sys.stderr.write("Make plots from transcript lengths\n")
  cmd = [args.rscript_path,udir+'/plot_transcript_lengths.r',
         args.tempdir+'/data/annotbest.txt.gz',
         args.tempdir+'/plots/transcript_distro.pdf']
  sys.stderr.write(" ".join(cmd)+"\n")
  mycall(cmd,args.tempdir+'/logs/transcript_distro_pdf')
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("make length distribution from annotations")
  # 7. Make length distributions for plotting
  sys.stderr.write("making length distributions from annotations\n")
  cmd = [udir+'/annotated_length_analysis.py',
         args.tempdir+'/data/best.sorted.gpd.gz',
         args.tempdir+'/data/annotbest.txt.gz',
         '-o',args.tempdir+'/data/annot_lengths.txt.gz']
  sys.stderr.write(" ".join(cmd)+"\n")
  annotated_length_analysis.external_cmd(cmd)
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("plot annot length distro png")
  # 8. Plot length distributions
  cmd = [args.rscript_path,udir+'/plot_annotation_analysis.r',
         args.tempdir+'/data/annot_lengths.txt.gz',
         args.tempdir+'/plots/annot_lengths.png']
  sys.stderr.write(" ".join(cmd)+"\n")
  mycall(cmd,args.tempdir+'/logs/annot_lengths_png')
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("plot annot length distro png")
  cmd = [args.rscript_path,udir+'/plot_annotation_analysis.r',
         args.tempdir+'/data/annot_lengths.txt.gz',
         args.tempdir+'/plots/annot_lengths.pdf']
  sys.stderr.write(" ".join(cmd)+"\n")
  mycall(cmd,args.tempdir+'/logs/annot_lengths_pdf')
  tlog.write(" ".join(cmd))
  tlog.stop()

  # 9. Get rarefraction curve data
  global rcnt
  if rcnt < 0:
    sys.stderr.write("Getting read count\n")
    rcnt = 0
    tinf = gzip.open(args.tempdir+'/data/lengths.txt.gz')
    for line in tinf:  rcnt += 1
    tinf.close()

  tlog.start("get rarefraction gene")
  sys.stderr.write("Writing rarefraction curves\n")
  cmd =  [udir+'/gpd_annotation_to_rarefraction.py',
          args.tempdir+'/data/annotbest.txt.gz',
          '--samples_per_xval',str(args.samples_per_xval),
          '--original_read_count',str(rcnt),'--threads',str(args.threads),
          '--gene','-o',args.tempdir+'/data/gene_rarefraction.txt']
  sys.stderr.write(" ".join(cmd)+"\n")
  gpd_annotation_to_rarefraction.external_cmd(cmd)
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("rarefraction transcript")
  cmd =  [udir+'/gpd_annotation_to_rarefraction.py',
          args.tempdir+'/data/annotbest.txt.gz','--samples_per_xval',
          str(args.samples_per_xval),'--original_read_count',str(rcnt),
          '--threads',str(args.threads),'--transcript','-o',
          args.tempdir+'/data/transcript_rarefraction.txt']
  sys.stderr.write(" ".join(cmd)+"\n")
  gpd_annotation_to_rarefraction.external_cmd(cmd)
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("rarefraction gene full")
  cmd =  [udir+'/gpd_annotation_to_rarefraction.py',
          args.tempdir+'/data/annotbest.txt.gz','--samples_per_xval',
          str(args.samples_per_xval),'--original_read_count',str(rcnt),
          '--threads',str(args.threads),'--full','--gene','-o',
          args.tempdir+'/data/gene_full_rarefraction.txt']
  sys.stderr.write(" ".join(cmd)+"\n")
  gpd_annotation_to_rarefraction.external_cmd(cmd)
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("rarefraction gene full")
  cmd =  [udir+'/gpd_annotation_to_rarefraction.py',
          args.tempdir+'/data/annotbest.txt.gz','--samples_per_xval',
          str(args.samples_per_xval),'--original_read_count',str(rcnt),
          '--threads',str(args.threads),'--full','--transcript','-o',
          args.tempdir+'/data/transcript_full_rarefraction.txt']
  sys.stderr.write(" ".join(cmd)+"\n")
  gpd_annotation_to_rarefraction.external_cmd(cmd)
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("plot multiple rarefractions")
  # 10. Plot the rarefraction curves
  for type in ['gene','transcript']:
    for ext in ['png','pdf']:
      cmd = [args.rscript_path,udir+'/plot_annotation_rarefractions.r',
             args.tempdir+'/plots/'+type+'_rarefraction.'+ext,type,
             args.tempdir+'/data/'+type+'_rarefraction.txt','#FF000088',
             args.tempdir+'/data/'+type+'_full_rarefraction.txt','#0000FF88']
      sys.stderr.write(" ".join(cmd)+"\n")
      mycall(cmd,args.tempdir+'/logs/plot_'+type+'_rarefraction_'+ext)
      tlog.write(" ".join(cmd))
  tlog.stop()

  if os.name == 'nt' or sys.platform == 'darwin': return 
  ## For the bias data we need to downsample
  ## Using some system utilities to accomplish this
  sys.stderr.write("downsampling mappings for bias calculation\n")
  cmd0 = 'zcat'
  if args.threads > 1:
    cmd1 = ['sort','-R','-S1G','-T',
            args.tempdir+'/temp','--parallel='+str(args.threads)]
  else:
    cmd1 = ['sort','-R','-S1G','-T',args.tempdir+'/temp']
  cmd2 = 'head -n '+str(args.max_bias_data)
  if args.threads > 1:
    cmd3 = ['sort','-k3,3','-k5,5n','-k','6,6n','-S1G','-T',
            args.tempdir+'/temp --parallel='+str(args.threads)]
  else:
    cmd3 = ['sort','-k3,3','-k5,5n','-k','6,6n','-S1G','-T',
            args.tempdir+'/temp']
  inf = open(args.tempdir+'/data/best.sorted.gpd.gz')
  of = gzip.open(args.tempdir+'/temp/best.random.sorted.gpd.gz','w')
  if os.name != 'nt':
    p0 = Popen(cmd0.split(),stdin=inf,stdout=PIPE)
    p1 = Popen(cmd1,stdin=p0.stdout,stdout=PIPE)
    p2 = Popen(cmd2.split(),stdin=p1.stdout,stdout=PIPE)
    p3 = Popen(cmd3,stdin=p2.stdout,stdout=PIPE)
  else:
    sys.stderr.write("WARNING: Windows OS detected. using shell.")
    p0 = Popen(cmd0,stdin=inf,stdout=PIPE,shell=True)
    p1 = Popen(" ".join(cmd1),stdin=p0.stdout,stdout=PIPE,shell=True)
    p2 = Popen(cmd2,stdin=p1.stdout,stdout=PIPE,shell=True)
    p3 = Popen(" ".join(cmd3),stdin=p2.stdout,stdout=PIPE,shell=True)
  for line in p3.stdout:
    of.write(line)
  p3.communicate()
  p2.communicate()
  p1.communicate()
  p0.communicate()
  of.close()
  inf.close()
  # now downsample annotations
  sys.stderr.write("Downsampling annotations for bias\n")
  inf = gzip.open(args.tempdir+'/temp/best.random.sorted.gpd.gz')
  rnames = set()
  for line in inf:
    f = line.rstrip().split("\t")
    rnames.add(f[0])
  inf.close()
  of = gzip.open(args.tempdir+'/temp/annotbest.random.txt.gz','w')
  inf = gzip.open(args.tempdir+'/data/annotbest.txt.gz')
  for line in inf:
    f = line.rstrip().split("\t")
    if f[1] in rnames: of.write(line)
  inf.close()
  of.close()
  tlog.start("use annotations to check for 5' to 3' biase")
  # 11. Use annotation outputs to check for  bias
  sys.stderr.write("Prepare bias data\n")
  cmd = [udir+'/annotated_read_bias_analysis.py',
        args.tempdir+'/temp/best.random.sorted.gpd.gz',args.gpd,
        args.tempdir+'/temp/annotbest.random.txt.gz','-o',
        args.tempdir+'/data/bias_table.txt.gz','--output_counts',
        args.tempdir+'/data/bias_counts.txt','--allow_overflowed_matches',
        '--threads',str(args.threads),'--specific_tempdir',args.tempdir+'/temp']
  sys.stderr.write(" ".join(cmd)+"\n")
  annotated_read_bias_analysis.external_cmd(cmd)
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("plot bias png")
  # 12. Plot bias
  cmd = [args.rscript_path,udir+'/plot_bias.r',
         args.tempdir+'/data/bias_table.txt.gz',
         args.tempdir+'/plots/bias.png']
  sys.stderr.write(" ".join(cmd)+"\n")
  mycall(cmd,args.tempdir+'/logs/bias_png.log')
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("plot bias pdf")
  cmd = [args.rscript_path,udir+'/plot_bias.r',
         args.tempdir+'/data/bias_table.txt.gz',
         args.tempdir+'/plots/bias.pdf']
  sys.stderr.write(" ".join(cmd)+"\n")
  mycall(cmd,args.tempdir+'/logs/bias_pdf.log')
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("Prepare junction variance data")
  # 13. Get distances of observed junctions from reference junctions
  cmd = [udir+'/gpd_to_junction_variance.py','-r',args.gpd,
        args.tempdir+'/temp/best.random.sorted.gpd.gz',
        '--specific_tempdir',args.tempdir+'/temp','-o',
        args.tempdir+'/data/junvar.txt','--threads',str(args.threads)]
  sys.stderr.write(" ".join(cmd)+"\n")
  gpd_to_junction_variance.external_cmd(cmd)
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("plot junvar png")
  # 14. Junction distances
  cmd = [args.rscript_path,udir+'/plot_junvar.r',
         args.tempdir+'/data/junvar.txt',
         args.tempdir+'/plots/junvar.png']
  sys.stderr.write(" ".join(cmd)+"\n")
  mycall(cmd,args.tempdir+'/logs/junvar_png.log')
  tlog.write(" ".join(cmd))
  tlog.stop()

  tlog.start("plot junvar pdf")
  cmd = [args.rscript_path,udir+'/plot_junvar.r',
         args.tempdir+'/data/junvar.txt',
         args.tempdir+'/plots/junvar.pdf']
  sys.stderr.write(" ".join(cmd)+"\n")
  mycall(cmd,args.tempdir+'/logs/junvar_pdf.log')
  tlog.write(" ".join(cmd))
  tlog.stop()

  return


def mycall(cmd,lfile):
  ofe = open(lfile+'.err','w')
  ofo = open(lfile+'.out','w')
  p = Popen(cmd,stderr=ofe,stdout=ofo)
  p.communicate()
  ofe.close()
  ofo.close()
  return

#def do_inputs():
#  # Setup command line inputs
#  parser=argparse.ArgumentParser(description="Create an output report",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#  parser.add_argument('input',help="INPUT FILE or '-' for STDIN")
#  parser.add_argument('-o','--output',help="OUTPUT Folder or STDOUT if not set")
#  parser.add_argument('--portable_output',help="OUTPUT file in a portable html format")
#  group1 = parser.add_mutually_exclusive_group(required=True)
#  group1.add_argument('-r','--reference',help="Reference Fasta")
#  group1.add_argument('--no_reference',action='store_true',help="No Reference Fasta")
#  parser.add_argument('--annotation',help="Reference annotation genePred")
#  parser.add_argument('--threads',type=int,default=1,help="INT number of threads to run. Default is system cpu count")
#  # Temporary working directory step 1 of 3 - Definition
#  parser.add_argument('--tempdir',required=True,help="This temporary directory will be used, but will remain after executing.")
#
#  ### Parameters for alignment plots
#  parser.add_argument('--min_aligned_bases',type=int,default=50,help="for analysizing alignment, minimum bases to consider")
#  parser.add_argument('--max_query_overlap',type=int,default=10,help="for testing gapped alignment advantage")
#  parser.add_argument('--max_target_overlap',type=int,default=10,help="for testing gapped alignment advantage")
#  parser.add_argument('--max_query_gap',type=int,help="for testing gapped alignment advantge")
#  parser.add_argument('--max_target_gap',type=int,default=500000,help="for testing gapped alignment advantage")
#  parser.add_argument('--required_fractional_improvement',type=float,default=0.2,help="require gapped alignment to be this much better (in alignment length) than single alignment to consider it.")
#
#  ### Parameters for locus analysis
#  parser.add_argument('--do_loci',action='store_true',help="This analysis is time consuming at the moment so don't do it unless necessary")
#  parser.add_argument('--min_depth',type=float,default=1.5,help="require this or more read depth to consider locus")
#  parser.add_argument('--min_coverage_at_depth',type=float,default=0.8,help="require at leas this much of the read be covered at min_depth")
#  parser.add_argument('--min_exon_count',type=int,default=2,help="Require at least this many exons in a read to consider assignment to a locus")
#  parser.add_argument('--locus_downsample',type=int,default=100,help="Only include up to this many long reads in a locus\n")
#
#  ### Params for alignment error plot
#  parser.add_argument('--alignment_error_scale',nargs=6,type=float,help="<ins_min> <ins_max> <mismatch_min> <mismatch_max> <del_min> <del_max>")
#  parser.add_argument('--alignment_error_max_length',type=int,default=100000,help="The maximum number of alignment bases to calculate error from")
#
#  ### Params for context error plot
#  parser.add_argument('--context_error_scale',nargs=6,type=float,help="<ins_min> <ins_max> <mismatch_min> <mismatch_max> <del_min> <del_max>")
#  parser.add_argument('--context_error_stopping_point',type=int,default=1000,help="Sample at least this number of each context")
#
#  ## Params for rarefraction plots
#  parser.add_argument('--samples_per_xval',type=int,default=500)
#
#  args = parser.parse_args()
#  # Temporary working directory step 2 of 3 - Creation
#  setup_tempdir(args)
#  return args

def setup_tempdir(args):
  if not os.path.exists(args.tempdir):
    os.makedirs(args.tempdir.rstrip('/'))
  if not os.path.exists(args.tempdir):
    sys.stderr.write("ERROR: Problem creating temporary directory\n")
    sys.exit()
  return

def external(args):
  main(args)

if __name__=="__main__":
  sys.stderr.write("excute as prepare all data as main\n")
  #do our inputs
  # Can disable calling as main
  #args = do_inputs()
  #main(args)
