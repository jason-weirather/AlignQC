# AlignQC
Generate a report on sequencing alignments to understand read alignments vs read sizes, error patterns in reads, annotations and rarefractions.

Ultimately this software should be suitable for assessing alignments of from a variety of sequencing platforms and a variety of sequence types.  The focus for the first version of this software is on the transcriptome analysis of third generation sequencing data outputs.

## Analysis
`alignqc analysis`

The analysis command is the most basic command for assessing an alignment.  It provides reports and plots in xhtml format.

### Inputs
#### 1. BAM format alignment file
The preferred format for transcriptome analysis is GMAP output (the 'samse') format.  Default output, or an output that can produce multiple alignment paths for a read is recommended if you want the ability to assess chimeric reads.

You can convert the SAM output of GMAP into BAM format using Samtools.
http://www.htslib.org/

Any properly formated BAM file should work however this software has only been tested with GMAP and hisat outputs at the moment.  
http://samtools.github.io/hts-specs/SAMv1.pdf

Please note that analyzing very large hiseq datasets has not been tested and memory requirements have not been optimized for this type of run.  If you want to check error patterns of HISEQ downsampling the data is advised.

#### (optional) 2. Genome fasta file
The reference genome these sequences were aligned to, in fasta format, can allows you to assess the error rates and error patterns in the alignments.

If you choose not to use a reference genome you must explicitly specify --no_reference

#### (optional) 3. GenePred format annotation file
Providing an annotation file provides context such as known transcripts, and exons, introns, and intergenic regions to help describe the data.  It is also necessary for rarefraction curves.

If you choose not to use a reference annotation you must explicitly specify --no_annotation

http://www.healthcare.uiowa.edu/labs/au/IDP/IDP_gpd_format.asp

The specific genepred format used here is "Gene Predictions and RefSeq Genes with Gene Names" genePred format described by UCSC.
https://genome.ucsc.edu/FAQ/FAQformat.html#format9

- geneName
- name
- chrom
- strand
- txStart
- txEnd
- cdsStart
- cdsEnd
- exoncount
- exonStarts
- exonEnds

### Outputs
At least one output format must be specified.

To view the output xhtml Mozilla Firefox or Google Chrome browser is recommend.

Since the recommended output type contains large URI data embedded in the xhtml page, Internet Explorer will likely not be compatible.  The memory requirements of the regular output may strain some systems.

If you want to share the visual results with others we recommend the `--portable_output` option because this version contains only the main text and png files.

If accessing the embedded data in the xhtml is a problem, you can output the data in a folder format `--output_folder`, which can provide you more convenient access.


#### Standard xhtml output
`-o` or `--output`

The recommended output of this software is a single xhtml file that contains all the relevant files from the analysis embeded as base64 encoded URI data.  

This means you can click any of the links in the document and the browser will download the data of interest (i.e. PDF format versions of a figure) from within the document.

Bed tracks compatible with the UCSC genome browser are also provided.

#### Portable xhtml output
`--portable_output`

This output is recommended if you want to email these results or share them over a bandwidth limited connection.  This format only has the png images and webpage text.  Links are disabled.

#### Output folder
`--output_folder`

Store an output folder that contains all the data and figures.  An xhtml file is still available in this folder.




