#!/usr/bin/python
import sys, re, os, inspect, argparse, gzip

#pre: genepred filename, an optional integer for smoothing size
#     The smoothing is done to combine adjacent "exons" within that size definition to get rid of small deletions or indels
#post: a bed file where each line has chromosome, index 1 start, and index 1 end coordinates and the gene_name from each psl line is listed

#bring in the folder to the path for our modules
#bring in the folder to the path for our modules

pythonfolder_loc = "../pythonlib"
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe() ))[0],pythonfolder_loc)))
if cmd_subfolder not in sys.path:
  sys.path.insert(0,cmd_subfolder)
from Bio.Format.GPD import GPDStream, GPD
#import GenePredBasics


def main(args):
  of = sys.stdout
  if args.output:
    if args.output[-3:]=='.gz':
      of = gzip.open(args.output,'w')
  color = '0,0,0'

  if args.color:
    if args.color == 'blue':
      color = '67,162,202'
    elif args.color == 'green':
      color = '49,163,84'
    elif args.color == 'orange':
      color = '254,178,76'
    elif args.color == 'purple':
      color = '136,86,167'
    elif args.color == 'red':
      color = '240,59,32'

  # set up the header if one is desired
  header = ''
  if not args.noheader:
    newname = 'longreads'
    m = re.search('([^\/]+)$',args.input)
    if m:
      newname = m.group(1)
    newname = re.sub('[\s]+','_',newname)
    if args.headername:
      newname = args.headername
    elif args.input == '-':
      newname = 'STDIN'
    header += "track\tname="+newname+"\t"
    description = newname+' GenePred Entries'
    if args.headerdescription:
       description = args.headerdescription
    header += 'description="'+description + '"'+"\t"
    header += 'itemRgb="On"'
    of.write(header+"\n")
  
  gpd_handle = sys.stdin
  if args.input != '-': 
    if args.input[-3:]=='.gz':
      gpd_handle = gzip.open(args.input)
    else:
      gpd_handle = open(args.input)
  gs = GPDStream(gpd_handle)
  #with gpd_handle as infile:
  for gpd in gs:
      #for line in infile:
      #if re.match('^#',line):
      #  continue
      #genepred_entry = GenePredBasics.line_to_entry(line)
      if args.minintron:
        gpd = GPD(gpd.smooth_gaps(args.minintron).get_gpd_line())
      exoncount = gpd.get_exon_count()
      ostring  = gpd.value('chrom') + "\t" 
      ostring += str(gpd.value('exonStarts')[0]) + "\t"
      ostring += str(gpd.value('exonEnds')[exoncount-1]) + "\t"
      if args.namefield == 1:
        ostring += gpd.value('gene_name') + "\t"
      else: 
        ostring += gpd.value('name')
      ostring += '1000' + "\t"
      ostring += gpd.value('strand') + "\t" 
      ostring += str(gpd.value('exonStarts')[0]) + "\t"
      ostring += str(gpd.value('exonEnds')[exoncount-1]) + "\t"      
      ostring += color+"\t"
      ostring += str(exoncount) + "\t"
      for i in range(0,exoncount):
        ostring += str(gpd.value('exonEnds')[i]-gpd.value('exonStarts')[i]) + ','
      ostring += "\t"
      for i in range(0,exoncount):
        ostring += str(gpd.value('exonStarts')[i]-gpd.value('exonStarts')[0])+','
      of.write(ostring+"\n")
      #for i in range(0,len(genepred_entry['exonStarts'])):
  gpd_handle.close()
  of.close()
def do_inputs():
  parser = argparse.ArgumentParser()
  parser.add_argument('--minintron',type=int,help='INT close gaps less than or equal to this, if set.')
  parser.add_argument('--namefield',type=int,default=1,choices=[1,2],help='INT[1 or 2] use the first or second field. Default (1)')
  parser.add_argument('--noheader',help='do not print any tack header')
  parser.add_argument('--headername',help='STRING name for the track. Default is the gpd file name')
  parser.add_argument('--headerdescription',help='STRING description for the track.')
  parser.add_argument('--color',choices=['blue','green','orange','purple','red'])
  parser.add_argument('input',help='FILENAME input genepred, use - for STDIN')
  parser.add_argument('-o','--output',help="output file or - for STDOUT")
  args = parser.parse_args()
  return args

def external_cmd(cmd):
  cache_argv = sys.argv
  sys.argv = cmd.split()
  args = do_inputs()
  main(args)
  sys.argv = cache_argv

if __name__=="__main__":
  args = do_inputs()
  main(args)
