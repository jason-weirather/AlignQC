""" Get data regarding platform specific results

  input table is our length data file

  1. read name
  2. alignment type
  3. best aligned length
  4. total aligned length
  5. read length
"""
import sys, argparse, re, gzip, inspect, os

from seqtools.statistics import average

def main(args):
  #define read name programs

  #ONT basecalls
  #ont matches a uuid4
  ont_prog = re.compile('^[a-f0-9]{8}-[a-f0-9]{4}-4[a-f0-9]{3}-[89ab][a-f0-9]{3}-[a-f0-9]{12}_Basecall_2D_(.*)$')
  pacbio_prog = re.compile('^(m[^\/]+)\/(\d+)\/(ccs|\d+_\d+)')

  inf = None
  if re.search('\.gz$',args.input):
    inf = gzip.open(args.input)
  else:
    inf = open(args.input)
  unclassified = {'aligned':0,'unaligned':0}
  classified = {}
  pb_cell = {}
  #pb_mol = set()
  for line in inf:
    f = line.rstrip().split("\t")
    name = f[0]
    m_ont = ont_prog.match(name)
    m_pb = pacbio_prog.match(name)
    if m_ont:
      if 'ONT' not in classified:
        classified['ONT'] = {}
      type = m_ont.group(1)
      if type not in classified['ONT']:
        classified['ONT'][type] = {'aligned':0,'unaligned':0}
      if f[1] != 'unaligned':
        classified['ONT'][type]['aligned']+=1
      else:
        classified['ONT'][type]['unaligned']+=1
    elif m_pb:
      cell = m_pb.group(1)
      mol = int(m_pb.group(2))
      if cell not in pb_cell: pb_cell[cell] = {'molecules':set(),'reads':0,'molecules_aligned':set()}
      pb_cell[cell]['molecules'].add(mol)
      pb_cell[cell]['reads'] += 1
      #pb_mol.add(mol)
      if 'PacBio' not in classified:
        classified['PacBio'] = {}
      type = 'ccs'
      if m_pb.group(3) != 'ccs': type = 'subread'
      if type not in classified['PacBio']:
        classified['PacBio'][type] = {'aligned':0,'unaligned':0}
      if f[1] != 'unaligned':
        classified['PacBio'][type]['aligned']+=1
        pb_cell[cell]['molecules_aligned'].add(mol)
      else:
        classified['PacBio'][type]['unaligned']+=1      
    else:
      if f[1] != 'unaligned':
        unclassified['aligned']+=1
      else:
        unclassified['unaligned']+=1
  inf.close()
  # Finished reading the reads now we can make a report
  of = open(args.output_base,'w')  
  if len(classified.keys()) > 0: 
    of.write("SP\n")
  for classification in sorted(classified.keys()):
    for subclass in sorted(classified[classification].keys()):
      dat = classified[classification][subclass]
      of.write("GN\t"+classification+"\t"+subclass+"\t"+str(sum(dat.values()))+"\t"+str(dat['aligned'])+"\t"+str(dat['unaligned'])+"\n")
  of.write("GN\tUnclassified\t\t"+str(sum(unclassified.values()))+"\t"+str(unclassified['aligned'])+"\t"+str(unclassified['unaligned'])+"\n")
  if 'PacBio' in classified:
    of.write("PB\tCell Count\t"+str(len(pb_cell.keys()))+"\n")
    of.write("PB\tMolecule Count\t"+str(sum([len(pb_cell[x]['molecules']) for x in pb_cell.keys()]))+"\n")
    of.write("PB\tAligned Molecule Count\t"+str(sum([len(pb_cell[x]['molecules_aligned']) for x in pb_cell.keys()]))+"\n")
    of.write("PB\tMax Reads Per Cell\t"+str(max([pb_cell[x]['reads'] for x in pb_cell.keys()]))+"\n")
    of.write("PB\tAvg Reads Per Cell\t"+str(average([pb_cell[x]['reads'] for x in pb_cell.keys()]))+"\n")
    of.write("PB\tMin Reads Per Cell\t"+str(min([pb_cell[x]['reads'] for x in pb_cell.keys()]))+"\n")
    of.write("PB\tMax Molecules Per Cell\t"+str(max([len(pb_cell[x]['molecules']) for x in pb_cell.keys()]))+"\n")
    of.write("PB\tAvg Molecules Per Cell\t"+str(average([len(pb_cell[x]['molecules']) for x in pb_cell.keys()]))+"\n")
    of.write("PB\tMin Molecules Per Cell\t"+str(min([len(pb_cell[x]['molecules']) for x in pb_cell.keys()]))+"\n")
    of.write("PB\tMax Aligned Molecules Per Cell\t"+str(max([len(pb_cell[x]['molecules_aligned']) for x in pb_cell.keys()]))+"\n")
    of.write("PB\tAvg Aligned Molecules Per Cell\t"+str(average([len(pb_cell[x]['molecules_aligned']) for x in pb_cell.keys()]))+"\n")
    of.write("PB\tMin Aligned Molecules Per Cell\t"+str(min([len(pb_cell[x]['molecules_aligned']) for x in pb_cell.keys()]))+"\n")
    mols = [[len(pb_cell[x]['molecules_aligned']),len(pb_cell[x]['molecules'])] for x in pb_cell.keys()]
    smols = sorted(mols,key=lambda x: x[0])
    of.write("PB\tMolecules Per Cell Distro\t"+",".join(['/'.join([str(x[0]),str(x[1])]) for x in smols])+"\n")
    of1 = open(args.output_base+'.pacbio','w')
    for val in smols:
      of1.write(str(val[0])+"\t"+str(val[1])+"\n")
    of1.close()
  of.close()

def do_args():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="input lengths table")
  parser.add_argument('output_base',help="output file basename")
  args = parser.parse_args()
  return args  

def external_cmd(cmd):
  cache_argv = sys.argv
  sys.argv = cmd
  args = do_args()
  main(args)
  sys.argv = cache_argv

if __name__=="__main__":
  args = do_args()
  main(args)
