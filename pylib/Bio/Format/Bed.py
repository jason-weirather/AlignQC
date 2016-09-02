import sys, uuid
import Bio.Structure
from Bio.Range import GenomicRange

# Bed format with 9 optional fields
class Bed12(Bio.Structure.Transcript):
  def __init__(self,bed_line):
    self._entry = self._line_to_entry(bed_line)
    self._line = bed_line.rstrip()
    self._range = None
    self.exons = []
    self.junctions = []
    self._payload = []
    self._direction = self.value('strand')
    self._transcript_name = None
    self._gene_name = None
    self._name = self.value('name')
    for i in range(0,self.value('blockCount')):
      ex = Bio.Structure.Exon(GenomicRange(self.value('chrom'),\
            self.value('chromStart')+self.value('blockStarts')[i]+1,\
            self.value('chromStart')+self.value('blockStarts')[i]+self.value('blockSizes')[i]))
      self.exons.append(ex)
    if self.value('blockCount') > 1:
      for i in range(0,self.value('blockCount')-1):
        l = GenomicRange(self.value('chrom'),\
             self.value('chromStart')+self.value('blockStarts')[i]+self.value('blockSizes')[i],\
             self.value('chromStart')+self.value('blockStarts')[i+1]+self.value('blockSizes')[i+1])
        r = GenomicRange(self.value('chrom'),\
             self.value('chromStart')+self.value('blockStarts')[i]+1,\
             self.value('chromStart')+self.value('blockStarts')[i+1]+1)
        junc = Bio.Structure.Junction(l,r)
        junc.set_exon_left(self.exons[i])
        junc.set_exon_right(self.exons[i+1])
        self.junctions.append(junc)
    self._range = GenomicRange(self.value('chrom'),self.value('chromStart')+1,self.value('chromEnd'))
    self._id = str(uuid.uuid4())
    self._sequence = None
  def __str__(self):
    return self.get_bed_line()  

  #output the original gpd line
  # Overrides Structure.Transcript
  def get_bed_line(self):
    return self._line

  def get_line(self):
    return self._line

  def value(self,key):
    return self._entry[key]

  def _line_to_entry(self,line):
    f = line.rstrip().split("\t")
    d = {}
    d['chrom'] = f[0]
    d['chromStart'] = int(f[1])
    d['chromEnd'] = int(f[2])
    d['name'] = f[3]
    d['score'] = int(f[4])
    d['strand'] = f[5]
    d['thickStart'] = int(f[6])
    d['thickEnd'] = int(f[7])
    d['itemRgb'] = [int(x) for x in f[8].rstrip(',').split(',')]
    d['blockCount'] = int(f[9])
    d['blockSizes'] = [int(x) for x in f[10].rstrip(',').split(',')]
    d['blockStarts'] = [int(x) for x in f[11].rstrip(',').split(',')]
    return d

