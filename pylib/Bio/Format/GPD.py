import uuid, sys, time
import Bio.Structure
from Bio.Range import GenomicRange
from subprocess import Popen, PIPE

# This whole format is a subclass of the Transcript subclass
class GPD(Bio.Structure.Transcript):
  def __init__(self,gpd_line):
    self._entry = self._line_to_entry(gpd_line)
    self._line = gpd_line.rstrip()
    self._range = None
    self.exons = []
    self.junctions = []
    self._payload = []
    self._direction = self.value('strand')
    self._gene_name = self.value('gene_name')
    self._transcript_name = self.value('name')
    self._name = None
    for i in range(0,self.value('exonCount')):
      ex = Bio.Structure.Exon(GenomicRange(self.value('chrom'),self.value('exonStarts')[i]+1,self.value('exonEnds')[i]))
      self.exons.append(ex)
    if self.value('exonCount') > 1:
      for i in range(0,self.value('exonCount')-1):
        l = GenomicRange(self.value('chrom'),self.value('exonEnds')[i],self.value('exonEnds')[i])
        r = GenomicRange(self.value('chrom'),self.value('exonStarts')[i+1]+1,self.value('exonStarts')[i+1]+1)
        junc = Bio.Structure.Junction(l,r)
        junc.set_exon_left(self.exons[i])
        junc.set_exon_right(self.exons[i+1])
        self.junctions.append(junc)
    self._range = GenomicRange(self.value('chrom'),self.value('exonStarts')[0]+1,self.value('exonEnds')[-1])
    self._id = str(uuid.uuid4())
    self._sequence = None
  def __str__(self):
    return self.get_gpd_line()  

  #output the original gpd line
  # Overrides Structure.Transcript
  def get_gpd_line(self):
    return self._line

  def get_line(self):
    return self._line

  def value(self,key):
    return self._entry[key]

  def _line_to_entry(self,line):
    f = line.rstrip().split("\t")
    d = {}
    d['gene_name'] = f[0]
    d['name'] = f[1]
    d['chrom'] = f[2]
    d['strand'] = f[3]
    d['txStart'] = int(f[4])
    d['txEnd'] = int(f[5])
    d['cdsStart'] = int(f[6])
    d['cdsEnd'] = int(f[7])
    d['exonCount'] = int(f[8])
    exonstarts = [int(x) for x in f[9].rstrip(",").split(",")]
    d['exonStarts'] = exonstarts
    exonends = [int(x) for x in f[10].rstrip(",").split(",")]
    d['exonEnds'] = exonends
    return d

class GPDStream:
  def __init__(self,fh):
    self.fh = fh

  def read_entry(self):
    ln = self.fh.readline()
    if not ln: return False
    gpd = GPD(ln)
    return gpd

  def __iter__(self):
    return self

  def next(self):
    r = self.read_entry()
    if not r:
      raise StopIteration
    else:
      return r

class SortedOutputFile:
  def __init__(self,filename,type='location',tempdir=None):
    if type not in ['location','name']:
      sys.stderr.write("ERROR: must be type location or name\n")
      sys.exit()
    self._gz = False
    self._fh = open(filename,'w')
    self._sh = None
    if filename[-3:] == '.gz':
      self._gz = True
    self._pipes  = []
    scmd = "sort -k1,1 -k2,2"
    if type == 'location':
      scmd = "sort -k3,3 -k5,5n -k6,6n -k4,4"
    if tempdir: scmd += " -T "+tempdir.rstrip('/')+'/'
    if self._gz:
      cmd1 = "gzip"
      p1 = Popen(cmd1.split(),stdout=self._fh,stdin=PIPE,close_fds=True)
      p2 = Popen(scmd.split(),stdout=p1.stdin,stdin=PIPE,close_fds=True)
      self._sh = p2.stdin
      self._pipes = [p2,p1]
    else:
      p = Popen(scmd.split(),stdout=self._fh,stdin=PIPE)
      self._sh = p.stdin
      self._pipes = [p]
  def write(self,value):
    self._sh.write(value)
  def close(self):
    #self._sh.flush()
    #self._sh.close()
    for p in self._pipes:
      #p.stdin.flush()
      #p.stdin.close()
      p.communicate()
    #self._pipes[0].stdin.flush()
    #self._pipes[0].stdin.close()
    #self._pipes[1].stdin.flush()
    #self._pipes[1].stdin.close()
    self._fh.close()
