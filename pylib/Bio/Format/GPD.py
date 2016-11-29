import uuid, sys, time, re
import Bio.Structure
from Bio.Range import GenomicRange
from subprocess import Popen, PIPE

# This whole format is a subclass of the Transcript subclass
class GPD(Bio.Structure.Transcript):
  def __init__(self,gpd_line):
    # Only store the line and ID at first.  
    self._line = gpd_line.rstrip()
    self._id = str(uuid.uuid4())
    m = re.match('[^\t]+\t[^\t]+\t([^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)',gpd_line)
    self._range = GenomicRange(m.group(1),int(m.group(2))+1,int(m.group(3)))
    self._initialized = False
    # Most of GPD has not been set yet.  Each method accessing GPD
    # will need to check to see if initialize has been run

  def _initialize(self): # Wait to initialize to speed up streaming
    if self._initialized: return # nothing to do if its done
    self._initialized = True
    self._entry = _line_to_entry(self._line)
    self._exons = []
    self._junctions = []
    self._payload = []
    self._direction = self.value('strand')
    self._gene_name = self.value('gene_name')
    self._transcript_name = self.value('name')
    self._name = None
    for i in range(0,self.value('exonCount')):
      ex = Bio.Structure.Exon(GenomicRange(self.value('chrom'),self.value('exonStarts')[i]+1,self.value('exonEnds')[i]))
      self._exons.append(ex)
    if self.value('exonCount') > 1:
      for i in range(0,self.value('exonCount')-1):
        l = GenomicRange(self.value('chrom'),self.value('exonEnds')[i],self.value('exonEnds')[i])
        r = GenomicRange(self.value('chrom'),self.value('exonStarts')[i+1]+1,self.value('exonStarts')[i+1]+1)
        junc = Bio.Structure.Junction(l,r)
        junc.set_exon_left(self._exons[i])
        junc.set_exon_right(self._exons[i+1])
        self._junctions.append(junc)
    self._sequence = None

  @property
  def junctions(self):
    self._initialize()
    return self._junctions
  @property
  def exons(self):
    self._initialize()
    return self._exons

  # override, we are garunteed to have the range since we initialize on reading a line
  def get_range(self):
    return self._range

  def __str__(self):
    return self.get_gpd_line()  

  #output the original gpd line
  # Overrides Structure.Transcript
  def get_gpd_line(self):
    return self._line

  def get_line(self):
    return self._line

  def value(self,key):
    self._initialize()
    return self._entry[key]

def _line_to_entry(line):
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
