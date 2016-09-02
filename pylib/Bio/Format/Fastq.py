import re
import Bio.Sequence

#Iterable Stream
class FastqHandle:
  def __init__(self,fh):
    self.fh = fh
  def __iter__(self):
    return self
  def next(self):
    v = self.get_entry()
    if not v:
      raise StopIteration
    else:
      return v
  def get_entry(self):
    line1 = self.fh.readline().rstrip()[1:]
    if not line1: return None
    line2 = self.fh.readline().rstrip()
    if not line2: return None
    line3 = self.fh.readline().rstrip()
    if not line3: return None
    line4 = self.fh.readline().rstrip()
    if not line4: return None
    return Fastq([line1,line2,line3,line4])

class Fastq(Bio.Sequence.Seq):
  def __init__(self,v):
    self.lines = v
    self.name = v[0]
    self.seq = v[1]
    self.qual = v[3]
  def __getitem__(self,key):
    if isinstance(key,slice):
      newseq = self.seq[key.start:min(key.stop,len(self.seq))]
      newqual = self.qual[key.start:min(key.stop,len(self.seq))]
      return Fastq([self.name,newseq,self.lines[2],newqual])
    return {'name':self.name,'seq':self.seq,'qual':self.lines[3]}[key]
  def rc(self):
    return Fastq([self.name,Bio.Sequence.rc(self.seq),self.lines[2],self.qual[::-1]])
  def copy(self):
    return Fastq([self.name,self.seq,self.lines[2],self.qual])
  def fastq(self):
    return '@'+"\n".join(self.lines)+"\n"
  def __str__(self):
    return self.fastq().rstrip()
