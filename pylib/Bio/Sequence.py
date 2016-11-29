import re, sys, base64, zlib
from string import maketrans

#Basic Sequence structure
class Seq:
  def __init__(self,seq=None,name=None):
    self.name = name
    self.seq = seq
  def __getitem__(self,key):
    #if its a slice deal with the sequence only
    if isinstance(key,slice):
      newseq = self.seq[key.start:min(key.stop,len(self.seq))]
      return Seq(newseq,self.name)
    return {'name':self.name,'seq':self.seq}[key]
  def __str__(self):
    return self.seq
  def __len__(self):
    return len(self.seq)
  def rc(self):
    return Seq(rc(self.seq),self.name)
  def copy(self):
    return Seq(self.seq,self.name)


  #Pre: seq and name are set
  #Post: string representation of the fasta entry
  def fasta(self):
    return '>'+self.name+"\n"+self.seq+"\n"

  #Pre: seq is set
  #Post: float with the gc fraction
  def gc_content(self):
    if len(self.seq) == 0: return None
    n_count = self.n_count()
    if len(self.seq) - n_count == 0: return None
    return float(self.seq.translate(maketrans('GCgc','GGGG')).count('G')-n_count)/float(len(self.seq)-n_count)
  def n_count(self):
    return self.seq.translate(maketrans('Nn','NN')).count('N')

def rc(seq):
  complement = maketrans('ACTGUNXactgunx','TGACANXtgacanx')
  return seq.translate(complement)[::-1]


def encode_name(conversion_string):
  compressed_string = zlib.compress(conversion_string,9)
  enc_string = base64.b32encode(compressed_string)
  return 'SZ_'+enc_string.rstrip('=')

def decode_name(safename):
  frag = safename.lstrip('SZ_')
  padding = (8-(len(frag) % 8)) % 8
  c = base64.b32decode(frag+'='*padding)
  return  zlib.decompress(c)

