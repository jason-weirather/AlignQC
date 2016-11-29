import os, re, gzip
import Bio.Sequence

#Iterable Stream
class FastaHandle:
  def __init__(self,fh,custom_buffer_size=10000000):
    self.fh = fh
    self.buffer_size = custom_buffer_size
    self.working_string = self.fh.read(self.buffer_size)
    self.buffered_results = []
    self.p = re.compile('>([^\n]+)\n([^>]+)')
    self.file_finished = False
    if not self.working_string: self.file_finished = True

  def __iter__(self):
    return self
  def next(self):
    v = self.get_entry()
    if not v:
      raise StopIteration
    else:
      return v

  def get_entry(self):
    if len(self.buffered_results) > 0:
      m = self.buffered_results.pop(0)
      m1 = re.match('(\S+)',m.group(1))
      return Bio.Sequence.Seq(m.group(2).rstrip(),m1.group(1))
    vals = [x for x in self.p.finditer(self.working_string)]
    while not self.file_finished:
      if vals: #have a match
        if vals[0].end() != len(self.working_string): # match is totally covered
          break
      chunk = self.fh.read(self.buffer_size)
      if not chunk: 
        self.file_finished = True
        break
      self.working_string += chunk
      vals = [x for x in self.p.finditer(self.working_string)]
    # If we've finished file clean up
    if self.file_finished:
      self.buffered_results = [x for x in self.p.finditer(self.working_string)]
      self.working_string = ''
    else:  # we have more to do so fix working string buffer
      self.buffered_results = [x for x in vals if x.end() != len(self.working_string)]
      self.working_string = self.working_string[self.buffered_results[-1].end():]
    if len(self.buffered_results) > 0:
      m = self.buffered_results.pop(0)
      m1 = re.match('(\S+)',m.group(1))
      return Bio.Sequence.Seq(m.group(2).rstrip(),m1.group(1))
    return None
# Slicable fast fasta
# It loses any additional header information in fasta header
# only the first non-whitespace is what we use
class FastaData:
  def __init__(self,data=None,file=None,dict=None):
    self._lengths = {}
    self._seqs = {}
    self._names = []
    # now populate
    if data:
      self._scan_data(data)
    elif file and re.search('\.gz$',file):
      self._scan_data(gzip.open(file,'rb').read())
    elif file:
      self._scan_data(open(file,'rb').read())
    elif dict:
      for name in dict.get_names(): 
        self._names.append(name)
        self._seqs = dict
        self._lengths[name] = len(dict[name])
  def clear(self):
    for k in self.keys():
      del self._seqs[k]

  def remove(self,key):
    del self._seqs[key]

  def keys(self):
    return self._seqs.keys()

  def __getitem__(self,key):
    return self._seqs[key]

  def get_sequence(self,chr=None,start=None,end=None,dir=None,rng=None):
    if rng: 
      chr = rng.chr
      start = rng.start
      end = rng.end
      dir = rng.direction
    if not start: start = 1
    if not end: end = self.fai[chr]['length']
    if not dir: dir = '+'
    if dir == '-':
      return Bio.Sequence.rc(self._seqs[chr][start-1:end])
    return self._seqs[chr][start-1:end]

  def _scan_data(self,dat):
    p = re.compile('>([^\n]+)\n([^>]+)')
    pos = 0
    for m in p.finditer(dat):
      m1 = re.match('(\S+)',m.group(1))
      self._names.append(m1.group(1))
      seq = m.group(2).replace("\n",'')
      self._lengths[m1.group(1)] = len(seq)
      self._seqs[m1.group(1)] = seq

# Do random access with an indexed Fasta File
# Creates the index if its not there already
# Pre: An uncompressed fasta file
#      Can be called by chromosome and location slices
#          Slices are same as array - zero indexed
# Post: Makes index if doesn't exist upon being called.
#       Can access sequence
# Modifies: File IO reads the fasta, and writes a fasta index file
class FastaFile:
  def __init__(self,fname,index=None):
    self.fname = fname
    self.index = index
    self.fai = {}
    if not self.index:
      if os.path.isfile(self.fname+'.fai'): self.index = self.fname+'.fai'
    if not self.index:
      sys.stderr.write("Warning no index trying to create\n")
      self._make_index()
    self._read_index()
    self.fh = open(fname)

  def __getitem__(self,key):
    chr = FastaFile.Chromosome(self,key)
    if chr.sliced: return chr
    return self.get_sequence(key)

  class Chromosome:
    def __init__(self,outer,chr):
      self.outer = outer
      self.chr = chr
      self.sliced = False

    def __getitem__(self,val):
      self.sliced = True
      if val.step:  
        sys.stderr.write("ERROR: FastaFile doesn't support step access\n")
        sys.exit()
      clen = self.outer.fai[self.chr]['length']
      return self.outer.get_sequence(self.chr,val.start+1,min(val.stop,clen))

    def __len__(self):
      return self.outer.fai[self.chr]['length']

  def get_sequence(self,chr=None,start=None,end=None,dir=None,rng=None):
    if rng: 
      chr = rng.chr
      start = rng.start
      end = rng.end
      dir = rng.direction
    if not start: start = 1
    if not end: end = self.fai[chr]['length']
    if not dir: dir = '+'
    [sblocks,srem] = divmod(start,self.fai[chr]['linebases'])
    missing_start = sblocks*(self.fai[chr]['linewidth']-self.fai[chr]['linebases'])
    [eblocks,erem] = divmod(end,self.fai[chr]['linebases'])
    missing_end = eblocks*(self.fai[chr]['linewidth']-self.fai[chr]['linebases'])
    pos_start = self.fai[chr]['offset']+start+missing_start-1
    pos_end = self.fai[chr]['offset']+end+missing_end
    self.fh.seek(pos_start)
    v = self.fh.read(pos_end-pos_start).replace("\n",'')
    if dir == '-':
      return Bio.Sequence.rc(v)
    return v

  def _read_index(self):
    with open(self.index) as inf:
      for line in inf:
        v = line.rstrip().split("\t")
        self.fai[v[0]] = {}
        self.fai[v[0]]['name'] = v[0]
        self.fai[v[0]]['length'] = int(v[1])
        self.fai[v[0]]['offset'] = int(v[2])
        self.fai[v[0]]['linebases'] = int(v[3])
        self.fai[v[0]]['linewidth'] = int(v[4])

  def _make_index(self):
    of = None
    try:
      of = open(self.fname+'.fai','w')
    except IOError:
      sys.stderr.write("ERROR: could not open file\n")
      sys.exit()
    self.index = self.fname+'.fai'
    # now we can parse the file
    p1 = re.compile('>([^\r\n]+)([\n\r]+)([^>]+)')
    p2 = re.compile('([^\r\n]+)($|[\n\r]+)')
    pos = 0
    for m1 in p1.finditer(open(self.fname).read()):
      name = m1.group(1)
      pos += len(m1.group(1))+1+len(m1.group(2))
      linewidth_bases = None
      linewidth_bytes = None
      seqlen = 0
      nextoffset = len(m1.group(3))
      for m2 in p2.finditer(m1.group(3)):
        if not linewidth_bases: linewidth_bases = len(m2.group(1))
        elif linewidth_bases != len(m2.group(1)) and len(m2.group(2)) > 0:
          sys.stderr.write("ERROR: irregular line breaks\n")
          sys.exit()
        if not linewidth_bytes: linewidth_bytes = len(m2.group(1))+len(m2.group(2))
        elif linewidth_bytes != len(m2.group(1))+len(m2.group(2)) and len(m2.group(2)) > 0:
          sys.stderr.write("ERROR: irregular bytes line breaks\n")
          sys.exit() 
        seqlen += len(m2.group(1))      
      of.write(name + "\t" + str(seqlen) + "\t"+str(pos)+"\t" + str(linewidth_bases) + "\t" + str(linewidth_bytes)+"\n")
      pos += nextoffset
    of.close()

