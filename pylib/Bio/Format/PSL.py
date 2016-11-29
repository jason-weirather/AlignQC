import Bio.Align
from Bio.Range import GenomicRange
from Bio.Sequence import rc

class PSL(Bio.Align.Alignment):
  #Pre: psl_line is a psl formated line
  #     reference is a dict/slice accessable sequence
  #     query_sequences is a dict/slice accessable sequence
  #     query_sequence is just the string that is the query
  def __init__(self,psl_line,reference=None,query_sequences=None,query_sequence=None,query_quality=None):
    self._line = psl_line.rstrip()
    self._query_sequences = query_sequences
    self._query_sequence = query_sequence
    self._query_quality = query_quality
    self._reference = reference
    self._target_range = None
    self._private_values = PSL.PrivateValues()
    self._parse_psl_line()
    # Private values holds entries
    self._alignment_ranges = None
    self._set_alignment_ranges()

  def __str__(self):
    return self._line

  def get_line(self):
    return self._line

  #Do our overrides of Bio.Alignment.Align functions
  #Overrides Bio.Alignment.Align.get_query_sequence()
  def get_query_sequence(self):
    if self._query_sequence: return self._query_sequence
    if not self._query_sequences: return None
    if self.value('qName') not in self._query_sequences: return None
    #if self.value('strand') == '-': return rc(self._query_sequences[self.value('qName')])
    return self._query_sequences[self.value('qName')]

  def get_query_quality(self):
    return self._query_quality

  #Overrides Bio.Alignment.Align.get_reference()
  def get_reference(self):
    return self._reference
  #Overrides Bio.Alignment.Align.get_query_length()
  def get_query_length(self):
    return self.value('qSize')
  # Override the obvious access
  def get_PSL(self):
    return self
  # Override 
  def get_strand(self):
    return self.value('strand')

  def _parse_psl_line(self):
    f = self._line.rstrip().split("\t")
    if len(f) != 21:
      sys.stderr.write("ERROR: PSL line must contain 21 entries\n")
    self._private_values.set_entry('matches',int(f[0]))
    self._private_values.set_entry('misMatches',int(f[1]))
    self._private_values.set_entry('repMatches',int(f[2]))
    self._private_values.set_entry('nCount',int(f[3]))
    self._private_values.set_entry('qNumInsert',int(f[4]))
    self._private_values.set_entry('qBaseInsert',int(f[5]))
    self._private_values.set_entry('tNumInsert',int(f[6]))
    self._private_values.set_entry('tBaseInsert',f[7])
    self._private_values.set_entry('strand',f[8])
    self._private_values.set_entry('qName',f[9])
    self._private_values.set_entry('qSize',int(f[10]))
    self._private_values.set_entry('qStart',int(f[11]))
    self._private_values.set_entry('qEnd',int(f[12]))
    self._private_values.set_entry('tName',f[13])
    self._private_values.set_entry('tSize',int(f[14]))
    self._private_values.set_entry('tStart',int(f[15]))
    self._private_values.set_entry('tEnd',int(f[16]))
    self._private_values.set_entry('blockCount',int(f[17]))
    self._private_values.set_entry('blockSizes',[int(x) for x in f[18].rstrip(',').split(',')])
    self._private_values.set_entry('qStarts',[int(x) for x in f[19].rstrip(',').split(',')])
    self._private_values.set_entry('tStarts',[int(x) for x in f[20].rstrip(',').split(',')])
    # Now we can set things
  # Set list of [target range, query range]
  def _set_alignment_ranges(self):
    self._target_range = GenomicRange(self.value('tName'),self.value('tStart'),self.value('tEnd'))
    self._alignment_ranges = []
    for i in range(0,len(self.value('blockSizes'))):
      trng = GenomicRange(self.value('tName'),self.value('tStarts')[i]+1,self.value('tStarts')[i]+self.value('blockSizes')[i])
      qrng = GenomicRange(self.value('qName'),self.value('qStarts')[i]+1,self.value('qStarts')[i]+self.value('blockSizes')[i])
      self._alignment_ranges.append([trng,qrng])
    return

  #Here is how we access value
  def value(self,key):
    return self._private_values.get_entry(key)    
  # Values from the orginal should just be accessed though functions for consistancy sake.
  # This class should remind us well that entires need to be accessed this way
  class PrivateValues:
    def __init__(self):
      self.__entries = {}
    def set_entries_dict(self,mydict): self.__entries = mydict
    def get_entry(self,key):
      if key not in self.__entries:
        sys.stderr.write("WARNING: key "+str(key)+"not in entries\n")
        return None
      return self.__entries[key]
    def is_entry_key(self,key):
      if key in self.__entries:  return True
      return False
    def set_entry(self,key,value): self.__entries[key] = value
