import struct, zlib, sys, re, os, gzip, random
import Bio.Align
#import Bio.Format.BamIndex as BamIndex
from Bio.Sequence import rc
from cStringIO import StringIO
from string import maketrans
from Bio.Range import GenomicRange
from subprocess import Popen, PIPE
_bam_ops = maketrans('012345678','MIDNSHP=X')
_bam_char = maketrans('abcdefghijklmnop','=ACMGRSVTWYHKDBN')
_bam_value_type = {'c':[1,'<b'],'C':[1,'<B'],'s':[2,'<h'],'S':[2,'<H'],'i':[4,'<i'],'I':[4,'<I']}
_sam_cigar_target_add = re.compile('[M=XDN]$')

# A sam entry
class SAM(Bio.Align.Alignment):
  def __init__(self,line,reference=None,reference_lengths=None):
    self._line = line.rstrip()
    self._reference = reference
    self._reference_lengths = None # reference would also cover this
    self._target_range = None
    self._private_values = SAM.PrivateValues()
    self._parse_sam_line()
    # Private values holds tags cigar and entries
    self._alignment_ranges = None
    self._set_alignment_ranges()
    return

  def __str__(self):
    if self._line:
      return self._line
    return self.get_line()

  # Get the length of the target sequence
  def get_target_length(self):
    if not self.is_aligned():
      sys.stderr.write("ERROR no length for reference when not aligned\n")
      sys.exit()
    if self._reference_lengths:
      if self.value('rname') in self._reference_lengths:
        return self._reference_lengths[self.value('rname')]
    elif self._reference:
      return len(self._reference[self.value('rname')])
    else:
      sys.stderr.write("ERROR some reference needs to be set to go from psl to bam\n")
      sys.exit()
    sys.stderr.write("ERROR reference found\n")
    sys.exit()
 
  #Overrides Bio.Alignment.Align.get_query_sequence()
  def get_query_sequence(self):
    if self.value('seq') == '*': return None
    if self.check_flag(0x10): return rc(self.value('seq'))
    return self.value('seq')
  #Overrides Bio.Alignment.Align.get_query_sequence()
  def get_query_quality(self):
    if not self.get_query_sequence(): return None
    if self.value('qual') == '*': return None
    if self.check_flag(0x10): return self.value('qual')[::-1]
    return self.value('qual')

  #Overrides Bio.Alignment.Align.get_query_length()
  def get_query_length(self):
    seq = self.value('seq')
    if seq != '*': return len(self.value('seq'))
    return sum([x[0] for x in self.get_cigar() if re.match('[MIS=X]',x[1])])

  # Similar to get_get_query_length, but it also includes
  # hard clipped bases
  # if there is no cigar, then default to trying the sequence
  def get_original_query_length(self):
    if not self.is_aligned():
      return self.get_query_length()    
    if self.get_cigar() == '*':
      return self.get_query_length()
    return sum([x[0] for x in self.get_cigar() if re.match('[HMIS=X]',x[1])])
  
  # This accounts for hard clipped bases 
  # and a query sequence that hasnt been reverse complemented
  def get_actual_original_query_range(self):
    l = self.get_original_query_length()
    a = self.get_alignment_ranges()
    qname = a[0][1].chr
    qstart = a[0][1].start
    qend = a[-1][1].end
    #rng = self.get_query_range()
    start = qstart
    end = qend
    if self.get_strand() == '-':
      end = l-(qstart-1)
      start = 1+l-(qend)
    return GenomicRange(qname,start,end,self.get_strand())

  #Overrides Bio.Alignment.Align.get_strand()
  #Which strand is the query aligned to
  def get_strand(self):
    if self.check_flag(0x10): return '-'
    return '+'

  #Overrides Bio.Alignment.Align.get_SAM()
  def get_SAM(self):
    return self

  def get_tag(self,key):
    return self._private_values.get_tags()[key]['value']

  #Overrides Bio.Alignment.Align._set_alignment_ranges()
  #[target, query]
  def _set_alignment_ranges(self):
    if not self.is_aligned(): 
      self._alignment_ranges = None
      return
    self._alignment_ranges = []
    cig = self.get_cigar()[:]
    target_pos = self.value('pos')
    query_pos = 1
    while len(cig) > 0:
      c = cig.pop(0)
      if re.match('[S]$',c[1]): # hard or soft clipping
        query_pos += c[0]
      elif re.match('[ND]$',c[1]): # deleted from reference
        target_pos += c[0]
      elif re.match('[I]$',c[1]): # insertion to the reference
        query_pos += c[0]
      elif re.match('[MI=X]$',c[1]): # keep it
        t_start = target_pos
        q_start = query_pos
        target_pos += c[0]
        query_pos += c[0]
        t_end = target_pos-1
        q_end = query_pos-1
        self._alignment_ranges.append([GenomicRange(self.value('rname'),t_start,t_end),GenomicRange(self.value('qname'),q_start,q_end)])
    return

  def _parse_sam_line(self):
    f = self._line.rstrip().split("\t")
    self._private_values.set_entry('qname',f[0])
    self._private_values.set_entry('flag',int(f[1]))
    self._private_values.set_entry('rname',f[2])
    if f[2] == '*':
      self._private_values.set_entry('pos',0)
    else: 
      self._private_values.set_entry('pos',int(f[3]))
    self._private_values.set_entry('mapq',int(f[4]))
    self._private_values.set_entry('cigar',f[5])
    self._private_values.set_entry('rnext',f[6])
    self._private_values.set_entry('pnext',int(f[7]))
    self._private_values.set_entry('tlen',int(f[8]))
    self._private_values.set_entry('seq',f[9])
    self._private_values.set_entry('qual',f[10])
    self._private_values.set_cigar([])
    if self.value('cigar') != '*':
      cig = [[int(m[0]),m[1]] for m in re.findall('([0-9]+)([MIDNSHP=X]+)',self.value('cigar'))]
      self._private_values.set_cigar(cig)
    tags = {}
    if len(f) > 11:
      for m in [[y.group(1),y.group(2),y.group(3)] for y in [re.match('([^:]{2,2}):([^:]):(.+)$',x) for x in f[11:]]]:
        if m[1] == 'i': m[2] = int(m[2])
        elif m[1] == 'f': m[2] = float(m[2])
        tags[m[0]] = {'type':m[1],'value':m[2]}
    self._private_values.set_tags(tags)

  # Necessary function for doing a locus stream
  # For the context of a SAM file we set this to be the target range
  def get_range(self):
    return self.get_target_range()

  def get_target_range(self):
    if not self.is_aligned(): return None
    if self._target_range: return self._target_range
    global _sam_cigar_target_add
    tlen = sum([x[0] for x in self.get_cigar() if _sam_cigar_target_add.match(x[1])])
    self._target_range = GenomicRange(self.value('rname'),self.value('pos'),self.value('pos')+tlen-1)
    return self._target_range
  def check_flag(self,inbit):
    if self.value('flag') & inbit: return True
    return False
  def is_aligned(self):
    return not self.check_flag(0x4)

  #assemble the line if its not there yet
  def get_line(self):
    if not self._line:
      chr = self.value('rname')
      rnext = self.value('rnext')
      if not self.is_aligned(): 
        chr = '*'
        rnext = '*'
      self._line = self.value('qname')+"\t"+str(self.value('flag'))+"\t"+chr+"\t"+str(self.value('pos'))+"\t"+str(self.value('mapq'))+"\t"+self.value('cigar')+"\t"+rnext+"\t"+str(self.value('pnext'))+"\t"+str(self.value('tlen'))+"\t"+self.value('seq')+"\t"+self.value('qual')
      if self.value('remainder'):
        self._line += "\t"+self.value('remainder')
    return self._line
  def value(self,key):
    return self._private_values.get_entry(key)
  def get_tags(self): 
    return self._private_values.get_tags()
  def get_cigar(self): 
    return self._private_values.get_cigar()

  #Bam files need a specific override to get_tags and get_cigar that would break other parts of the class if we 
  # access the variables other ways
  #Force tags and cigars to be hidden so we don't accidently change them.
  class PrivateValues:
    def __init__(self):
      self.__tags = None
      self.__cigar = None
      self.__entries = {}
    def set_tags(self,tags): self.__tags=tags
    def get_tags(self): return self.__tags
    def set_cigar(self,cigar): self.__cigar=cigar
    def get_cigar(self): return self.__cigar
    def set_entries_dict(self,mydict): self.__entries = mydict # set the entire dictionary at once
    def get_entry(self,key): 
      if key not in self.__entries:
        sys.stderr.write("WARNING: key "+str(key)+"not in entries\n")
        return None
      return self.__entries[key]
    def is_entry_key(self,key):  
      if key in self.__entries:  return True
      return False
    def set_entry(self,key,value): self.__entries[key] = value

# Very much like a sam entry but optimized for access from a bam
# Slows down for accessing things that need more decoding like
# sequence, quality, cigar string, and tags
class BAM(SAM):
  def __init__(self,bin_data,ref_names,fileName=None,blockStart=None,innerStart=None,ref_lengths=None,reference=None,line_number=None):
    part_dict = _parse_bam_data_block(bin_data,ref_names)
    #self._bamfileobj = bamfileobj #this is most like our parent
    self._line = None
    self._line_number = line_number # the line number in the bam file
    self._reference = reference
    self._target_range = None
    self._alignment_ranges = None #should be accessed by method because of BAM
    self._ref_lengths = ref_lengths
    self._file_position = {'fileName':fileName,'blockStart':blockStart,'innerStart':innerStart} # The most special information about the bam
    self._private_values = BAM.PrivateValues() # keep from accidently accessing some variables other than by methods
    self._private_values.set_entries_dict(part_dict)
    #self._set_alignment_ranges()
    return

  def get_alignment_ranges(self):
    if not self._alignment_ranges:
      self._set_alignment_ranges()
    return self._alignment_ranges

  def get_line_number(self):
    return self._line_number
  def get_target_length(self):
    return self._ref_lengths[self.value('rname')]
  def get_filename(self):
    return self._file_position['fileName']
  def get_coord(self):
    return [self._file_position['blockStart'],self._file_position['innerStart']]
  def get_block_start(self):
    return self._file_position['blockStart']
  def get_inner_start(self):
    return self._file_position['innerStart']
  def get_file_position_string(self):
    return 'fileName: '+self._file_position['fileName']+" "\
           'blockStart: '+str(self._file_position['blockStart'])+" "\
           'innerStart: '+str(self._file_position['innerStart'])
  def get_tag(self,key): 
    cur = self._private_values.get_tags()
    if not cur:
      v1,v2 = _bin_to_extra(self.value('extra_bytes'))
      self._private_values.set_tags(v1) #keep the cigar array in a special palce
      self._private_values.set_entry('remainder',v2)
    return self._private_values.get_tags()[key]['value']
  def get_cigar(self): 
    cur = self._private_values.get_cigar()
    if not cur:
      v1,v2 = _bin_to_cigar(self.value('cigar_bytes'))
      self._private_values.set_cigar(v1) #keep the cigar array in a special palce
      self._private_values.set_entry('cigar',v2)
    return self._private_values.get_cigar()
  #def indexed_as_primary_alignment(self):
  #  ind = self._bamfileobj.index
  #  if not ind:
  #    sys.stderr.write("ERROR: to access indexed as primary alignment, and index must have been loaded\n")
  #    sys.exit()
  #  e = self._bamfileobj.index.get_index_line(self._line_number)
  #  if e['flag'] & 2304: return False
  #  return True
  def value(self,key):
    if not self._private_values.is_entry_key(key):
      if key == 'seq':
        v = _bin_to_seq(self.value('seq_bytes'))
        if not v: v = '*'
        self._private_values.set_entry('seq',v)
        return v
      elif key == 'qual':
        v = _bin_to_qual(self.value('qual_bytes'))
        if not v: v = '*'
        self._private_values.set_entry('qual',v)
        return v
      elif key == 'cigar':
        v1,v2 = _bin_to_cigar(self.value('cigar_bytes'))
        self._private_values.set_cigar(v1) #keep the cigar array in a special palce
        self._private_values.set_entry('cigar',v2)
        return v2
      elif key == 'remainder':
        v1,v2 = _bin_to_extra(self.value('extra_bytes'))
        self._private_values.set_tags(v1) #keep the cigar array in a special palce
        self._private_values.set_entry('remainder',v2)
        return v2
    return self._private_values.get_entry(key)


class SAMHeader:
  def __init__(self,header_text):
    self._text = header_text
    self.tags = []
    self._sequence_lengths = {}
    for line in self._text.split("\n"):
      if len(line) == 0: continue
      tag = line[1:3]
      rem = line[4:]
      self.tags.append({'tag':tag,'info':{}})
      for c in [{'field':x[0:2],'value':x[3:]} for x in rem.split("\t")]:
        self.tags[-1]['info'][c['field']] = c['value']
    for v in [x['info'] for x in self.tags if x['tag'] == 'SQ']:
      self._sequence_lengths[v['SN']] = int(v['LN'])
    return
  def get_sequence_names(self):
    return self._sequence_lengths.keys()
  #dictionary to get sequence lengths
  def get_sequence_lengths(self):
    return self._sequence_lengths
  def get_sequence_length(self,sname):
    return self._sequence_lengths[sname]

# reference is a dict
class BAMFile:
  #def __init__(self,filename,blockStart=None,innerStart=None,cnt=None,index_obj=None,index_file=None,reference=None):
  def __init__(self,filename,blockStart=None,innerStart=None,cnt=None,reference=None):
    self.path = filename
    self._reference = reference # dict style accessable reference
    self.fh = BGZF(filename)
    self._line_number = 0 # entry line number ... after header.  starts with 1
    # start reading the bam file
    self.header_text = None
    self._header = None
    self.n_ref = None
    self._read_top_header()
    self.ref_names = []
    self.ref_lengths = {}
    self._output_range = None
    #self.index = index_obj
    self._read_reference_information()
    # prepare for specific work
    if self.path and blockStart is not None and innerStart is not None:
      self.fh.seek(blockStart,innerStart)
      #if self.index:
      #  lnum = self.index.get_coord_line_number([blockStart,innerStart])
      #  if lnum:
      #    self._line_number = lnum-1 #make it zero indexed
  def close(self):
    self.fh.close()
  # return a string that is the header
  def get_header(self):
    if not self._header:
      self._header = SAMHeader(self.header_text)
      return self._header
    return self._header
  #def has_index(self):
  #  if self.index: return True
  #  return False

  # Index file is a gzipped TSV file with these fields:
  # 1. qname
  # 2. target range
  # 3. bgzf file block start
  # 4. bgzf inner block start
  # 5. aligned base count
  # 6. flag
  #def write_index(self,index_file,verbose=False):
  #  BamIndex.write_index(self.path,index_file,verbose=verbose)
  #  #_write_index(self.path,index_file,verbose=verbose)

  #def read_index(self,index_file=None):
  #  #prepare index
  #  if index_file: 
  #    self.index = BamIndex.BAMIndex(index_file)
  #    return True
  #  elif os.path.exists(self.path+'.bgi'):
  #    self.index = BamIndex.BAMIndex(self.path+'.bgi')
  #    return True
  #  return False

  def __iter__(self):
    return self

  def read_entry(self):
    e = self.read_entry2()
    #print e
    if self._output_range: # check and see if we are past out put range
      if not e.is_aligned(): 
        e = None
      else:
        rng2 = e.get_target_range()
        if self._output_range.chr != rng2.chr: e = None 
        if self._output_range.cmp(rng2) == 1: e = None
    if not e:
      return None
    else: return e

  def next(self):
    e = self.read_entry()
    if not e:
      raise StopIteration
    else: return e

  def read_entry2(self):
    bstart = self.fh.get_block_start()
    innerstart = self.fh.get_inner_start()
    b = self.fh.read(4) # get block size bytes
    if not b: return None
    block_size = struct.unpack('<i',b)[0]
    #print 'block_size '+str(block_size)
    self._line_number += 1
    bam = BAM(self.fh.read(block_size),self.ref_names,fileName=self.path,blockStart=bstart,innerStart=innerstart,ref_lengths=self.ref_lengths,reference=self._reference,line_number = self._line_number)
    return bam

  def _set_output_range(self,rng):
    self._output_range = rng
    return

  #def fetch_random(self):
  #  cnt = self.index.get_length()
  #  num = random.randint(0,cnt-1)
  #  iline = self.index.get_index_line(num+1)
  #  bf2 = BAMFile(self.path,blockStart=iline['filestart'],innerStart=iline['innerstart'])
  #  bam = bf2.read_entry()
  #  bf2.close()
  #  return bam

  #def fetch_by_range(self,rng):
  #  coord = self.index.get_range_start_coord(rng)
  #  line_number = self.index.get_range_start_line_number(rng)
  #  if not coord: return None
  #  b2 = BAMFile(self.path,blockStart=coord[0],innerStart=coord[1],index_obj=self.index,reference=self._reference)
  #  b2._set_output_range(rng)
  #  return b2

  # A special way to access via bam
  #def fetch_by_query(self,name):
  #  bams = []
  #  for coord in self.index.get_coords_by_name(name):
  #    b2 = BAMFile(self.path,blockStart=coord[0],innerStart=coord[1],index_obj=self.index,reference=self._reference)
  #    bams.append(b2.read_entry())
  #    b2.close()
  #  return bams
    
  # only get a single
  def fetch_by_coord(self,coord):
    #b2 = BAMFile(self.path,blockStart=coord[0],innerStart=coord[1],index_obj=self.index,reference=self._reference)
    b2 = BAMFile(self.path,blockStart=coord[0],innerStart=coord[1],reference=self._reference)
    bam = b2.read_entry()
    b2.close()
    b2 = None
    return bam

  def fetch_starting_at_coord(self,coord):
    #b2 = BAMFile(self.path,blockStart=coord[0],innerStart=coord[1],index_obj=self.index,reference=self._reference)
    b2 = BAMFile(self.path,blockStart=coord[0],innerStart=coord[1],reference=self._reference)
    return b2

  def _read_reference_information(self):
    for n in range(self.n_ref):
      l_name = struct.unpack('<i',self.fh.read(4))[0]
      name = self.fh.read(l_name).rstrip('\0')
      l_ref = struct.unpack('<i',self.fh.read(4))[0]
      self.ref_lengths[name] = l_ref
      self.ref_names.append(name)
  def _read_top_header(self):
    magic = self.fh.read(4)
    l_text = struct.unpack('<i',self.fh.read(4))[0]
    self.header_text = self.fh.read(l_text).rstrip('\0')
    self.n_ref = struct.unpack('<i',self.fh.read(4))[0]

def _parse_bam_data_block(bin_in,ref_names):
  v = {}
  data = StringIO(bin_in)
  rname_num = struct.unpack('<i',data.read(4))[0]
  v['rname'] = ref_names[rname_num] #refID to check in ref names
  v['pos'] = struct.unpack('<i',data.read(4))[0] + 1 #POS
  bin_mq_nl = struct.unpack('<I',data.read(4))[0]
  bin =  bin_mq_nl >> 16 
  v['mapq'] = (bin_mq_nl & 0xFF00) >> 8 #mapq
  l_read_name = bin_mq_nl & 0xFF #length of qname
  flag_nc = struct.unpack('<I',data.read(4))[0] #flag and n_cigar_op
  v['flag'] = flag_nc >> 16
  n_cigar_op = flag_nc & 0xFFFF
  l_seq = struct.unpack('<i',data.read(4))[0]
  rnext_num = struct.unpack('<i',data.read(4))[0]
  if rnext_num == -1:
    v['rnext'] = '*'
  else:
    v['rnext'] = ref_names[rnext_num] #next_refID in ref_names
  v['pnext'] = struct.unpack('<i',data.read(4))[0]+1 #pnext
  tlen = struct.unpack('<i',data.read(4))[0]
  v['tlen'] = tlen
  v['qname'] = data.read(l_read_name).rstrip('\0') #read_name or qname
  #print 'n_cigar_op '+str(n_cigar_op)
  v['cigar_bytes'] = data.read(n_cigar_op*4)
  #print 'cigar bytes '+str(len(v['cigar_bytes']))
  v['seq_bytes'] = data.read((l_seq+1)/2)
  v['qual_bytes'] = data.read(l_seq)
  v['extra_bytes'] = data.read()
  #last second tweak
  if v['rnext'] == v['rname']: v['rnext'] = '='
  return v

def _bin_to_qual(qual_bytes):
  if len(qual_bytes) == 0: return '*'
  if struct.unpack('<B',qual_bytes[0])[0] == 0xFF: return '*'
  #print qual_bytes
  #try:
  qual = ''.join([chr(struct.unpack('<B',x)[0]+33) for x in qual_bytes])
  #except:
  #  return '*'
  return qual

def _bin_to_seq(seq_bytes):
  if len(seq_bytes) == 0: return None
  global _bam_char
  #print len(seq_bytes)
  seq = ''.join([''.join([''.join([chr(z+97).translate(_bam_char) for z in  [y>>4,y&0xF]]) for y in struct.unpack('<B',x)]) for x in seq_bytes]).rstrip('=')
  return seq

def _bin_to_cigar(cigar_bytes):
  global _bam_ops
  if len(cigar_bytes) == 0: return [[],'*']
  cigar_packed = [struct.unpack('<I',x)[0] for x in \
             [cigar_bytes[i:i+4] for i in range(0,len(cigar_bytes),4)]]
  cigar_array = [[c >> 4, str(c &0xF).translate(_bam_ops)] for c in cigar_packed]
  cigar_seq = ''.join([''.join([str(x[0]),x[1]]) for x in cigar_array])
  return [cigar_array,cigar_seq]

#Pre all the reamining bytes of an entry
#Post an array of 
# 1. A dict keyed by Tag with {'type':,'value':} where value is a string unless type is i
# 2. A string of the remainder
def _bin_to_extra(extra_bytes):
  global _bam_value_type
  extra = StringIO(extra_bytes)
  tags = {}
  rem = ''
  while extra.tell() < len(extra_bytes):
    tag = extra.read(2)
    val_type = extra.read(1)
    if val_type == 'Z':
      rem += tag+':'
      rem += val_type+':'
      p = re.compile('([!-~])')
      m = p.match(extra.read(1))
      vre = ''
      while m:
        vre += m.group(1)
        c = extra.read(1)
        #print c
        m = p.match(c)
      rem += vre+"\t"
      tags[tag] = {'type':val_type,'value':vre}
    elif val_type == 'A':
      rem += tag+':'
      rem += val_type+':'
      vre = extra.read(1)
      rem += vre+"\t"      
      tags[tag] = {'type':val_type,'value':vre}      
    elif val_type in _bam_value_type:
      rem += tag+':'
      rem += 'i'+':'
      val = struct.unpack(_bam_value_type[val_type][1],extra.read(_bam_value_type[val_type][0]))[0]
      rem += str(val)+"\t"
      tags[tag] = {'type':val_type,'value':val}
    elif val_type == 'B':
      sys.sterr.write("WARNING array not implmented\n")
      continue
      rem += tag+':'
      rem += val_type+':'
      array_type = _bam_value_type[extra.read(1)]
      element_count = struct.unpack('<I',extra.read(4))[0]
      array_bytes = extra.read(element_count*_bam_value_type[array_type][0])
      for by in [array_bytes[i:i+_bam_value_type[array_type][0]] for i in range(0,len(array_bytes),_bam_value_type[array_type][0])]:
        aval = struct.unpack(_bam_value_type[array_type][1],by)
  return [tags,rem.rstrip("\t")]


class BGZF:
  # Methods adapted from biopython's bgzf.py
  def __init__(self,filename,blockStart=None,innerStart=None):
    self.path = filename
    self.fh = open(filename,'rb')
    if blockStart: self.fh.seek(blockStart)
    self._block_start = 0
    #self.pointer = 0
    #holds block_size and data
    self._buffer = self._load_block()
    self._buffer_pos = 0
    if innerStart: self._buffer_pos = innerStart
  def close(self):
    self.fh.close()
  def get_block_start(self):
    return self._block_start
  def get_inner_start(self):
    return self._buffer_pos
  def seek(self,blockStart,innerStart):
    self.fh.seek(blockStart)
    self._buffer_pos = 0
    self._buffer = self._load_block()
    self._buffer_pos = innerStart
  def read(self,size):
    done = 0 #number of bytes that have been read so far
    v = ''
    while True:
      if size-done < len(self._buffer['data']) - self._buffer_pos:
        v +=  self._buffer['data'][self._buffer_pos:self._buffer_pos+(size-done)]
        self._buffer_pos += (size-done)
        #self.pointer += size
        return v
      else: # we need more buffer
        vpart = self._buffer['data'][self._buffer_pos:]
        self._buffer = self._load_block()
        v += vpart
        self._buffer_pos = 0
        if len(self._buffer['data'])==0: return v
        done += len(vpart)

  def _load_block(self):
    #pointer_start = self.fh.tell()
    if not self.fh: return {'block_size':0,'data':''}
    self._block_start = self.fh.tell()
    magic = self.fh.read(4)
    if len(magic) < 4:
      #print 'end?'
      #print len(self.fh.read())
      return {'block_size':0,'data':''}
    gzip_mod_time, gzip_extra_flags, gzip_os,extra_len = struct.unpack("<LBBH",self.fh.read(8))
    pos = 0
    block_size = None
    #get block_size
    while pos < extra_len:
      subfield_id = self.fh.read(2)
      subfield_len = struct.unpack("<H",self.fh.read(2))[0]
      subfield_data = self.fh.read(subfield_len)
      pos += subfield_len+4
      if subfield_id == 'BC':
        block_size = struct.unpack("<H",subfield_data)[0]+1
    #block_size is determined
    deflate_size = block_size - 1 - extra_len - 19
    d = zlib.decompressobj(-15)
    data = d.decompress(self.fh.read(deflate_size))+d.flush()
    expected_crc = self.fh.read(4)
    expected_size = struct.unpack("<I",self.fh.read(4))[0]
    if expected_size != len(data):
      sys.stderr.write("ERROR unexpected size\n")
      sys.exit()
    crc = zlib.crc32(data)
    if crc < 0:  crc = struct.pack("<i",crc)
    else:  crc = struct.pack("<I",crc)
    if crc != expected_crc:
      sys.stderr.write("ERROR crc fail\n")
      sys.exit()
    return {'block_size':block_size, 'data':data}

class SamStream:
  #  minimum_intron_size greater than zero will only show sam entries with introns (junctions)
  #  minimum_overhang greater than zero will require some minimal edge support to consider an intron (junction)
  def __init__(self,fh=None,minimum_intron_size=0,minimum_overhang=0,reference=None):
    self.previous_line = None
    self.in_header = True
    self._reference = reference
    self.minimum_intron_size = minimum_intron_size
    self.minimum_overhang = minimum_overhang
    if minimum_intron_size <= 0:
      self.junction_only = False
    else:
      self.junction_only = True
      self.minimum_intron_size = minimum_intron_size
    self._header = None
    self.header_text = ''
    if fh:
      self.fh = fh
      self.assign_handle(fh)
      self.get_header()

  # return a string that is the header
  def get_header(self):
    if not self._header:
      self._header = SAMHeader(self.header_text)
      return self._header
    return self._header


  def set_junction_only(self,mybool=True):
    self.junction_only = mybool

  def assign_handle(self,fh):
    if self.in_header:
      while True:
        self.previous_line = fh.readline()
        if is_header(self.previous_line):
          self.header_text += self.previous_line
          #self.header.append(self.previous_line)
        else:
          self.in_header = False
          self.previous_line = self.previous_line
          break
      # make sure our first line is
      if self.junction_only:
        while True:
          if not self.previous_line: break
          if is_junction_line(self.previous_line,self.minimum_intron_size,self.minimum_overhang): break
          self.previous_line = self.fh.readline()

  def __iter__(self):
    return self

  def next(self):
    r = self.read_entry()
    if not r:
      raise StopIteration
    else:
      return r

  def read_entry(self):
    if not self.previous_line: return False
    out = self.previous_line
    self.previous_line = self.fh.readline()
    if self.junction_only:
      while True:
        if not self.previous_line: break
        if is_junction_line(self.previous_line,self.minimum_intron_size,self.minimum_overhang): break
        self.previous_line = self.fh.readline()
    if out:
      s = SAM(out,reference=self._reference)
      s.get_range()
      return s
    return None

def is_junction_line(line,minlen=68,minoverhang=0):
  prog = re.compile('([0-9]+)([NMX=])')
  f = line.rstrip().split("\t")
  v = prog.findall(f[5])
  #get the indecies of introns
  ns = [i for i in range(0,len(v)) if v[i][1]=='N' and int(v[i][0]) >= minlen]
  if len(ns) == 0: return False
  if minoverhang==0: return True
  good_enough = False
  for intron_index in ns:
    left = sum([int(x[0]) for x in v[0:intron_index] if x[1] != 'N'])
    right = sum([int(x[0]) for x in v[intron_index+1:] if x[1] != 'N'])
    worst = min(left,right)
    if worst >= minoverhang: good_enough = True
  if good_enough: return True
  return False


#pre: a flag from a sam file, in integer format
#     a bit to convert, given as a hex number ie 0x10
#post: returns true if the flag is set on
def is_header(line):
  if re.match('^@',line):
    f = line.rstrip().split("\t")
    if(len(f) > 9):
      return False
    return True
  return False


class SamtoolsBAMStream(SamStream):
  def __init__(self,path,minimum_intron_size=0,minimum_overhang=0,reference=None):
    self.previous_line = None
    self.in_header = True
    self._reference = reference
    self.minimum_intron_size = minimum_intron_size
    self.minimum_overhang = minimum_overhang
    if minimum_intron_size <= 0:
      self.junction_only = False
    else:
      self.junction_only = True
      self.minimum_intron_size = minimum_intron_size
    self.header_text = ''
    self._header = None
    self.path = path
    cmd = 'samtools view -h '+self.path
    self.fh_orig = Popen(cmd.split(),stdout=PIPE)
    self.fh = self.fh_orig.stdout
    self.assign_handle(self.fh)
    self.get_header()
  def close(self):
    self.fh_orig.communicate()
  #def write_index(self,opath,verbose=False):
  #  BamIndex.write_index(self.path,opath,verbose=verbose,samtools=True)
  #  #_write_index(self.path,opath,verbose=verbose,samtools=True)

def sort_header(header_text):
  #sort the chromosomes in a header text
  lines = header_text.rstrip().split("\n")
  rlens = {}
  for ln in lines:
    m = re.match('@SQ\tSN:(\S+)\tLN:(\S+)',ln)
    if m:
      rlens[m.group(1)] = m.group(2)
  output = ''
  done_lens = False
  for ln in lines:
    if re.match('@SQ\tSN:',ln):
      if not done_lens:
        done_lens = True
        for chr in sorted(rlens.keys()):
          output += "@SQ\tSN:"+chr+"\tLN:"+str(rlens[chr])+"\n"
    else:
      output += ln.rstrip("\n")+"\n"
  return output

def check_flag(flag,inbit):
  if flag & inbit: return True
  return False

