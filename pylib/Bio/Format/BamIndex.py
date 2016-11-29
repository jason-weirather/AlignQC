import gzip, sys, random, os
from Bio.Range import GenomicRange
from Bio.Format.Sam import BAMFile

# Index file is a gzipped TSV file with these fields:
# 1. qname
# 2. target range
# 3. bgzf file block start
# 4. bgzf inner block start
# 5. aligned base count
# 6. flag

# Usage:
# name_to_num is used to get all the names at random
# get_longest_target_alignment_coords_by_name is used to get the best random
# coord hash is import for random access
# There are some inactive methods because the datastructures
# they needed were not getting used and were memory intensive.
# subsequent updates could put them back or even better only use them
# when the methods requring them are called the first time
# This class is actually incredibly bulky for working with a big index
# > 1M reads.  I think some more specific cases may need to be written
class BAMIndex:
  def __init__(self,index_file):
    self.index_file = index_file
    self._name_to_num = {} #name index to line number
    self._lines = []
    self._coords = {} # get the one indexed line number from coordinates
    inf = gzip.open(self.index_file)
    linenum = 0
    for line in inf:
        f = line.rstrip("\n").split("\t")
        name = f[0]
        if name not in self._name_to_num:
          self._name_to_num[name] = []
        self._name_to_num[name].append(linenum)
        self._lines.append({'qname':f[0],'rng_str':f[1],'filestart':int(f[2]),'innerstart':int(f[3]),'basecount':int(f[4]),'flag':int(f[5])})
        linenum+=1
        if int(f[2]) not in self._coords: self._coords[int(f[2])] = {}
        self._coords[int(f[2])][int(f[3])] = linenum
    inf.close()
    return

  def destroy(self):
    self._name_to_num = None
    self._lines = None
    self._coords = None
    return

  # Pre: nothing
  # Post: True if each chromosome is listed together as a chunk and if the range starts go from smallest to largest
  #       otherwise false
  def check_ordered(self): 
    sys.stderr.write("error unimplemented check_ordered\n")
    sys.exit()
    seen_chrs = set()
    curr_chr = None
    prevstart = 0
    for l in self._lines:
      if not l['rng']: continue
      if l['rng'].chr != curr_chr:
        prevstart = 0
        if l['rng'].chr in seen_chrs:
          return False
        curr_chr = l['rng'].chr
        seen_chrs.add(curr_chr)
      if l['rng'].start < prevstart:  return False
      prevstart = l['rng'].start
    return True

  # Return how many entries have been indexed
  def get_length(self):
    return len(self._lines)

  def get_names(self):
    return self._name_to_num.keys()

  def get_coords_by_name(self,name):
    sys.stderr.write("error unimplemented get_coords_by_name\n")
    sys.exit()
    return [[self._lines[x]['filestart'],self._lines[x]['innerstart']] for x in self._queries[self._name_to_num[name]]]

  def get_longest_target_alignment_coords_by_name(self,name):
    longest = -1
    coord = None
    #for x in self._queries[self._name_to_num[name]]:
    for line in [self._lines[x] for x in self._name_to_num[name]]:
      if line['flag'] & 2304 == 0: 
        return [line['filestart'],line['innerstart']]
    return None
    sys.stderr.write("ERROR: no primary alignment set in index\n")
    sys.exit()

  # Tak the 1-indexed line number and return its index information
  def get_index_line(self,lnum):
    if lnum < 1: 
      sys.stderr.write("ERROR: line number should be greater than zero\n")
      sys.exit()
    elif lnum > len(self._lines):
      sys.stderr.write("ERROR: too far this line nuber is not in index\n")
      sys.exit()  
    return self._lines[lnum-1]

  # return the one-indexed line number given the coordinates
  def get_coord_line_number(self,coord):
    if coord[0] in self._coords:
      if coord[1] in self._coords[coord[0]]:
        return self._coords[coord[0]][coord[1]]
    return None

  def get_unaligned_lines(self):
    sys.stderr.write("error unimplemented get_unaligned_lines\n")
    sys.exit()
    return [self._lines[x-1] for x in self._unaligned]
    #return [x for x in self._lines if x['flag'] & 4]

  def get_unaligned_start_coord(self):
    sys.stderr.write("error unimplemented get_unaligned_start_coord\n")
    sys.exit()
    if len(self._unaligned)==0: return None
    return [self._lines[self._unaligned[0]-1]['filestart'],self._lines[self._unaligned[0]-1]['innerstart']]

  def get_range_start_coord(self,rng):
    sys.stderr.write("error unimplemented get_range_start_coord\n")
    sys.exit()
    if rng.chr not in self._chrs: return None
    for l in [self._lines[x-1] for x in self._chrs[rng.chr]]:
      ####
      y = l['rng']
      c = y.cmp(rng)
      if c > 0: return None
      if c == 0:
        x = y.get_payload()
        return [x[1],x[2]] # don't need the name
    return None
  
  # return the line number 1-indexed of the first occurance after range
  def get_range_start_line_number(self,rng):
    sys.stderr.write("error unimplemented get_range_start_line\n")
    sys.exit()
    for i in range(0,len(self._lines)):
      if rng.cmp(self._lines[i]['rng'])==0: return i+1
    return None

# The best index class will read an index file and only provide access
# to primary alignment coordinates
class BAMIndexRandomAccessPrimary:
  def __init__(self,index_file=None,alignment_file=None,verbose=False):
    self.verbose=verbose
    self.alignment_file = None
    if alignment_file: self.alignment_file = alignment_file
    elif os.path.exists(input_index):
      if os.path.exists(input_index[:-4]):
        self.alignment_file = input_index[:-4]
    self.index_file = None
    if index_file: self.index_file = index_file
    elif self.alignment_file:
      if os.path.exists(self.alignment_file+'.bgi'):
        self.index_file = self.alignment_file+'.bgi'
    if not index_file:
      sys.stderr.write("ERROR: Someway and somehow you need to define an index file.  Either through an alignment with one or directly or both\n")
      sys.exit()
    fh = gzip.open(index_file)
    self.bests = []
    z = 0
    tot = 0
    for line in fh:
      z += 1      
      f = line.rstrip().split("\t")
      if check_flag(int(f[5]),2304):
        continue # only process primary alignments
      tot += 1
      self.bests.append([int(f[2]),int(f[3])])
      if self.verbose and tot % 1000 == 0:
        sys.stderr.write(str(tot)+'/'+str(z)+' primary alignments read'+"\r")
    if self.verbose:
      sys.stderr.write("\n")
    return
  def destroy(self):
    self.bests = []
    return
  def get_random_coord(self):
    return random.choice(self.bests)
  #def get_alignment(self):
  #  if not self.alignment_file:
  #    sys.stderr.write("ERROR: alignment file needs to be defined on initialization for this method\n")
  #    sys.exit()
  #  v = random.choice(self.bests)
  #  bf = BAMFile(self.alignment_)


def check_flag(flag,inbit):
  if flag & inbit: return True
  return False
      
# Index file is a gzipped TSV file with these fields:
# 1. qname
# 2. target range
# 3. bgzf file block start
# 4. bgzf inner block start
# 5. aligned base count
# 6. flag
def write_index(path,index_file,verbose=False,samtools=False):
  if verbose:
    sys.stderr.write("scanning for primaries\n")
  reads = {}
  z = 0
  # force use of primary alignment flag if its not already used
  # require one and only one primary alignment for each read (or mate)
  fail_primary = False
  b2 = None
  if samtools:
    b2 = SamtoolsBAMStream(path)
  else:
    b2 = BAMFile(path)

  for e in b2:
    z+=1
    if verbose:
      if z %1000==0: sys.stderr.write(str(z)+"\r")
    name = e.value('qname')
    if name not in reads:
      reads[name] = {}
    type = 'u'
    if e.check_flag(64):
      type = 'l' #left mate
    elif e.check_flag(128):
      type = 'r' #right mate
    if not e.check_flag(2304):
      if type not in reads[name]: reads[name][type] = 0
      reads[name][type] += 1 # we have one
      if reads[name][type] > 1: 
        fail_primary = True
        break # too many primaries set to be useful
  # see if we have one primary set for each read
  for name in reads:
    for type in reads[name]:
      if reads[name][type] != 1:
        fail_primary = True
  if verbose:
    sys.stderr.write("\n")
  if fail_primary:
    sys.stderr.write("Failed to find a single primary for each read (or each mate).  Reading through bam to find best.\n")
    best = {}
    # must find the primary for each

    b2 = None
    if samtools:
      b2 = SamtoolsBAMStream(path)
    else:
      b2 = BAMFile(path)
    z = 0
    for e in b2:
      z += 1
      if verbose:
        if z %1000==0: sys.stderr.write(str(z)+"\r")
      name = e.value('qname')
      type = 'u'
      if e.check_flag(64):
        type = 'l' #left mate
      elif e.check_flag(128):
        type = 'r' #right mate
      if name not in best: best[name] = {}
      # get length
      l = 0
      if e.is_aligned():
        l = e.get_aligned_bases_count()
      if type not in best[name]: best[name][type] = {'line':z,'bpcnt':l}
      if l > best[name][type]['bpcnt']: 
        best[name][type]['bpcnt'] = l
        best[name][type]['line'] = z
    bestlinenumbers = set()
    for name in best:
      for type in best[name]:
        bestlinenumbers.add(best[name][type]['line'])
    if verbose:
      sys.stderr.write("\n")
  of = None
  try:
    of = gzip.open(index_file,'w')
  except IOError:
    sys.sterr.write("ERROR: could not find or create index\n")
    sys.exit()


  b2 = None
  if samtools:
    b2 = SamtoolsBAMStream(path)
  else:
    b2 = BAMFile(path)
  z = 0
  for e in b2:
    z+=1
    if verbose:
      if z%1000==0:
        sys.stderr.write(str(z)+" reads indexed\r")
    myflag = e.value('flag')
    if fail_primary: # see if this should be a primary
      if z not in bestlinenumbers:
        myflag = myflag | 2304
    rng = e.get_target_range()
    if rng: 
      l = e.get_aligned_bases_count()
      of.write(e.value('qname')+"\t"+rng.get_range_string()+"\t"+str(e.get_block_start())+"\t"+str(e.get_inner_start())+"\t"+str(l)+"\t"+str(myflag)+"\n")
    else: of.write(e.value('qname')+"\t"+''+"\t"+str(e.get_block_start())+"\t"+str(e.get_inner_start())+"\t"+'0'+"\t"+str(myflag)+"\n")
  sys.stderr.write("\n")
  of.close()
