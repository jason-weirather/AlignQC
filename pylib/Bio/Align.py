import re, sys
from Bio.Sequence import rc
from Bio.Range import GenomicRange
from Bio.Structure import Transcript, Exon, Junction

from string import maketrans
# Basic class for common elements of alignments
# You don't have to have a query sequence and a reference sequence to do an alignment
# But 
class Alignment:
  def __init__(self):
    self._alignment_ranges = None #access through function because of BAM
    self._query_sequence = None
    self._query_quality = None
    self._query_length = 0
    self._target_length = 0
    self._reference = None
    self._set_alignment_ranges()
    return

  def get_aligned_bases_count(self):
    return sum([x[0].length() for x in self.get_alignment_ranges()])

  # These methods need to be overridden by an alignment type
  # target, query
  def _set_alignment_ranges(self):
    self._alignment_ranges = None
    sys.stderr.write("ERROR: needs overridden\n")
    sys.stderr.exit()
  # Post: Return sequence string, or None if not set
  def get_query_sequence(self):
    sys.stderr.write("ERROR: needs overridden\n")
    sys.stderr.exit()
    return self._query_sequence
  def set_query_sequence(self,seq):
    self._query_sequence = seq
  def get_target_length(self):
    sys.stderr.write("ERROR: needs overridden\n")
    sys.stderr.exit()
  def get_query_length(self):
    sys.stderr.write("ERROR: needs overridden\n")
    sys.stderr.exit()
  def get_strand(self):
    sys.stderr.write("ERROR: needs overridden\n")
    sys.stderr.exit()
    return self._query_direction

  # can be overridden
  # Post: Return quality string, or None if not set
  def get_query_quality(self):
    return self._query_quality

  def get_target_range(self):
    a = self.get_alignment_ranges()
    return GenomicRange(a[0][0].chr,a[0][0].start,a[-1][0].end)
  
  ## Range on the query string ... is the reverse complemented query if its on the negative strand
  #def get_query_range(self):
  #  a = self._alignment_ranges
  #  return GenomicRange(a[0][1].chr,a[0][1].start,a[-1][1].end,self.get_strand())

  # This is the actual query range for the positive strand
  def get_actual_query_range(self):
    a = self.get_alignment_ranges()
    #return GenomicRange(a[0][1].chr,a[0][1].start,a[-1][1].end,self.get_strand())
    if self.get_strand() == '+':
      return GenomicRange(a[0][1].chr,a[0][1].start,a[-1][1].end,self.get_strand())
    #must be - strand
    return GenomicRange(a[0][1].chr,self.get_query_length()-a[-1][1].end+1,self.get_query_length()-a[0][1].start+1,self.get_strand())

  def get_reference(self):
    return self._reference
  def set_reference(self,ref):
    self._reference = ref
  
  # Returns an array of alignment ranges
  def get_alignment_ranges(self):
    return self._alignment_ranges

  # Process the alignment to get information like
  # the alignment strings for each exon
  def get_alignment_strings(self,min_intron_size=68):
    qseq = self.get_query_sequence()
    if not qseq:
      sys.exit("ERROR: Query sequence must be accessable to get alignment strings\n")
      sys.exit()
    ref = self.get_reference()
    qual = self.get_query_quality()
    if not qual: 
      qual = 'I'*len(qseq) # for a placeholder quality
    if self.get_strand() == '-': 
      qseq = rc(qseq)
      qual = qual[::-1]
    tarr = []
    qarr = []
    yarr = []
    tdone = ''
    qdone = ''
    ydone = '' #query quality
    for i in range(len(self.get_alignment_ranges())):
      [t,q] = self.get_alignment_ranges()[i]
      textra = ''
      qextra = ''
      yextra = ''
      if i >= 1:
        dift = t.start-self.get_alignment_ranges()[i-1][0].end-1
        difq = q.start-self.get_alignment_ranges()[i-1][1].end-1
        if dift < min_intron_size:
          if dift > 0:
            textra = ref[t.chr][t.start-dift-1:t.start-1].upper()
            qextra = '-'*dift
            yextra = '\0'*dift
          elif difq > 0:
            textra = '-'*difq
            qextra = qseq[q.start-difq-1:q.start-1].upper()
            yextra = qual[q.start-difq-1:q.start-1]
        else:
          tarr.append(tdone)
          qarr.append(qdone)
          yarr.append(ydone)
          tdone = ''
          qdone = ''
          ydone = ''
      tdone += textra+ref[t.chr][t.start-1:t.end].upper()
      qdone += qextra+qseq[q.start-1:q.end].upper()
      ydone += yextra+qual[q.start-1:q.end]
    if len(tdone) > 0: 
      tarr.append(tdone)
      qarr.append(qdone)
      yarr.append(ydone)
    if self.get_query_quality() == '*': yarr = [x.replace('I',' ') for x in yarr]
    #query, target, query_quality
    return [qarr,tarr,yarr]

  def _analyze_alignment(self,min_intron_size=68):
    [qstrs,tstrs,ystrs] = self.get_alignment_strings(min_intron_size=min_intron_size)
    matches = sum([x[0].length() for x in self.get_alignment_ranges()]) 
    misMatches = 0
    for i in range(len(qstrs)):
      misMatches += sum([int(qstrs[i][j]!=tstrs[i][j] and qstrs[i][j]!='-' and tstrs[i][j]!='-' and tstrs[i][j]!='N') for j in range(len(qstrs[i]))])
    nCount = sum([len([y for y in list(x) if y == 'N'])  for x in tstrs])
    qNumInsert = sum([len(re.findall('[-]+',x)) for x in tstrs])
    qBaseInsert = sum([len(re.findall('[-]',x)) for x in tstrs])
    tNumInsert = sum([len(re.findall('[-]+',x)) for x in qstrs])
    tBaseInsert = sum([len(re.findall('[-]',x)) for x in qstrs])
    matches = matches - misMatches - nCount
    return {'matches':matches,\
            'misMatches':misMatches,\
            'nCount':nCount,\
            'qNumInsert':qNumInsert,\
            'qBaseInsert':qBaseInsert,\
            'tNumInsert':tNumInsert,\
            'tBaseInsert':tBaseInsert}

  # Pre: Have the alignment strings
  #      have get_query_sequence()
  #       and get_reference()
  def print_alignment(self,chunk_size=40,min_intron_size=68):
    has_qual = True
    if not self.get_query_quality(): has_qual = False
    trantab = maketrans('01',' *')
    [qstrs,tstrs,ystrs] = self.get_alignment_strings(min_intron_size=min_intron_size)
    print 'Alignment for Q: '+str(self.get_alignment_ranges()[0][1].chr)
    for i in range(len(qstrs)):
      print 'Exon '+str(i+1)
      #+' T: '+self._alignment_ranges[i][0].get_range_string()+' Q: '+str(self._alignment_ranges[i][1].start)+'-'+str(self._alignment_ranges[i][1].end)
      mm = ''.join([str(int(qstrs[i][j]!=tstrs[i][j] and qstrs[i][j]!='-' and tstrs[i][j]!='-' and tstrs[i][j]!='N')) for j in range(len(qstrs[i]))]).translate(trantab)
      t =  tstrs[i] #target
      q = qstrs[i]  #query
      s = ystrs[i] #quality
      for y in [[mm[x:x+chunk_size],t[x:x+chunk_size],q[x:x+chunk_size],s[x:x+chunk_size]] for x in range(0,len(mm),chunk_size)]:
        print '  '+y[0]
        print 'T '+y[1]
        print 'Q '+y[2]
        if has_qual: print 'Y '+y[3]        
        print ''

  # These methods may be overrriden by an alignment type to just return themself
  # clearly this should be over written by the PSL type to just give itself
  def get_PSL(self,min_intron_size=68):
    from Bio.Format.PSL import PSL
    if not self.get_alignment_ranges(): return None
    matches = sum([x[0].length() for x in self.get_alignment_ranges()]) # 1. Matches - Number of matching bases that aren't repeats
    misMatches = 0 # 2. Mismatches - Number of baess that don't match
    repMatches = 0 # 3. repMatches - Number of matching baess that are part of repeats
    nCount = 0 # 4. nCount - Number of 'N' bases
    qNumInsert = 0 # 5. qNumInsert - Number of inserts in query
    qBaseInsert = 0 # 6. qBaseInsert - Number of bases inserted into query
    tNumInsert = 0 # 7. Number of inserts in target
    tBaseInsert = 0 # 8. Number of bases inserted into target
    sub = self.get_query_sequence()
    ref = self.get_reference()
    if sub and ref:
      v = self._analyze_alignment(min_intron_size=min_intron_size)
      matches = v['matches']
      misMatches = v['misMatches'] # 2. Mismatches - Number of baess that don't match
      nCount = v['nCount'] # 4. nCount - Number of 'N' bases
      qNumInsert = v['qNumInsert'] # 5. qNumInsert - Number of inserts in query
      qBaseInsert = v['qBaseInsert'] # 6. qBaseInsert - Number of bases inserted into query
      tNumInsert = v['tNumInsert'] # 7. Number of inserts in target
      tBaseInsert = v['tBaseInsert'] # 8. Number of bases inserted into target
    strand = self.get_strand() # 9. strand 
    qName = self.get_alignment_ranges()[0][1].chr # 10. qName - Query sequence name
    qSize = self.get_query_length()
    qStart = self.get_alignment_ranges()[0][1].start-1
    qEnd = self.get_alignment_ranges()[-1][1].end
    tName = self.get_alignment_ranges()[0][0].chr
    tSize = self.get_target_length()
    tStart = self.get_alignment_ranges()[0][0].start-1
    tEnd = self.get_alignment_ranges()[-1][0].end
    blockCount = len(self.get_alignment_ranges())
    blockSizes = ','.join([str(x[0].length()) for x in self.get_alignment_ranges()])+','
    qStarts = ','.join([str(x[1].start-1) for x in self.get_alignment_ranges()])+','
    tStarts = ','.join([str(x[0].start-1) for x in self.get_alignment_ranges()])+','

    psl_string = str(matches)+"\t"+\
    str(misMatches)+"\t"+\
    str(repMatches)+"\t"+\
    str(nCount)+"\t"+\
    str(qNumInsert)+"\t"+\
    str(qBaseInsert)+"\t"+\
    str(tNumInsert)+"\t"+\
    str(tBaseInsert)+"\t"+\
    strand+"\t"+\
    qName+"\t"+\
    str(qSize)+"\t"+\
    str(qStart)+"\t"+\
    str(qEnd)+"\t"+\
    tName+"\t"+\
    str(tSize)+"\t"+\
    str(tStart)+"\t"+\
    str(tEnd)+"\t"+\
    str(blockCount)+"\t"+\
    blockSizes+"\t"+\
    qStarts+"\t"+\
    tStarts
    return PSL(psl_string,query_sequence=self.get_query_sequence(),reference=self.get_reference(),query_quality=self.get_query_quality())

  #clearly this should be overwritten by the SAM class to give itself
  def get_SAM(self,min_intron_size=68):
    from Bio.Format.Sam import SAM
    #ar is target then query
    qname = self.get_alignment_ranges()[0][1].chr
    flag = 0
    if self.get_strand() == '-': flag = 16
    rname = self.get_alignment_ranges()[0][0].chr
    pos = self.get_alignment_ranges()[0][0].start
    mapq = 255
    cigar = self.construct_cigar(min_intron_size)
    rnext = '*'
    pnext = 0
    tlen = self.get_target_range().length()
    seq = self.get_query_sequence()
    if not seq: seq = '*'
    qual = self.get_query_quality()
    if not qual: qual = '*'
    #seq = '*'
    #qual = '*'
    if self.get_strand() == '-':
      seq = rc(seq)
      qual = qual[::-1]
    ln = qname + "\t" + str(flag) + "\t" + rname + "\t" + \
         str(pos) + "\t" + str(mapq) + "\t" + cigar + "\t" + \
         rnext + "\t" + str(pnext) + "\t" + str(tlen) + "\t" + \
         seq + "\t" + qual
    return SAM(ln,reference=self._reference)

  def construct_cigar(self,min_intron_size=68):
    # goes target query
    ar = self.get_alignment_ranges()
    cig = ''
    if ar[0][1].start > 1: # soft clipped
      cig += str(ar[0][1].start-1)+'S'
    for i in range(len(ar)):
      exlen = ar[i][0].length()
      cig += str(exlen)+'M'
      if i < len(ar)-1:
        # we can look at distances
        dt = ar[i+1][0].start-ar[i][0].end-1
        dq = ar[i+1][1].start-ar[i][1].end-1
        if dq > 0: cig += str(dq)+'I'
        if dt >= min_intron_size:
          cig += str(dt)+'N'
        elif dt > 0: cig += str(dt)+'D'
        elif dq <= 0:
          sys.stderr.write("ERROR cant form alignment\n")
          sys.exit()

    if ar[-1][1].end < self.get_query_length(): # soft clipped
      cig += str(self.get_query_length()-ar[-1][1].end)+'S'
    return cig

  def get_target_transcript(self,min_intron=1):
    if min_intron < 1: 
      sys.stderr.write("ERROR minimum intron should be 1 base or longer\n")
      sys.exit()
    tx = Transcript()
    rngs = [self.get_alignment_ranges()[0][0].copy()]
    rngs[0].direction = None
    for i in range(len(self.get_alignment_ranges())-1):
      dist = self.get_alignment_ranges()[i+1][0].start - rngs[-1].end-1
      #print 'dist '+str(dist)
      if dist >= min_intron:
        rngs.append(self.get_alignment_ranges()[i+1][0].copy())
        rngs[-1].direction = None
      else:
        rngs[-1].end = self.get_alignment_ranges()[i+1][0].end
    tx.set_exons_and_junctions_from_ranges(rngs)
    tx.set_range()
    tx.set_strand(self.get_strand())
    tx.set_transcript_name(self.get_alignment_ranges()[0][1].chr)
    tx.set_gene_name(self.get_alignment_ranges()[0][1].chr)
    return tx
