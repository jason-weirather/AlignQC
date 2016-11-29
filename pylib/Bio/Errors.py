from Bio.Sequence import rc
import sys
### Error Analysis ####
# I am to describe errors at several levels
# 
# Errors in the query sequence
# 
# 1. Is a query base an error or not?
#    Probability - Sometimes it can be ambiguous which base is in error
# 
# 2. What is the basic type of error?
#    Mismatch
#    Insertion
#      Total insertion
#      Homopolymer insertion
#    Deletion
#      Total deletion
#        Before
#        After
#      Homopolymer deletion
#   sum of probabilities should add up to 1.)
#
# 3. What is the more specific error?
#    Mismatch type
#    insertion/deletion - Base, Length
#
valid_types = set(['match','mismatch','total_insertion','total_deletion','homopolymer_insertion','homopolymer_deletion'])

class ErrorProfileFactory:
  def __init__(self):
    self._alignment_errors = []
    self._target_context_errors = None
    self._query_context_errors = None
    self._general_errors = GeneralErrorStats()
    return

  def close(self):
    self._target_context_errors = None
    self._query_context_errors = None
    self._general_errors = None
    for ae in self._alignment_errors:
      ae.close()
    self._alignment_errors = None
      

  def add_alignment_errors(self,ae):
    self._target_context_errors = None
    self._query_context_errors = None
    self._alignment_errors.append(ae)
    self._general_errors.add_alignment_errors(ae)

  def add_alignment(self,align):
    self._target_context_errors = None
    self._query_context_errors = None
    ae = AlignmentErrors(align)
    self._alignment_errors.append(ae)
    self._general_errors.add_alignment_errors(ae)

  def get_alignment_errors(self):
    return self._general_errors

  def get_target_context_error_report(self):
    report = {}
    report['header'] = ['before','after','reference','query','fraction']
    report['data'] = []
    r = self.get_target_context_errors()
    for b in sorted(r.keys()):
      for a in sorted(r[b].keys()):
        for t in sorted(r[b][a]):
          for q in sorted(r[b][a]):
            v = 0
            if r[b][a][t]['total'] > 0:
              v = float(r[b][a][t]['types'][q])/float(r[b][a][t]['total'])
            report['data'].append([b,a,t,q,v])
    return report

  def get_min_context_count(self,context_type):
    cnt = 10000000000
    bases = ['A','C','G','T']
    basesplus = ['A','C','G','T','-']
    r = None
    if context_type == 'target':
      r = self.get_target_context_errors()
    elif context_type == 'query':
      r = self.get_query_context_errors()
    else:
      sys.stderr.write("ERROR incorrect context type\n")
      sys.exit()
    for b1 in bases:
      for b2 in bases:
        for b3 in basesplus:
          if r[b1][b2][b3]['total'] < cnt: cnt = r[b1][b2][b3]['total']
    return cnt

  def write_context_error_report(self,file,context_type):
    if context_type == 'target':
      r = self.get_target_context_error_report()
    elif context_type == 'query':
      r = self.get_query_context_error_report()
    else:
      sys.stderr.write("ERROR invalid type must be target or query\n")
      sys.exit()
    of = open(file,'w')
    of.write("\t".join(r['header'])+"\n")
    for row in r['data']:
      of.write("\t".join([str(x) for x in row])+"\n")
    return

  def get_query_context_error_report(self):
    report = {}
    report['header'] = ['before','after','reference','query','fraction']
    report['data'] = []
    r = self.get_query_context_errors()
    for b in sorted(r.keys()):
      for a in sorted(r[b].keys()):
        for t in sorted(r[b][a]):
          for q in sorted(r[b][a]):
            v = 0
            if r[b][a][q]['total'] > 0:
              v = float(r[b][a][q]['types'][t])/float(r[b][a][q]['total'])
            report['data'].append([b,a,t,q,v])
    return report

  def get_target_context_errors(self):
    if not self._target_context_errors:
      self.combine_context_errors()
    return self._target_context_errors

  def get_query_context_errors(self):
    if not self._query_context_errors:
      self.combine_context_errors()
    return self._query_context_errors

  def combine_context_errors(self):
    r = {}
    if self._target_context_errors: r = self._target_context_errors
    for k in [x.get_context_target_errors() for x in self._alignment_errors]:
      for b in k:
        if b not in r: r[b] = {}
        for c in k[b]:
          if c not in r[b]: r[b][c] = {}
          for a in k[b][c]:
            if a not in r[b][c]: 
              r[b][c][a] = {}
              r[b][c][a]['total'] = 0
              r[b][c][a]['types'] = {}
            r[b][c][a]['total'] += k[b][c][a]['total']
            for type in k[b][c][a]['types']:
              if type not in r[b][c][a]['types']: r[b][c][a]['types'][type] = 0
              r[b][c][a]['types'][type] += k[b][c][a]['types'][type]
    self._target_context_errors = r
    r = {}
    if self._query_context_errors: r = self._query_context_errors
    for k in [x.get_context_query_errors() for x in self._alignment_errors]:
      for b in k:
        if b not in r: r[b] = {}
        for c in k[b]:
          if c not in r[b]: r[b][c] = {}
          for a in k[b][c]:
            if a not in r[b][c]: 
              r[b][c][a] = {}
              r[b][c][a]['total'] = 0
              r[b][c][a]['types'] = {}
            r[b][c][a]['total'] += k[b][c][a]['total']
            for type in k[b][c][a]['types']:
              if type not in r[b][c][a]['types']: r[b][c][a]['types'][type] = 0
              r[b][c][a]['types'][type] += k[b][c][a]['types'][type]
    self._query_context_errors = r

  def __str__(self):
    return self.get_string()

  def get_string(self):
    ostr = ''
    ostr += str(len(self._alignment_errors))+" Alignments\n"
    ostr += 'Target: '+"\n"
    totbases = sum([len(x.get_target_sequence()) for x in self._alignment_errors])
    ostr += '  '+str(totbases)+" Target Bases\n"
    adjerror = sum([sum([y.get_error_probability() for y in x.get_target_errors()]) for x in self._alignment_errors])
    ostr += '  '+str(adjerror)+" Approximate error count\n"
    ostr += '  '+str(float(adjerror)/float(totbases))+" Error rate\n"
    ostr += 'Query: '+"\n"
    totbases = sum([len(x.get_query_sequence()) for x in self._alignment_errors])
    ostr += '  '+str(totbases)+" Query Bases\n"
    adjerror = sum([sum([y.get_error_probability() for y in x.get_query_errors()]) for x in self._alignment_errors])
    ostr += '  '+str(adjerror)+" Approximate error count\n"
    ostr += '  '+str(float(adjerror)/float(totbases))+" Error rate\n"
    return ostr

class BaseError():
  def __init__(self,type):
    self._type = type
    if type != 'query' and type != 'target':
      sys.stderr.write("ERROR specify type as query or target\n")
      sys.exit()
    self._unobservable = BaseError.UnobservableError(self._type)
    self._observable = BaseError.ObservableError(self._type)
    return
  def get_homopolymer(self):
    return self._observable.get_homopolymer()  
  # Is there any possible way to attribute this call to an error?
  def is_any_error(self):
    if self.get_observable_error_probability() > 0:
      return True
    if self.get_unobservable_error_probability() > 0:
      return True
    return False

  def get_observable(self):
    return self._observable
  def get_unobservable(self):
    return self._unobservable

  # This means for the base we are talking about how many errors between 0 and 1 do we attribute to it?
  # For the 'unobserved' errors, these can only count when one is adjacent to base
  def get_error_probability(self):
    a = self._observable.get_error_probability()
    b = self._unobservable.get_error_probability()
    return a+(1-a)*b

  def get_observable_error_probability(self):
    return self._observable.get_error_probability()
  def get_unobservable_error_probability(self):
    return self._unobservable.get_error_probability()

  def set_observable(self,tseq,qseq):
    tnt = None
    qnt = None
    if len(tseq) > 0: tnt = tseq[0]
    if len(qseq) > 0: qnt = qseq[0]
    self._observable.set(len(tseq),len(qseq),tnt,qnt)

  def set_unobserved_before(self,tlen,qlen,nt,p):
    self._unobservable.set_before(tlen,qlen,nt,p)
  def set_unobserved_after(self,tlen,qlen,nt,p):
    self._unobservable.set_after(tlen,qlen,nt,p)

  def get_adjusted_error_count(self):
    p1 =  self._observable.get_attributable_length()
    p1 += self._unobservable.get_attributable_length()
    return p1

  def get_base(self):
    if self._type == 'query':
      return self._observable.get_query_base()
    return self._observable.get_target_base()

  #def get_type(self):
  #  otype = self._observable.get_type()
  #  if otype[0] != 'match': return otype
  #  before = self._unobservable.get_before_type()
  #  after = self._unobservable.get_after_type()
  #  if before: return before
  #  if after: return after
  #  return otype    

  def __str__(self):
    return self.get_string()

  def get_string(self):
    ostr = ''
    ostr += 'BaseError for ['+self._type+'] base: '+self.get_base()+"\n"
    if self._observable.get_error_probability() > 0:
      ostr += '  Homopolymer set:'+"\n"
      ostr += '    '+str(self.get_homopolymer())+"\n"
      ostr += '  Observable:'+"\n"
      ostr += '    type is: '+str(self.get_observable_type())+"\n"
      ostr += '    P(error): '+str(self._observable.get_error_probability())+"\n"
      ostr += '    Elength:  '+str(self._observable.get_attributable_length())+"\n"
    before = self._unobservable.get_before_type()
    after = self._unobservable.get_after_type()
    if before or after:
      ostr += '  Unobservable '
      if self._type == 'query': ostr += 'deletion:'+"\n"
      else: ostr += 'insertion:'+"\n"
      ostr += '    P(error): '+str(self._unobservable.get_error_probability())+"\n"
      ostr += '    Elength:  '+str(self._unobservable.get_attributable_length())+"\n"
      if before:
        ostr += '    before: '+"\n"
        ostr += '      P(error): '+str(self._unobservable.get_before_probability())+"\n"
        ostr += '      '+str(before)+"\n"
      if after:
        ostr += '    after:'+"\n"
        ostr += '      P(error): '+str(self._unobservable.get_after_probability())+"\n"
        ostr += '      '+str(after)+"\n"
    return ostr


  # Unobservable error is a deletion for a query base
  #                     an insertion for a target base
  # A non base error has a probability of occuring before a base
  # and a probability of occuring after
  class UnobservableError:
   def __init__(self,type):
     # Type is the perspective
     self._type = type # Type is 'query' or 'target'
     if type != 'query' and type != 'target':
       sys.stderr.write("ERROR specify type as query or target\n")
       sys.exit()
     self._before_prob = 0  # probability specific base should have some missing call before it
     self._after_prob = 0   # probability specific base should have some missing call after it
     self._before = {'qlen':0,'tlen':0,'nt':None}
     self._after = {'qlen':0,'tlen':0,'nt':None}
   def get_error_probability(self):
     # P(before or after)
     return self._before_prob+(1-self._before_prob)*self._after_prob
   def set_after(self,tlen,qlen,nt,p):
     self._after = {'qlen':qlen,'tlen':tlen,'nt':nt}
     self._after_prob = float(p)
   def set_before(self,tlen,qlen,nt,p):
     self._before = {'qlen':qlen,'tlen':tlen,'nt':nt}
     self._before_prob = float(p)
   def get_before_type(self):
     if self._before_prob > 0 and self._type == 'query':
       return ['deletion','total_deletion',[\
                                            [self._before['nt'],self._before['tlen']],\
                                            [self._before['nt'],0]\
                                           ]\
              ]
     if self._before_prob > 0 and self._type == 'target':
       return ['insertion','total_insertion',[\
                                              [self._before['nt'],0],\
                                              [self._before['nt'],self._before['qlen']]\
                                             ]\
              ]
     return None
   def get_after_type(self):
     if self._after_prob > 0 and self._type == 'query':
       return ['deletion','total_deletion',[\
                                            [self._after['nt'],self._after['tlen']],\
                                            [self._after['nt'],0]\
                                           ]\
              ]
     if self._after_prob > 0 and self._type == 'target':
       return ['insertion','total_insertion',[\
                                              [self._after['nt'],0],\
                                              [self._after['nt'],self._after['qlen']]\
                                             ]\
              ]
     return None
   def get_after_probability(self):
     return self._after_prob
   def get_before_probability(self):
     return self._before_prob
   def get_attributable_length(self):
     bef = self._before_prob*abs(self._before['qlen']-self._before['tlen'])
     af = self._after_prob*abs(self._after['qlen']-self._after['tlen'])
     return bef+af

  # Class to describe a homopolymer error or an observable
  # insertion or deletion
  class ObservableError:
    def __init__(self,type):
      self._type = type # Type is 'query' or 'target'
      if type != 'query' and type != 'target':
        sys.stderr.write("ERROR specify type as query or target\n")
        sys.exit()
      self._prob = 0 # proportion of error we will attribute to this base
      self._details = {'qlen':0,'tlen':0,'qnt':None,'tnt':None}
    def set(self,tlen,qlen,tnt,qnt):
      self._details = {'qlen':qlen,'tlen':tlen,'qnt':qnt,'tnt':tnt}
      # can figure out probability now
      if qlen == tlen and qnt == tnt:
        self._prob = float(0)
      elif qnt != tnt:
        self._prob = float(1)
      else:
        delta = self.get_changed_length()
        if self._type == 'query':
          if delta > qlen:
            self._prob = float(1) # we could ascribe one or more insertions to each base so call them all an error
          else:
            self._prob = float(delta)/float(qlen)
        elif self._type == 'target':
          if delta > tlen:
            self._prob = float(1)
          else:
            self._prob = float(delta)/float(tlen)
        else:
          sys.stderr.write("unknown perspective type\n")
          sys.exit()
    def get_homopolymer(self):
      tnt = ''
      qnt = ''
      if self._details['tnt']: tnt = self._details['tnt']
      if self._details['qnt']: qnt = self._details['qnt']
      return {'tseq':self._details['tlen']*tnt,'qseq':self._details['qlen']*qnt}
    # Return the basic type of observable error
    def get_type(self):
      if self._details['tlen'] == self._details['qlen'] and\
         self._details['tnt'] == self._details['qnt']:
        return ['match','match',[[self._details['tnt'],1],[self._details['qnt'],1]]]
      if self._details['tlen'] == self._details['qlen'] and\
         self._details['tnt'] != self._details['qnt']:
        return ['mismatch','mismatch',[[self._details['tnt'],1],[self._details['qnt'],1]]]
      if self._details['tlen'] > self._details['qlen']:
        if self._details['qlen'] == 0:
          return ['deletion','total_deletion',[\
                                                [self._details['tnt'],self._details['tlen']],\
                                                [self._details['tnt'],0]\
                                              ]\
                 ]
        return ['deletion','homopolymer_deletion',[\
                                                   [self._details['tnt'],self._details['tlen']],\
                                                   [self._details['qnt'],self._details['qlen']]\
                                                  ]\
               ]
      if self._details['qlen'] > self._details['tlen']:
        if self._details['tlen'] == 0:
          return ['insertion','total_insertion',[\
                                                 [self._details['qnt'],0],\
                                                 [self._details['qnt'],self._details['qlen']]\
                                                ]\
                 ]
        return ['insertion','homopolymer_insertion',[\
                                                     [self._details['tnt'],self._details['tlen']],\
                                                     [self._details['qnt'],self._details['qlen']]\
                                                    ]\
               ]
      return 'UNKNOWN'

    def get_query_base(self):
      return self._details['qnt']        
    def get_target_base(self):
      return self._details['tnt']        
    def get_error_probability(self):
      return self._prob
    # For calculating total error counts
    def get_attributable_length(self):
      delta = self.get_changed_length()
      # calculate extra
      extra = 0
      if self._type == 'query' and self._details['qlen'] < delta:
        remainder = delta - self._details['qlen']
        extra = float(remainder)/float(self._details['qlen'])
      elif self._type == 'target' and self._details['tlen'] < delta:
        remainder = delta - self._details['tlen']
        extra = float(remainder)/float(self._details['tlen'])        
      return self._prob+extra

    def get_changed_length(self): #if we evenly distribute the length of the change in size, how much goes with this
      return abs(self._details['qlen']-self._details['tlen'])

class AlignmentErrors:
  # Pre: Take an alignment between a target and query
  #      Uses get_strand from alignment to orient the query
  #      All results are on the positive strand of the query
  #      (meaning may be the reverse complement of target if negative)
  def __init__(self,alignment,min_intron_size=68):
    #self._alns = []
    self._min_intron_size=min_intron_size
    self._aligned_query = None
    self._hpas = []
    self._has_quality = False # can be changed when add_alignment uses one that has quality
    self._alignment = alignment
    self._quality_distro = None # gets set by analyze_quality
    self._deletion_type = None
    self._query_errors = None
    self._target_errors = None
    self._context_query_errors = None
    self._context_target_errors = None
    astrings = self._alignment.get_alignment_strings(min_intron_size=self._min_intron_size)
    if self._alignment.get_query_quality(): self._has_quality = True
    if len(astrings) == 0: return None
    alns = []
    for i in range(len(astrings[0])):
      if self._alignment.get_strand() == '+':
        alns.append({'query':astrings[0][i],'target':astrings[1][i],'query_quality':astrings[2][i]})
      else:
        alns.insert(0,{'query':rc(astrings[0][i]),'target':rc(astrings[1][i]),'query_quality':astrings[2][i][::-1]})
    #if self._alignment.get_strand() == '-':
    #  alns = alns[::-1]
    #get homopolymer alignments
    self._hpas = self._misalign_split(alns) # split alignment into homopolymer groups
    self._query_hpas = []
    self._target_hpas = []
    qi = 0
    for i in range(len(self._hpas)):
      prev = None
      if i > 0: prev = self._hpas[i-1]
      foll = None
      if i + 1 < len(self._hpas): foll = self._hpas[i+1]
      qlen = len(self._hpas[i].get_query())
      for j in range(0,qlen):
        self._query_hpas.append({'hpa':self._hpas[i],'pos':j,'prev-hpa':prev,'next-hpa':foll})
      qi+=qlen
    ti = 0
    for i in range(len(self._hpas)):
      prev = None
      if i > 0: prev = self._hpas[i-1]
      foll = None
      if i + 1 < len(self._hpas): foll = self._hpas[i+1]
      tlen = len(self._hpas[i].get_target())
      for j in range(0,tlen):
        self._target_hpas.append({'hpa':self._hpas[i],'pos':j,'prev-hpa':prev,'next-hpa':foll})
      ti+=tlen
    self._target_errors = self.get_target_errors()
    self._query_errors = self.get_query_errors()  
    self._context_target_errors = self.get_context_target_errors()

  def close(self):
      self._min_intron_size= None
      self._aligned_query = None
      self._hpas = None
      self._has_quality = None # can be changed when add_alignment uses one that has quality
      self._alignment = None
      self._quality_distro = None # gets set by analyze_quality
      self._deletion_type = None
      self._query_errors = None
      self._target_errors = None
      self._context_query_errors = None
      self._context_target_errors = None
      self._hpas = None # split alignment into homopolymer groups
      self._query_hpas = None
      self._target_hpas = None
      self._target_errors = None
      self._query_errors = None
      self._context_target_errors = None
      return

  def get_HPAGroups(self):
    return self._hpas

  # way to accumulate totals of error types
  # General error report will be relative to to the total alignment length
  # error rate = mismatches + insertions + deletions / alignment length
  def get_general_errors(self):
    r = GeneralErrorStats()
    r.add_alignment_errors(self)

  def get_context_target_errors(self):
    if self._context_target_errors:  return self._context_target_errors
    if len(self._query_errors) < 3: return {}
    nts = ['A','C','G','T']
    poss = ['A','C','G','T','-']
    r = {}
    for i in nts:
      if i not in r:  r[i] = {}
      for j in nts:
        if j not in r[i]: r[i][j] = {}
        for k in poss:
          if k not in r[i][j]: 
            r[i][j][k] = {}
            r[i][j][k]['types'] = {}
            r[i][j][k]['total'] = 0
          for l in poss:
            if l not in r[i][j][k]['types']: r[i][j][k]['types'][l] = 0
    # now r is initialized
    for i in range(1,len(self._target_errors)-1):

      tobs = self._target_errors[i].get_observable()
      tunobs = self._target_errors[i].get_unobservable()
      otype = tobs.get_type()
      op = tobs.get_error_probability()
      before = tunobs.get_before_type()
      bp = tunobs.get_before_probability()
      after = tunobs.get_after_type()
      ap = tunobs.get_after_probability()
      if otype[2][0][0] == 'N': continue
      if otype[2][1][0] == 'N': continue
      if before:
        if before[2][0][0] == 'N': continue
        if before[2][1][0] == 'N': continue
      if after:
        if after[2][0][0] == 'N': continue
        if after[2][1][0] == 'N': continue

      tbefore = self._target_errors[i-1].get_base()
      t = self._target_errors[i].get_base()
      tafter = self._target_errors[i+1].get_base()

      if tbefore == 'N' or tafter == 'N' or t == 'N': continue
      r[tbefore][t]['-']['total'] += 0.5
      r[t][tafter]['-']['total'] += 0.5
      r[tbefore][tafter][t]['total'] += 1

      # We know we made an observation
      if otype[0] == 'mismatch':
        tb = otype[2][0][0]
        qb = otype[2][1][0]
        r[tbefore][tafter][t]['types'][qb] += op
      elif otype[0] == 'match':
        tb = otype[2][0][0]
        qb = otype[2][1][0]
        r[tbefore][tafter][t]['types'][qb] += float(1)
        #print op  
      elif otype[0] == 'deletion':
        tb = otype[2][0][0]
        qb = otype[2][1][0]
        r[tbefore][tafter][t]['types']['-'] += op
        r[tbefore][tafter][t]['types'][qb] += (1-op)
      # make sure our insertion can't be bigger than 1
      hp_insert_before = 0
      hp_insert_after = 0
      if otype[0] == 'insertion':
        tb = otype[2][0][0]
        qb = otype[2][1][0]
        r[tbefore][tb]['-']['types'][qb] += op/2
        r[tb][tafter]['-']['types'][qb] += op/2
        #homopolymer ... so we do have the correct base
        r[tbefore][tafter][t]['types'][qb] += 1

      # now take care of total insertions
      total_bp = 0
      total_ap = 0
      if before:
        qb = before[2][1][0]
        r[tbefore][t]['-']['types'][qb] += bp
      if after:
        qb = after[2][1][0]
        r[t][tafter]['-']['types'][qb] += ap
      #r[tbefore][t]['-']['types'][qb] += total_bp+(1-total_bp)*hp_insert_before
      #r[t][tafter]['-']['types'][qb] += total_ap+(1-total_ap)*hp_insert_after

      ##type = self._target_errors[i].get_type()
      #p = self._target_errors[i].get_error_probability()
      #if p > 0:
      #  if type[0] not in r[tbefore][t][tafter]['types']:
      #    r[tbefore][t][tafter]['types'][type[0]] = 0
      #  r[tbefore][t][tafter]['types'][type[0]] += p
      #r[tbefore][t][tafter]['total']+=1
    for b in r:
      for a in r:
        val = sum([r[b][a]['-']['types'][q] for q in nts])
        r[b][a]['-']['types']['-'] = r[b][a]['-']['total'] - val
    return r

  def get_context_query_errors(self):
    if self._context_query_errors:  return self._context_query_errors
    if len(self._query_errors) < 3: return {}
    nts = ['A','C','G','T']
    poss = ['A','C','G','T','-']
    r = {}
    for i in nts:
      if i not in r:  r[i] = {}
      for j in nts:
        if j not in r[i]: r[i][j] = {}
        for k in poss:
          if k not in r[i][j]: 
            r[i][j][k] = {}
            r[i][j][k]['types'] = {}
            r[i][j][k]['total'] = 0
          for l in poss:
            if l not in r[i][j][k]['types']: r[i][j][k]['types'][l] = 0
    # now r is initialized
    for i in range(1,len(self._query_errors)-1):

      tobs = self._query_errors[i].get_observable()
      tunobs = self._query_errors[i].get_unobservable()
      otype = tobs.get_type()
      op = tobs.get_error_probability()
      before = tunobs.get_before_type()
      bp = tunobs.get_before_probability()
      after = tunobs.get_after_type()
      ap = tunobs.get_after_probability()
      if otype[2][0][0] == 'N': continue
      if otype[2][1][0] == 'N': continue
      if before:
        if before[2][0][0] == 'N': continue
        if before[2][1][0] == 'N': continue
      if after:
        if after[2][0][0] == 'N': continue
        if after[2][1][0] == 'N': continue

      tbefore = self._query_errors[i-1].get_base()
      t = self._query_errors[i].get_base()
      tafter = self._query_errors[i+1].get_base()

      if tbefore == 'N' or tafter == 'N' or t == 'N': continue
      r[tbefore][t]['-']['total'] += 0.5
      r[t][tafter]['-']['total'] += 0.5
      r[tbefore][tafter][t]['total'] += 1

      # We know we made an observation
      if otype[0] == 'mismatch':
        tb = otype[2][0][0]
        qb = otype[2][1][0]
        r[tbefore][tafter][t]['types'][tb] += op
      elif otype[0] == 'match':
        tb = otype[2][0][0]
        qb = otype[2][1][0]
        r[tbefore][tafter][t]['types'][tb] += float(1)
        #print op  
      elif otype[0] == 'insertion':
        tb = otype[2][0][0]
        qb = otype[2][1][0]
        r[tbefore][tafter][t]['types']['-'] += op
        r[tbefore][tafter][t]['types'][tb] += (1-op)
      # make sure our deletion can't be bigger than 1
      hp_deletion_before = 0
      hp_deletion_after = 0
      if otype[0] == 'deletion':
        tb = otype[2][0][0]
        qb = otype[2][1][0]
        r[tbefore][tb]['-']['types'][tb] += op/2
        r[tb][tafter]['-']['types'][tb] += op/2
        #homopolymer ... so we do have the correct base
        r[tbefore][tafter][t]['types'][tb] += 1

      # now take care of total deletions
      if before:
        tb = before[2][0][0]
        r[tbefore][t]['-']['types'][tb] += bp
      if after:
        tb = after[2][0][0]
        r[t][tafter]['-']['types'][tb] += ap

      ##type = self._target_errors[i].get_type()
      #p = self._target_errors[i].get_error_probability()
      #if p > 0:
      #  if type[0] not in r[tbefore][t][tafter]['types']:
      #    r[tbefore][t][tafter]['types'][type[0]] = 0
      #  r[tbefore][t][tafter]['types'][type[0]] += p
      #r[tbefore][t][tafter]['total']+=1
    for b in r:
      for a in r:
        val = sum([r[b][a]['-']['types'][q] for q in nts])
        r[b][a]['-']['types']['-'] = r[b][a]['-']['total'] - val
    return r

  def get_query_errors(self):
    if self._query_errors: return self._query_errors
    v = []
    for i in range(len(self._query_hpas)):
      v.append(self.get_query_error(i))
    return v

  # Pre:  given an index in the aligned query
  # Post: return the error description for that base
  def get_query_error(self,i):
    x = self._query_hpas[i]
    h = x['hpa']
    pos = x['pos']
    prob = 0
    be = BaseError('query')
    be.set_observable(h.get_target(),h.get_query())
    #print pos
    if i != 0 and pos == 0: # check for a total deletion before
      prev = x['prev-hpa']
      if len(prev.get_query()) == 0: # total deletion
        be.set_unobserved_before(len(prev.get_target()),0,prev.get_target()[0],0.5)
    if i != len(self._query_hpas)-1 and pos == len(h.get_query())-1: # check for a total deletion before
      if x['next-hpa']:
        foll = x['next-hpa']
        if len(foll.get_query()) == 0: # total deletion
          be.set_unobserved_after(len(foll.get_target()),0,foll.get_target()[0],0.5)
    #else: print h
    return be

  def get_target_errors(self):
    if self._target_errors: return self._target_errors
    v = []
    for i in range(len(self._target_hpas)):
      v.append(self.get_target_error(i))
    return v

  # Pre:  given an index in the aligned query
  # Post: return the error description for that base
  def get_target_error(self,i):
    x = self._target_hpas[i]
    h = x['hpa']
    pos = x['pos']
    prob = 0
    be = BaseError('target')
    be.set_observable(h.get_target(),h.get_query())
    #print pos
    if i != 0 and pos == 0: # check for a total deletion before
      prev = x['prev-hpa']
      if len(prev.get_target()) == 0: # total insertion
        ilen = len(prev.get_query())
        be.set_unobserved_before(0,len(prev.get_query()),prev.get_query()[0],0.5)
    if i != len(self._target_hpas)-1 and pos == len(h.get_target())-1: # check for a total deletion before
      if x['next-hpa']:
        foll = x['next-hpa']
        if len(foll.get_target()) == 0: # total insertion
          be.set_unobserved_after(0,len(foll.get_query()),foll.get_query()[0],0.5)
    return be

  def get_query_sequence(self):
    return ''.join([x['hpa'].get_query()[0] for x in self._query_hpas])
  def get_target_sequence(self):
    return ''.join([x['hpa'].get_target()[0] for x in self._target_hpas])

  # Go through HPAGroups and store the distro of ordinal values of quality scores
  def analyze_quality(self):
    res = {}
    for h in self._hpas:
      if h.type() not in res: res[h.type()]={}
      for c in h.get_quality():
        if c not in res[h.type()]: res[h.type()][c] = 0
        res[h.type()][c]+=1
    self._quality_distro = res
  def get_quality_report_string(self):
    if not self._quality_distro:
      self.analyze_quality()
    ostr = ""
    for type in sorted(self._quality_distro.keys()):
      total = sum([ord(x)*self._quality_distro[type][x] for x in self._quality_distro[type]])
      cnt = sum([self._quality_distro[type][x] for x in self._quality_distro[type]])
      if cnt == 0: continue
      print 'type: '+type+' '+str(cnt)+' '+str(float(total)/float(cnt))
    return ostr

  def has_quality(self):
    return self._has_quality
  #Pre: alignment strings have been set so for each exon we have
  #     query, target and query_quality
  #     _has_quality will specify whether or not the quality is meaningful
  def _misalign_split(self,alns):
    total = []
    z = 0
    for x in alns:
      z += 1
      exon_num = z
      if self._alignment.get_strand() == '-':
        exon_num = (len(alns)-z)+1
      buffer = {'query':x['query'][0],'target':x['target'][0],'query_quality':x['query_quality'][0],'exon':exon_num}
      if buffer['query'] == '-': buffer['nt'] = buffer['target']
      elif buffer['target'] == '-': buffer['nt'] = buffer['query']
      elif buffer['query'] == buffer['target']: buffer['nt'] = buffer['query']
      elif buffer['query'] != buffer['target']: buffer['nt'] = '*'
      else:
	sys.stderr.write("WARNING unkonwn case\n")
      for i in range(1,len(x['query'])):
        qchar = x['query'][i]
        tchar = x['target'][i]
        qualchar = x['query_quality'][i]
        if qchar != tchar and (qchar != '-' and tchar != '-'):
          #classic mismatch
          #print 'mismatch'
          #print buffer
          total.append(buffer)
          buffer = {'query':qchar,'target':tchar,'query_quality':qualchar,'exon':exon_num}
          buffer['nt'] = '*'
        elif qchar == buffer['nt'] or tchar == buffer['nt']:
          # its a homopolymer match
          buffer['query'] += qchar

          buffer['target'] += tchar
          buffer['query_quality'] += qualchar
          #print 'homopoly'
        else:
          #print 'new thing'
          #print buffer
          total.append(buffer)
          buffer = {'query':qchar,'target':tchar,'query_quality':qualchar,'exon':exon_num}
          if qchar == '-': buffer['nt'] = tchar
          else:  buffer['nt'] = qchar
      total.append(buffer)
    result = [AlignmentErrors.HPAGroup(self,y) for y in total]
    return result

  #Homopolymer alignment group
  class HPAGroup:
    # takes a chunk of homopolymer alignment
    # as a dictionary with 'query' and 'target' sequences set
    # query should always be positive strand
    def __init__(self,parent,mydict):
      self._error_profile = parent
      self._data = mydict
      self._qseq = self._data['query'].replace('-','')
      self._tseq = self._data['target'].replace('-','')
      self._nt = self._data['nt'] # the nulceotide or * for mismatch
      self._qquality = self._data['query_quality'].replace('\0','')
      self._exon_number = self._data['exon']
      self._type = None
      if self._qseq == self._tseq:
        self._type = 'match'
      ### handle mismatches
      elif self._nt == '*':
        self._type = 'mismatch'
        self._code = self._tseq+'>'+self._qseq
      # Total deletion
      elif len(self._qseq) == 0:
        self._type = 'total_deletion'
      # Total insert
      elif len(self._tseq) == 0:
        self._type = 'total_insertion'
      elif len(self._qseq) < len(self._tseq):
        self._type = 'homopolymer_deletion'
      elif len(self._qseq) > len(self._tseq):
        self._type = 'homopolymer_insertion'
      else:
	sys.stderr.write("ERROR unsupported type\n")
        sys.exit()
    def get_nt(self): return self._data['nt']
    # always + strand
    def get_query(self):  return self._qseq
    # could be + or - strand
    def get_target(self):  return self._tseq
    def get_exon(self):  return self._exon_number
    def get_length(self):
      return {'query':len(self._qseq),'target':len(self._tseq)}
    def __str__(self):
      return self.get_string()
    def get_string(self):
      ostr = ''
      ostr += 'Target:  '+self._tseq+"\n"
      ostr += 'Query:   '+self._qseq+"\n"
      if self._error_profile.has_quality(): ostr += 'Quality: '+self._qquality+"\n"
      ostr += 'Type: '+str(self._type)+"\n"
      return ostr
    def has_quality(self):
      return self._error_profile.has_quality()
    def get_quality(self):
      if not self.has_quality(): return False
      return self._qquality
    def type(self):
      return self._type

# Keep track of general errors across the length of an alignment
class GeneralErrorStats:
  def __init__(self):
    self.alignment_count = 0 #number of alignments
    self.alignment_length = 0 #total bp
    self.mismatches = 0
    self.matches = 0
    self.deletions = {}
    self.deletions['total'] = 0
    self.deletions['specific'] = 0
    self.deletions['homopolymer'] = 0
    self.insertions = {}
    self.insertions['total'] = 0
    self.insertions['specific'] = 0
    self.insertions['homopolymer'] = 0
    possible = ['-','A','C','G','T']
    nts = ['A','C','G','T']
    #matrix holds the ref -> query changes, and the rates they occur
    self.matrix = {}
    for p1 in possible:
      self.matrix[p1] = {}
      for p2 in possible:
        self.matrix[p1][p2] = 0

  def __str__(self):
    return self.get_string()

  def get_string(self):
    ostr = ''
    errtotal = self.deletions['total']+self.insertions['total']+self.mismatches
    ostr += 'from '+str(self.alignment_length)+' bp of alignment'+"\n"
    ostr += '  '+str(float(errtotal)/float(self.alignment_length))+" error rate\n"
    ostr += '    '+str(float(self.mismatches)/float(self.alignment_length))+ " mismatches\n"
    ostr += '    '+str(float(self.deletions['total'])/float(self.alignment_length))+ " deletions\n"
    ostr += '      '+str(float(self.deletions['specific'])/float(self.alignment_length))+ " total deletions\n"
    ostr += '      '+str(float(self.deletions['homopolymer'])/float(self.alignment_length))+ " homopolymer deletions\n"
    ostr += '    '+str(float(self.insertions['total'])/float(self.alignment_length))+ " insertions\n"
    ostr += '      '+str(float(self.insertions['specific'])/float(self.alignment_length))+ " total insertions\n"
    ostr += '      '+str(float(self.insertions['homopolymer'])/float(self.alignment_length))+ " homopolymer insertions\n"
    ostr += '  More specific errors'+"\n"
    poss = ['-','A','C','G','T']
    ostr += '  -    A    C    G    T'+"\n"
    t = 0
    for p1 in poss:
      ostr += p1
      for p2 in poss:
        val = float(self.matrix[p1][p2])/float(self.alignment_length)
        ostr += " "+str(round(val,3))
        t += val
      ostr += "\n"
    ostr += "\n"
    return ostr

  def get_stats(self):
    ostr = ''
    errtotal = self.deletions['total']+self.insertions['total']+self.mismatches
    ostr += "ALIGNMENT_COUNT\t"+str(self.alignment_count)+"\n"
    ostr += "ALIGNMENT_BASES\t"+str(self.alignment_length)+"\n"
    ostr += "ANY_ERROR\t"+str(errtotal)+"\n"
    ostr += "MISMATCHES\t"+str(self.mismatches)+"\n"
    ostr += "ANY_DELETION\t"+str(self.deletions['total'])+"\n"
    ostr += "COMPLETE_DELETION\t"+str(self.deletions['specific'])+"\n"
    ostr += "HOMOPOLYMER_DELETION\t"+str(self.deletions['homopolymer'])+"\n"
    ostr += "ANY_INSERTION\t"+str(self.insertions['total'])+"\n"
    ostr += "COMPLETE_INSERTION\t"+str(self.insertions['specific'])+"\n"
    ostr += "HOMOPOLYMER_INSERTION\t"+str(self.insertions['homopolymer'])+"\n"
    return ostr

  def get_report(self):
    ostr = ''
    ostr += "target\tquery\tcnt\ttotal\n"
    poss = ['-','A','C','G','T']
    for target in poss:
      for query in poss:
        ostr += target+ "\t"+query+"\t"+str(self.matrix[target][query])+"\t"+str(self.alignment_length)+"\n"
    return ostr

  def add_alignment_errors(self,ae):
    self.alignment_count += 1
    for v in ae.get_HPAGroups():
      self._add_HPAGroup(v)

  def _add_HPAGroup(self,v):
    # Skip over N stuff
    if v.get_target():
      if v.get_target()[0] == 'N': return
    if v.get_query():
      if v.get_query()[0] == 'N': return

    l = max(v.get_length().values())
    self.alignment_length += l
    if v.type() == 'match':
      self.matches += l
      self.matrix[v.get_target()[0]][v.get_query()[0]]+=l
    elif v.type() == 'mismatch':
      self.mismatches += l
      self.matrix[v.get_target()[0]][v.get_query()[0]]+=l
    elif v.type() == 'total_deletion':
      self.deletions['total'] += v.get_length()['target']
      self.deletions['specific'] += v.get_length()['target']
      self.matrix[v.get_target()[0]]['-'] += v.get_length()['target']
    elif v.type() == 'homopolymer_deletion':
      self.deletions['total'] += v.get_length()['target']-v.get_length()['query']
      self.deletions['homopolymer'] += v.get_length()['target']-v.get_length()['query']
      self.matches += v.get_length()['query']
      self.matrix[v.get_target()[0]]['-'] += v.get_length()['target']-v.get_length()['query']
      self.matrix[v.get_target()[0]][v.get_query()[0]] += v.get_length()['query'] 
    elif v.type() == 'total_insertion':
      self.insertions['total'] += v.get_length()['query']
      self.insertions['specific'] += v.get_length()['query']
      self.matrix['-'][v.get_query()[0]] += v.get_length()['query']
    elif v.type() == 'homopolymer_insertion':
      self.insertions['total'] += v.get_length()['query']-v.get_length()['target']
      self.insertions['homopolymer'] += v.get_length()['query']-v.get_length()['target']
      self.matches += v.get_length()['target']
      self.matrix['-'][v.get_query()[0]] += v.get_length()['query']-v.get_length()['target']
      self.matrix[v.get_target()[0]][v.get_query()[0]] += v.get_length()['target'] 
    else:
      sys.stderr.write("ERROR unexpected error type: "+str(v.type())+"\n")
      sys.exit()
    return    
