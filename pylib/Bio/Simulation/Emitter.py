from Bio.Simulation.RandomSource import RandomSource

#Give it a transcriptome definition and a reference genome for it
#initialy give it uniform probability
class TranscriptomeEmitter:
  def __init__(self,transcriptome,seed=None,rand=None):
    if rand: self.random = rand
    elif seed: self.random = RandomSource(seed)
    else: self.random = RandomSource()

    self._transcriptome = transcriptome
    ######
    tcnt = len(self._transcriptome.get_transcripts())
    self._weights = [float(i+1)/float(tcnt) for i in range(0,tcnt)]
    ## _log stores what we are emitting ##
    self._log = []

  def emit_transcript(self):
    i = self.random.get_weighted_random_index(self._weights)
    return self._transcriptome.get_transcripts()[i]

  # input: an array of weights <<txname1> <weight1>> <<txname2> <weight2>>...
  def set_weights_by_dict(self,weights):
    self._weights = []
    txnames = [x.get_transcript_name() for x in self._transcriptome.get_transcripts()]
    for txname in txnames:
      if txname in weights:
        self._weights.append(float(weights[txname]))
      else:
        self._weights.append(float(0))
    return
