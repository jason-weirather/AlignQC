from math import sqrt 

def average(arr):
  if len(arr) == 0:
    sys.stderr.write("ERROR: no content in array to take average\n")
    sys.exit()
  if len(arr) == 1:  return arr[0]
  return float(sum(arr))/float(len(arr))

def median(arr):
  if len(arr) == 0:
    sys.stderr.write("ERROR: no content in array to take average\n")
    sys.exit()
  if len(arr) == 1: return arr[0]
  quot = len(arr)/2
  rem = len(arr)%2
  if rem != 0:
    return sorted(arr)[quot]
  return float(sum(sorted(arr)[quot-1:quot+1]))/float(2)

def standard_deviation(arr):
  return sqrt(variance(arr))

def variance(arr):
  avg = average(arr)
  return sum([(float(x)-avg)**2 for x in arr])/float(len(arr)-1)

def N50(arr):
  if len(arr) == 0:
    sys.stderr.write("ERROR: no content in array to take N50\n")
    sys.exit()
  tot = sum(arr)
  half = float(tot)/float(2)
  cummulative = 0
  for l in sorted(arr):
    cummulative += l
    if float(cummulative) > half: 
      return l
  sys.stderr.write("ERROR: problem finding M50\n")
  sys.exit()
