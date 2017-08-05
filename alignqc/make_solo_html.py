#!/usr/bin/env python
import sys, argparse, base64, re, os


def main(args):
  of = sys.stdout
  if args.output:
    of = open(args.output,'w')
  filecontents = open(args.input).read()
  imgs = []
  for m in re.finditer('(<\s*.*img.*>)',filecontents):
    imgs.append([m.start(),m.end()])
  prev = 0
  newcontents = ''
  for i in range(len(imgs)):
    newcontents += filecontents[prev:imgs[i][0]]
    newcontents += do_image_tag(filecontents[imgs[i][0]:imgs[i][1]],args)
    prev = imgs[i][1]
  newcontents += filecontents[prev:]

  filecontents = newcontents

  if args.all:
    a = []
    for m in re.finditer('(<a.*>)',filecontents):
      a.append([m.start(),m.end()])
    prev = 0
    newcontents = ''
    for i in range(len(a)):
      newcontents += filecontents[prev:a[i][0]]
      newcontents += do_a_tag(filecontents[a[i][0]:a[i][1]],args)
      prev = a[i][1]
    newcontents += filecontents[prev:]


  filecontents = newcontents
  styles = []
  for m in re.finditer('(<\s*.*type.*text/css.*>)',filecontents):
    styles.append([m.start(),m.end()])
  prev = 0    
  for i in range(len(styles)):
    of.write(filecontents[prev:styles[i][0]]+"\n")
    of.write(do_style_sheet(filecontents[styles[i][0]:styles[i][1]],args)+"\n")
    prev = styles[i][1]
  of.write(filecontents[prev:]+"\n")
  of.close()

def do_style_sheet(style_sheet,args):
  m=re.match('^(.*)(href\s*=\s*["\'][^"\']*["\'])(.*)$',style_sheet)
  if not m: 
    return style_sheet #cant replace for some reason
  start = m.group(1)
  finish = m.group(3)
  src_full = m.group(2)
  m = re.match('href\s*=\s*["\']([^"\']+)["\']',src_full)
  if not m:
    return style_sheet #cant replace for some reason
  srcpathpart = m.group(1)
  srcpath = os.path.dirname(args.input)+'/'+srcpathpart
  if not re.search('\.css',srcpath): return style_sheet
  if not os.path.isfile(srcpath): return style_sheet
  encoded = base64.b64encode(open(srcpath,'rb').read())
  disabler = "\n"
  if not args.all:
    disabler = "a {\n  pointer-events: none;\n}\n"
  return '<style type="text/css">'+disabler+"\n"+open(srcpath,'rb').read()+'</style>'

def do_a_tag(a_tag,args):
  m=re.match('^(.*)(href\s*=\s*["\'][^"\']*["\'])(.*)$',a_tag)
  #print m.group(2)
  if not m: 
    return a_tag #cant replace for some reason
  start = m.group(1)
  finish = m.group(3)
  src_full = m.group(2)
  m = re.match('href\s*=\s*["\']([^"\']+)["\']',src_full)
  if not m:
    return a_tag #cant replace for some reason
  srcpathpart = m.group(1)
  srcpath = os.path.dirname(args.input)+'/'+srcpathpart
  #if not re.search('\.png',srcpath): return img_tag
  if not os.path.isfile(srcpath): return a_tag
  #sys.stderr.write(srcpath+"\n")
  #sys.exit()
  encoded = base64.b64encode(open(srcpath,'rb').read())
  if re.search('\.pdf$',srcpath):  
    return start+' href="data:application/pdf;base64,'+encoded+'" '+finish
  if re.search('\.gz$',srcpath):  
    return start+' href="data:application/x-gzip;base64,'+encoded+'" '+finish
  return start+' href="data:text/plain;base64,'+encoded+'" '+finish

def do_image_tag(img_tag,args):
  m=re.match('^(.*)(src\s*=\s*["\'][^"\']*["\'])(.*)$',img_tag)
  if not m: 
    return img_tag #cant replace for some reason
  start = m.group(1)
  finish = m.group(3)
  src_full = m.group(2)
  m = re.match('src\s*=\s*["\']([^"\']+)["\']',src_full)
  if not m:
    return img_tag #cant replace for some reason
  srcpathpart = m.group(1)
  srcpath = os.path.dirname(args.input)+'/'+srcpathpart
  if not re.search('\.png',srcpath): return img_tag
  if not os.path.isfile(srcpath): return img_tag
  encoded = base64.b64encode(open(srcpath,'rb').read())
  return start+' src="data:image/png;base64,'+encoded+'" '+finish

def do_inputs():
  parser = argparse.ArgumentParser(description="Put css style sheets and PNG images into html file",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Input file")
  parser.add_argument('-o','--output',help="output file")
  parser.add_argument('-a','--all',action="store_true",help="replace other files as well")
  args = parser.parse_args()
  return args  

# for calling from another python script
def external(args):
  main(args)

def external_cmd(cmd):
  cache_argv = sys.argv
  sys.argv = cmd
  args = do_inputs()
  main(args)
  sys.argv = cache_argv

if __name__=="__main__":
  args = do_inputs()
  main(args)
