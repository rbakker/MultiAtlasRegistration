import re
import os
import os.path as op
from collections import OrderedDict

def basename(f):
  return re.search('(.*).nii(.gz)*',op.basename(f)).group(1)

## Read files from directory, using include and exclude regular expressions
def dir2files(d,include=None,exclude=None):
  matches = []
  for root,dirs,files in os.walk(d):
    matches.extend([ op.join(root,v) for v in files if (include and re.search(include,v)) and not (exclude and re.search(exclude,v)) ])
  return matches

## Convert a list of files to a list of basenames, using a regular expression pattern that MUST contain at least one (group).
def files2basenames(files,pattern):
  return [ re.search(pattern,op.basename(f)).group(1) for f in files ]

def sortedAtlasesFromDir(atlasdir,scaninclude,scanexclude,labelinclude,labelexclude):
  scanfiles = dir2files(atlasdir,include=scaninclude,exclude=scanexclude);
  labelfiles = dir2files(atlasdir,include=labelinclude,exclude=labelexclude);
  basenames_scan = files2basenames(scanfiles, pattern=scaninclude)
  basenames_label = files2basenames(labelfiles, pattern=labelinclude)
  diff = set(basenames_scan) ^ set(basenames_label)
  if diff:
    raise AssertionError('Some files in "{}" do not have a matching scan or label file:\n{}\nScan files: {}\nLabel files: {}'.format(atlasdir,diff,scanfiles,labelfiles))
  # sort by basenames
  _,scanfiles = zip(*sorted(zip(basenames_scan,scanfiles)))
  _,labelfiles = zip(*sorted(zip(basenames_label,labelfiles)))
  return (scanfiles,labelfiles)


class ElastixParamFile:
  def __init__(self,paramfile):
    with open(paramfile, "r") as fp:
      lines = fp.readlines()
    params = OrderedDict()
    for i,s in enumerate(lines):
      ans = re.match(r'^\s*\((\w+)\s(.*)\)\s+$',s)
      if ans: 
        params[ans.group(1)] = ans.group(2)
      else:
        params['#{}'.format(i)] = s
    self.params = params
      
  def __getitem__(self,key):
    return self.params[key]

  def __setitem__(self,key,value):
    self.params[key] = value
    
  def saveAs(self,newname, overwrite=False):
    if not overwrite and op.isfile(newname):
      raise FileExistsError('Elastix parameter file "{}" already exists.'.format(newname))
    with open(newname, "w") as fp:
      for k,v in self.params.iteritems():
        if k[0] == '#':
          fp.write(v)               
        else:
          fp.write('({} {})\n'.format(k,v))
#endclass

