from fancypipe import *
import os.path as op
import numpy
import nibabel
import collections

class CalculateVolume(FancyTask):
  title="Calculate volume of each labelled region."
  
  inputs = {
    'infile':  { 'type':assertFile, 'help': "The label file to calculate the volumes and mean of" },
    'labelnames': { 'type':assertList, 'default':[],
      'help':'List that maps label indices to label names. If omitted, label name for index X is "labelX"'
    }
  }
  
  def main(self,infile,labelnames):
    # initialize volumes based on the first label image
    nii = nibabel.load(infile)
    img = nii.get_data()
    labels = numpy.unique(img)
    if not len(labels) == len(labelnames):
      cc = collections.Counter(img.ravel())
      labels = [k for k,v in cc.items() if v > 1]
    if not labelnames:
      labelnames = ['label{}'.format(b) for b in labels]
    volume = { labelnames[b]:numpy.count_nonzero(img==b) for b in labels }

    return FancyDict(
      volume = volume
    )
#endclass

class CalculateVolumesMean(CalculateVolume):
  title="Calculate mean volume of each labelled region."
  
  inputs = {
    'infiles':  {'type':assertList, 'help': "The label files to calculate the volumes and mean of"},
    'labelnames': CalculateVolume.inputs['labelnames']
  }
  
  def main(self,infiles,labelnames):
    N = len(infiles)
    assert N>=1,\
      'Need more than {} label files to compute mean in node {}.'.format(N,self)
    nii = nibabel.load(infiles[0])
    hdr = nii.get_header()
    q = hdr.get_best_affine()
    V = CalculateVolume.main(self,infiles[0],labelnames)['volume']
    volumes = { k:[v] for k,v in V.items() }

    for i in range(1,N):
      nii = nibabel.load(infiles[i])
      hdr = nii.get_header()

      try:
        assert (q[1:3,1:3] == hdr.get_best_affine()[1:3,1:3]).all(),\
         'Label files must all have the same affine transform in node {}:\nQ1 {}\nQ2 {}'.format(self,q,hdr.get_best_affine())
      except AssertionError:
        print('WARNING: the affine matrices of the various label files are the same.')
      
      V = CalculateVolume.main(self,infiles[i],labelnames)['volume']
      for lbl in labelnames:
        volumes[lbl].append(V[lbl])
    
    meanVolume = { k:numpy.mean(v) for k,v in volumes.items() }
    medianVolume = { k:numpy.median(v) for k,v in volumes.items() }
    
    return FancyDict(
      meanVolume = meanVolume,
      medianVolume = medianVolume
    )
#endclass
  
if __name__ == '__main__':
  CalculateVolumesMean.fromCommandLine().run()
