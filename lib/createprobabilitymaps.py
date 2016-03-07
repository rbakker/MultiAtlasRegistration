from fancypipe import *
import os, os.path as op
import nibabel
import numpy
import re

class InspectProbabilityMaps(FancyTask):
  title = 'Create various inspection variables and files to evaluate probability maps.'
  
  def main(self,pmfile):
    import scipy.misc
    nii = nibabel.load(pmfile)
    img = nii.get_data()
    cut = (
      numpy.argmax(img.sum(axis=1).sum(axis=1)),
      numpy.argmax(img.sum(axis=0).sum(axis=1)),
      numpy.argmax(img.sum(axis=0).sum(axis=0))
    )
    slices = (
      img[cut[0],:,:],
      img[:,cut[1],:],
      img[:,:,cut[2]]
    )
    sliceFiles = []
    basename = re.sub('\.nii(\.gz)?$','',op.basename(pmfile))
    dims = ('i','j','k')
    for d,dim in enumerate(dims):
      fname = self.tempfile('{}_{}{}.png'.format(basename,dims[d],cut[d]))
      scipy.misc.toimage(slices[d].squeeze()).save(fname)
      sliceFiles.append(fname)
    return FancyDict(
      sliceFiles = sliceFiles
    )
#endclass

  
class CreateProbabilityMaps(FancyTask):
  title = 'Construct a probabilistic atlas from multiple registrations of a brain region.'

  inputs = {
    'infiles' : {'type':assertList, 'default':None,
      'help':'List of input label nifti files (label registered to target space).'
    },
    'labelnames': {'type':assertList, 'default':[],
      'help':'List that maps label indices to label names. If omitted, label name for index X is "labelX"'
    },
    'outfilebase' : {
      'help':'Output file base name, it will be appended with _labelname for each label'
    },
    'groundtruth' : {'type':assertFile, 'default':None,
      'help':'Nifti file that contains manual segmentation for the target scan'
    }
  }
  
  def main(self,infiles,labelnames,outfilebase,groundtruth):
    # prepare
    N = len(infiles)
    assert N>1,\
      'Need more than {} label files to create probability map in node {}.'.format(N,self)
    nii = nibabel.load(infiles[0])
    hdr = nii.get_header()
    q = hdr.get_best_affine()
    img = nii.get_data()
    labels = numpy.unique(img)
    if not labelnames:
      labelnames = ['label{}'.format(b) for b in labels]
    probmap = { labelnames[b]:(img==b).astype('uint8') for b in labels }
    pixcount = { labelnames[b]:[numpy.count_nonzero(img==b)] for b in labels}
    
    for i in range(1,N):
      nii = nibabel.load(infiles[i])
      hdr = nii.get_header()
      assert (q == hdr.get_best_affine()).all(),\
        'Label files must all have the same affine transform in node {}.'.format(self)
      img = nii.get_data()
      for b in labels:
        probmap[labelnames[b]] += (img==b).astype('uint8')
        pixcount[labelnames[b]].append(numpy.count_nonzero(img==b))

    histograms = {}
    probCutoff = {}
    falsePosNeg = {}
    if groundtruth:
      pixcount_truth = {}
      falsePos_truth = {}
      falseNeg_truth = {}
      diceIndex_truth = {}
      nii = nibabel.load(groundtruth)
      img_truth = nii.get_data()
    for b,label in enumerate(labelnames):
      pm = probmap[label]
      # H contains voxel count for each discrete probability value
      H = numpy.array([numpy.count_nonzero(pm==p) for p in range(1,N+1)])
      # H times prob, total probability in each probability zone
      prob = numpy.arange(1,N+1,dtype=numpy.float)/N
      HxP = H*prob
      Vol = HxP.sum()
      # H times 1-prob, total probability in each probability zone
      Hx1_P = H*(1.-prob)

      # false negatives: pixels wrongly excluded if bins 0..i are discarded
      fn = HxP.cumsum()-HxP
      # false positives: pixels wrongly included if bins 0..i are discarded
      fp = Hx1_P[::-1].cumsum()[::-1]
      # interpolate such that fn == fp
      i = 0; 
      while fn[i+1]<fp[i+1]: i+=1
      di = (fn[i]-fp[i]) / (fn[i]-fp[i]+fp[i+1]-fn[i+1])
      Vol_i = H[::-1].cumsum()[::-1][i]
      Vol_i1 = H[::-1].cumsum()[::-1][i+1]
      pCutoff = prob[i]+di*(prob[i+1]-prob[i])
      fpCutoff = fp[i]+di*(fp[i+1]-fp[i])
      histograms[label] = H
      probCutoff[label] = pCutoff
      falsePosNeg[label] = fpCutoff
      if groundtruth:
        truth = (img_truth==b)
        region_i1 = (pm>i+1)
        region_i = (pm>i)
        fp_i1 = numpy.count_nonzero(numpy.logical_and(region_i1,truth==False))
        fn_i1 = numpy.count_nonzero(numpy.logical_and(truth,region_i1==False))
        fp_i = numpy.count_nonzero(numpy.logical_and(region_i,truth==False))
        fn_i = numpy.count_nonzero(numpy.logical_and(truth,region_i==False))
        Vol_i1 = numpy.count_nonzero(region_i1)
        Vol_i = numpy.count_nonzero(region_i)
        Vol = Vol_i+di*(Vol_i1-Vol_i)
        pixcount_truth[label] = numpy.count_nonzero(truth)
        falsePos_truth[label] = fp_i+di*(fp_i1-fp_i)
        falseNeg_truth[label] = fn_i+di*(fn_i1-fn_i)
        diceIndex_truth[label] = (2*(pixcount_truth[label]-falseNeg_truth[label]))/(pixcount_truth[label]+numpy.mean(pixcount[label]))
        # Equivalent formula:
        # diceIndex_truth[label] = (2*(numpy.mean(pixcount[label])-falsePos_truth[label]))/(pixcount_truth[label]+numpy.mean(pixcount[label]))

    volumes = odict(
      ('pixcount_mean', { label:numpy.mean(pixcount[label]) for label in labelnames }),
      ('pixcount_std', { label:numpy.std(pixcount[label],ddof=1) for label in labelnames }),
      ('pixcount_median', { label:numpy.median(pixcount[label]) for label in labelnames }),
      ('falsePosNeg_est', falsePosNeg),
      ('falsePosNeg_cutoff', probCutoff)
    )
    if groundtruth:
      volumes['pixcount_truth'] = pixcount_truth
      volumes['falsePos_truth'] = falsePos_truth
      volumes['falseNeg_truth'] = falseNeg_truth
      volumes['diceIndex_truth'] = diceIndex_truth
      
    outfiles = {}
    inspect = FancyDict()
    for label,pm in probmap.items():
      outfiles[label] = '{}_{}.nii.gz'.format(outfilebase,label)
      pm = pm.astype(numpy.float32)/N
      nibabel.nifti1.save(nibabel.nifti1.Nifti1Image(pm,q),outfiles[label])
      if label != 'background':
        inspect[label] = InspectProbabilityMaps().setInput(outfiles[label]).requestOutput('sliceFiles')
    return FancyDict(
      outfiles = outfiles,
      histograms = histograms,
      volumes = volumes,
      inspect = inspect
    )
#endclass


if __name__ == '__main__':
  CreateProbabilityMaps.fromCommandLine().run()  
