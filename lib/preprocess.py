import sys
sys.path.append('../fancypipe')
from fancypipe import *
import os.path as op
import nibabel
import numpy

# preprocessing part of the hippocampus segmentation pipeline

def px2mm(niifile):
  nii = nibabel.load(niifile)
  q = nii.get_affine()
  return numpy.linalg.det(q)


class OrientToRAS(FancyTask):
  title = 'Rearrange data such that +i,+j,+k become (close to) Right,Anterior,Superior'
  inputs = {
    'inp': { 'type':assertFile, 'help':'Input image.' },
    'out': { 'help':'Output image.' }
  }
  
  def main(self,inp,out):
    nii = nibabel.load(inp)
    nii = nibabel.as_closest_canonical(nii)
    # remove trailing singleton dimension
    hdr = nii.get_header()
    shape = hdr.get_data_shape()
    if shape[-1] == 1:
      img = nii.get_data()
      img = img.squeeze(-1)
      nii = nibabel.nifti1.Nifti1Image(img,hdr.get_best_affine())
    nibabel.nifti1.save(nii, out)
    return FancyDict(
      out = out,
      dtype = nii.get_data_dtype()
    )
#endclass


class RegisterMask(FancyTask):
  title = "Register the MNI152 brain mask to a target image using NiftyReg."

  def main(self,inp,out,outoverlap,outaff,outf3d, mnit1,mnimask):
    title = 'Register MNI T1 to N4-corrected scan: affine'
    resultfile = self.tempfile('mni2target_aladin.nii.gz') if 'mni2target_aladin' in self.requests else ''
    mni2target_aladin = FancyExec().setProg('reg_aladin').setInput(
      **{
        '-flo' : mnit1,
        '-ref' : inp,
        '-ln'  : '3',
        '-lp'  : '2',
        '-res' : FancyOutputFile(resultfile),
        '-aff' : FancyOutputFile(outaff)
      }
    )
    
    title = 'Register MNI T1 to N4-corrected scan, free-form deformation'
    resultfile = self.tempfile('mni2target_f3d.nii.gz') if 'mni2target_f3d' in self.requests else ''
    mni2target_f3d = FancyExec().setProg('reg_f3d').setInput(
      **{
        '-flo' : mnit1,
        '-ref' : inp,
        '-ln'  : '3',
        '-lp'  : '2',
        '-res' : FancyOutputFile(resultfile),
        '-cpp' : FancyOutputFile(outf3d),
        '-aff' : mni2target_aladin.requestOutput('-aff')
      }
    )

    title = 'Apply same transformation to MNI mask'
    mnimask2target = FancyExec().setProg('reg_resample').setInput(
      '-LIN',
      **{
        '-flo' : mnimask,
        '-ref' : inp,
        '-res' : FancyOutputFile(out),
        '-cpp' : mni2target_f3d.requestOutput('-cpp')
      }
    )
    
    title = 'Apply same transformation to MNI bounding box'
    nii = nibabel.load(mnimask)
    hdr = nii.get_header()
    img = nii.get_data()
    img[:] = 1
    nii = nibabel.nifti1.Nifti1Image(img,hdr.get_best_affine())
    mnioverlap = self.tempfile('mnioverlap.nii.gz')
    nibabel.nifti1.save(nii,mnioverlap)
    mnioverlap2target = FancyExec().setProg('reg_resample').setInput(
      '-LIN',
      **{
        '-flo' : mnioverlap,
        '-ref' : inp,
        '-res' : FancyOutputFile(outoverlap),
        '-cpp' : mni2target_f3d.requestOutput('-cpp')
      }
    )
    
    return FancyDict(
      mnimask2target = mnimask2target.requestOutput('-res'),
      mnioverlap2target = mnioverlap2target.requestOutput('-res'),
      mni2target_aladin = mni2target_aladin.requestOutput('-res'),
      mni2target_f3d = mni2target_f3d.requestOutput('-res'),
      mniaffine = mni2target_aladin.requestOutput('-aff')
    )
#endclass

# create a mask using niftyreg
class RegAladin(FancyExec):
  title = 'Rough affine registration of the MNI152 template'
#endclass

class RunN4(FancyTask):
  title = 'N4 non-uniformity correction'
  description = 'Nick Tustison\'s N4 non-uniformity correction, applied to mask (i.e. brain tissue)'
  
  def main(self,inp,out,maxiter,convthr,dist,shrink,mask,mnit1,mnimask):
    # store input datatype
    nii = nibabel.load(inp)
    inp_dtype = nii.get_data_dtype()

    if not mask:
      mni2target = RegAladin().setProg('reg_aladin').setInput(
        **{
          '-flo' : mnit1,
          '-ref' : inp,
          '-ln'  : '3',
          '-lp'  : '2',
          '-res' : '',
          '-aff' : FancyOutputFile(self.tempfile('mni2target_aladin.txt'))
        }
      )
  
      title = 'Apply transformation to MNI152 mask'
      mnimask2target = FancyExec().setProg('reg_resample').setInput(
        '-LIN',
        **{
          '-flo' : str(mnimask),
          '-ref':str(inp),
          '-aff' : mni2target.requestOutput('-aff'),
          '-res' : FancyOutputFile(self.tempfile('mnimask2target.nii.gz'))
        }
      )
    #endif

    title = 'Apply the N4 bias field correction, using the supplied or generated brain mask'
    n4bfc = FancyExec().setProg('N4BiasFieldCorrection').setInput(
      **{
        '-i' : str(inp),
        '-o' : FancyOutputFile(self.tempfile('n4target.nii.gz')),
        '-c' : '[%s,%s]' % (maxiter,convthr),
        '-b' : '[%s]' % dist,
        '-s' : str(shrink),
        '-v' : '2',
        '-x' : str(mask) if mask else mnimask2target.requestOutput('-res')
      }
    )
    
    # retype data
    retype = RescaleRetype().setInput(
      inp   = n4bfc.requestOutput('-o'),
      out   = FancyOutputFile(out),
      mean  = None,
      std   = None,
      dtype = inp_dtype,
      mask  = None
    )

    return FancyDict(
      out = retype.requestOutput('out')
    )
#endclass

class RescaleRetype(FancyTask):
  title = 'Rescale image and/or retype data'
  inputs = {
    'inp': { 'type':assertFile, 'help':'Input image.' },
    'out': { 'type':assertOutputFile,'help':'Output image.' },
    'mean': { 'type':float, 'default':None, 'help':'New mean intensity.' },
    'std': { 'type':float, 'default':None, 'help':'New standard deviation.' },
    'dtype': { 'help':'Output datatype.', 'default': None, 'help':'New datatype (None to leave as is)' },
    'mask': { 'type':assertFile, 'help':'Use mask to compute mean and variance.' }
  }
  
  def main(self,inp,out,mean,std,dtype,mask):
    nii = nibabel.load(inp)
    hdr = nii.get_header()
    img = nii.get_data()
    if img.dtype.kind in 'ui':
      img = img.astype(numpy.float64)
      
    if mask:
      nii_mask = nibabel.load(mask)
      img_mask = nii.get_data().astype(numpy.bool)
      old_mean = img[img_mask].mean()
      old_std = img[img_mask].std()
    else:    
      old_mean = img.mean()
      old_std = img.std()
    
    if mean is not None and std is not None:
      # Scale image to new mean and std
      img -= old_mean
      img *= std/old_std
      img += mean
    
    # Convert to dtype and save
    if dtype is not None:
      if not isinstance(dtype,numpy.dtype): dtype = numpy.dtype(dtype)
      if dtype.kind in 'ui':
        mn = numpy.iinfo(dtype).min
        mx = numpy.iinfo(dtype).max
        img[img<mn] = mn
        img[img>mx] = mx
      img = img.astype(dtype)
    
    nii = nibabel.nifti1.Nifti1Image(img,hdr.get_best_affine())
    nibabel.nifti1.save(nii, out)
    
    return FancyDict(
      out = out
    )
#endclass
    

class EqualizeHistogram(FancyTask):
  title = 'Equalize histogram inside mask'
  inputs = dict(
    img = {'type':assertFile,
      'help':'Image file (nifti).'
    },
    mask = {'type':assertFile,
      'help':'Mask file (nifti).'
    },
    outfile = {'type':assertOutputFile,
      'help':'Output file (nifti).'
    },
    numBins = {'type':int,
      'default':256,
      'help':'Number of equalized histogram bins.'
    }
  )
  
  def main(self,imgfile,maskfile,outfile,numBins):
    nii = nibabel.load(imgfile)
    hdr = nii.get_header()
    img = nii.get_data()
    nii = nibabel.load(maskfile)
    mask = nii.get_data().astype(bool)
    data = img[mask]
    data = numpy.sort(data)
    if numBins<8: raise RuntimeError('Number of bins must be at least 8.')
    binWidth = (len(data)-1.0)/(numBins-2.0)
    bins = (numpy.array(range(0,numBins-1))*binWidth).astype(int)
    bins = data[bins]
    img = numpy.reshape((numpy.digitize(img.ravel(),bins,False)+numpy.digitize(img.ravel(),bins,True))/2,img.shape)
    dtype = numpy.uint8 if numBins<256 else numpy.uint16
    img = img.astype(dtype)
    nii =  nibabel.nifti1.Nifti1Image(img,hdr.get_best_affine())
    if not outfile: outfile = self.tempfile('{}_eq.nii.gz'.format(fu.basename(imgfile)))
    nibabel.save(nii,outfile)
    return FancyDict(outfile=outfile)
#endclass


class PreProcess(FancyTask):
  title = 'Apply non-uniformity correction and create brain mask'
  description = None
  inputs = odict(
    ('inp',{'type':assertFile, 'help':'Input image.'}),
    ('out',{'type':assertOutputFile,'help':'Output image.'}),
    ('outeq',{'type':assertOutputFile,'help':'Output image, histogram equalized.'}),
    ('outmask',{'type':assertOutputFile,'help':'Output brain mask.'}),
    ('outoverlap',{'type':assertOutputFile,'help':'Output mni overlap mask.'}),
    ('outaff',{'type':assertOutputFile,'help':'Output mni affine transform.'}),
    ('outf3d',{'type':assertOutputFile,'help':'Output mni f3d transform.'}),
    ('maxiter',{'type':int,'default': 200, 'help':'Maximum number of iterations.'}),
    ('convthr',{'type':float,'default':0.0001, 'help':'Convergence threshold, ratio of intensity changes between iterations.'}),
    ('dist',{'type':float,'default':50.0,'help':'Bspline control point spacing in mm.'}),
    ('shrink',{'type':int,'default':4,'help':'Shrink factor, amount of subsampling of the input image.'}),
    ('mask',{'type':assertFile,'default':None,'help':'Mask image to restrict N4 correction domain.'}),
    ('mnit1',{'type':assertFile,'help':'MN152 T1, used if no mask is specified.'}),
    ('mnimask',{'type':assertFile,'help':'MNI52 mask, used if no mask is specified.'}),
    ('zmuv',{'type':assertDict,'default':{},'help':'Mean/standard deviation/dtype/mask settings for RescaleRetype.'}),
    ('donuc',{'type':bool, 'default':True,'help':'Apply N4 non-uniformity correction - set to False if input is already N4 corrected.'})
  )
  
  def main(self,inp,out,outeq,outmask,outoverlap,outaff,outf3d, maxiter,convthr,dist,shrink,mask,mnit1,mnimask,zmuv,donuc):
    # orient to RAS
    orientras = OrientToRAS().setInput(
      inp     = inp,
      out     = self.tempfile('rasimage.nii.gz'),
    )

    # run N4
    if donuc is True:
      # apply N4
      runn4 = RunN4().setInput(
        inp     = orientras.requestOutput('out'),
        out     = FancyOutputFile(self.tempfile('n4image.nii.gz')),
        maxiter = maxiter,
        convthr = convthr,
        dist    = dist,
        shrink  = shrink,
        mask    = None,
        mnit1   = mnit1,
        mnimask = mnimask
      )
    else:
      # skip N4
      runn4 = orientras

    # register mask
    regmask = RegisterMask().setInput(
      inp        = runn4.requestOutput('out'),
      out        = FancyOutputFile(outmask),
      outoverlap = FancyOutputFile(outoverlap),
      outaff     = FancyOutputFile(outaff),
      outf3d     = FancyOutputFile(outf3d),
      mnit1      = mnit1,
      mnimask    = mnimask
    )
        
    # normalize and retype data
    rescale = RescaleRetype().setInput(
      **zmuv
    ).setInput(
      inp     = runn4.requestOutput('out'),
      out     = FancyOutputFile(out),
      dtype   = orientras.requestOutput('dtype'),
      mask    = regmask.requestOutput('mnioverlap2target')
    )

    # compute histogram equalized scan image
    eqhist = EqualizeHistogram().setInput(
      imgfile = runn4.requestOutput('out'),
      maskfile = regmask.requestOutput('mnimask2target'),
      outfile = outeq,
      numBins = 256
    )

    return FancyDict(
      out  = rescale.requestOutput('out'),
      outeq  = eqhist.requestOutput('outfile'),
      brainmask = regmask.requestOutput('mnimask2target'),
      mnioverlap = regmask.requestOutput('mnioverlap2target'),
      mni_aladin = regmask.requestOutput('mni2target_aladin'),
      mni_f3d = regmask.requestOutput('mni2target_f3d'),
      mniaffine = regmask.requestOutput('mniaffine')
    )
#endclass

if __name__ == '__main__':
  PreProcess.fromCommandLine().run()
