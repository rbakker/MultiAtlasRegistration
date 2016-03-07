import sys
sys.path.append('../fancypipe')
sys.path.append('../lib')
from fancypipe import *
import util
import os.path as op
import re
import numpy, nibabel

class RegAladin2Elastix(FancyTask):
  def outputFile(self,key):
    if key == 'tpfile': return self.tempfile('tpfile_regaladin.txt')
  
  title = "Convert (niftireg) affine matrix to Elastix/ITK transformation param. file."
  def affine2transform(self,aff,ctr):
    m = 3 # fixed image dimension
    L = numpy.diag([-1,-1,1]) # ras2lps = lps2ras
    Ab_nr = aff
    A_nr = Ab_nr[0:m,0:m]
    b_nr = Ab_nr[0:m,m]
    A_ex = L.dot(A_nr).dot(L)
    c_ex = ctr
    t_ex = L.dot(b_nr)+A_ex.dot(c_ex)-c_ex
    return '{}'.format(numpy.append(A_ex,t_ex))[1:-1].replace('\n',' ')
  
  # tpbase is the Elastix transformation parameter file in which to
  #   insert the reg_aladin transform,
  # afffile must contain the reg_aladin affine transformation matrix.
  def main(self,tpbase,afffile):
    tp = util.ElastixParamFile(tpbase)
    ctr = numpy.fromstring(tp['CenterOfRotationPoint'],sep=' ')
    tp['TransformParameters'] = self.affine2transform(numpy.loadtxt(afffile),ctr)
    tp['Transform'] = '"AffineTransform"'
    tp['NumberOfParameters'] = '12'
    tpfile = self.outputFile('tpfile')
    tp.saveAs(tpfile, overwrite=True)    
    return FancyDict(
      tpfile = tpfile
    )
#endclass


class NiftyReg(FancyExec):
  title="Call <reg_aladin> with optional mask and initial transform."
  inputs = {
    'ln': { 'default':'3' },
    'lp': { 'default':'2' },
    'outfile': { 'type': assertOutputFile },
    'fixedmask': { 'default':None },
    'movingmask': { 'default':None },
    'inaff': { 'default':None }
  }

  def main(self,fixed,moving,ln,lp,outfile,fixedmask,movingmask,inaff):
    kwargs = {
      '-ref' : fixed,
      '-flo' : moving,
      '-aff' : FancyOutputFile(outfile),
      '-ln'  : ln,
      '-lp'  : lp,
      '-res' : ''
    }
    if fixedmask: kwargs['-rmask'] = fixedmask
    if movingmask: kwargs['-fmask'] = movingmask
    if inaff: kwargs['-inaff'] = inaff
    return FancyExec.main(
      self,
      'reg_aladin',
      **kwargs
    )
#endclass


class Elastix(FancyExec):
  title= "Call <elastix> with multiple parameter files."
  inputs = {
    'fixedmask': { 'default':None },
    'movingmask': { 'default':None },
    'initialtransform': { 'default':None },
    'numThreads': { 'default':None }
  }
  
  def outputFile(self,key):
    if key=='log':
      return self.tempfile('elastix.log')
    else:
      m = re.match('tp(\d+)',key)
      if m:
        return self.tempfile('TransformParameters.{}.txt'.format(m.group(1)))
    
  def main(self,fixed,moving,paramfiles, fixedmask,movingmask,initialtransform,numThreads):
    args = []
    try:
      rp0 = util.ElastixParamFile(paramfiles[0])
      t0 = rp0['InitialTransformParametersFileName']
    except:
      t0 = initialtransform
    for pf in paramfiles:
      args.extend(['-p',pf])
    kwargs = {
      '-f': fixed,
      '-m': moving,
      '-out': self.tempdir()
    }
    if t0: kwargs['-t0'] = t0
    if fixedmask: kwargs['-fMask'] = fixedmask
    if movingmask: kwargs['-mMask'] = movingmask
    if numThreads: kwargs['-threads'] = numThreads
    FancyExec.main(
      self,
      'elastix',
      *args,
      **kwargs
    )

    output = FancyDict(log = self.outputFile('log'))
    for p in range(len(paramfiles)):
      output['tp{}'.format(p)] = self.outputFile('tp{}'.format(p))
    return output
#endclass


class BackgroundFromMask(FancyTask):
  def main(self,imgfile,maskfile,outfile, bg=0):
    nii = nibabel.load(imgfile)
    hdr = nii.get_header()
    img = nii.get_data()
    nii = nibabel.load(maskfile)
    mask = nii.get_data().astype(bool)
    img[numpy.logical_not(mask)] = bg
    nii = nibabel.nifti1.Nifti1Image(img,hdr.get_best_affine())
    nibabel.nifti1.save(nii, outfile)
    return FancyDict(outfile=outfile)
#endclass
  

class MultiRegistration(FancyTask):
  title = "Register multiple atlases to a target scan."
  
  def main(self,fixedimg,movingimgs,paramfiles,prog,fixedmasks,movingmasks,finalmask=-1,bgfrommask=None):
    tpfiles = []
    logfiles = []
    tpinfo = []
    step = 0
    if prog == 'elastix':
      params = util.ElastixParamFile(str(paramfiles[0]))
      try: useRegAladin = params['UseRegAladin']
      except: useRegAladin = False
      if useRegAladin:
        paramfiles = list(paramfiles) # copy paramfiles before using .pop()
        pf = paramfiles.pop(0)
        elastix = [None] * len(movingimgs)
        for i,movingimg in enumerate(movingimgs):
          elastix[i] = Elastix().setInput(
            fixed = fixedimg,
            moving = movingimg,
            paramfiles = [pf]
          )
          
        # call reg_aladin and use resulting affine matrix as input to elastix
        for m,fmask in enumerate(fixedmasks):
          regaladin = [None] * len(movingimgs)
          tp = [None] * len(movingimgs)
          mmasks = movingmasks[m] if movingmasks else None
          for i,movingimg in enumerate(movingimgs):
            regaladin[i] = NiftyReg().setInput(
              fixed = fixedimg,
              moving = movingimg,
              outfile = self.tempfile('regaladin{}_{}.txt'.format(m,util.basename(movingimg))),
              fixedmask = fmask,
              inaff = regaladin[i].requestOutput('outfile') if regaladin[i] else None
            )
            if mmasks:
              regaladin[i].setInput(movingmask = mmasks[i])
            regaladin2elastix = RegAladin2Elastix().setInput(
              tpbase = elastix[i].requestOutput('tp0'),
              afffile = regaladin[i].requestOutput('outfile')
            )
            # create elastix transformation file based on reg_aladin output
            tp[i] = regaladin2elastix.requestOutput('tpfile')
          tpfiles.append(FancyList(tp))
          tpinfo.append({
           'paramfile': pf,
           'fmask': fmask,
           'mmasks': mmasks
          })
      
      if len(paramfiles) > 0:
        # Elastix uses only the final mask (default: last available mask)
        fmask = fixedmasks[finalmask]
        if bgfrommask is not None:
          fixedimg = BackgroundFromMask().setInput(fixedimg,fixedmasks[bgfrommask],self.tempfile('fixed_bgfrommask.nii.gz')).requestOutput('outfile')
        mmasks = movingmasks[finalmask] if movingmasks else None
        if movingmasks and bgfrommask is not None:
          movingimgs = list(movingimgs)
          for i,movingimg in enumerate(movingimgs):            
            movingimgs[i] = BackgroundFromMask().setInput(movingimg,movingmasks[bgfrommask][i],self.tempfile('moving{}_bgfrommask.nii.gz'.format(i))).requestOutput('outfile')
        elastix = [None] * len(movingimgs)
        for i,movingimg in enumerate(movingimgs):            
          elastix[i] = Elastix().setInput(
            fixed = fixedimg,
            moving = movingimg,
            paramfiles = paramfiles,
            fixedmask = fmask
          )
          if mmasks:
            elastix[i].setInput(movingmask = mmasks[i])
          if len(tpfiles):
            elastix[i].setInput(initialtransform = tpfiles[-1][i])
          logfiles.append(elastix[i].requestOutput('log'))
        for p,pf in enumerate(paramfiles):
          tp = [None] * len(movingimgs)
          for i,movingimg in enumerate(movingimgs):
            tp[i] = elastix[i].requestOutput('tp{}'.format(p))
          tpfiles.append(FancyList(tp))
          tpinfo.append({
           'paramfile': pf,
           'fmask': fmask,
           'mmasks': mmasks
          })
          
    elif prog == 'niftyreg':
      # includes support for multiple masks
      for m,fmask in enumerate(fixedmasks):
        tp = [None] * len(movingimgs)
        mmasks = movingmasks[m] if movingmasks else None
        regaladin = [None] * len(movingimgs)
        for i,movingimg in enumerate(movingimgs):
          regaladin[i] = NiftyReg().setInput(
            fixed = fixedimg,
            moving = movingimg,
            outfile = self.tempfile('regaladin{}_{}.txt'.format(m,util.basename(movingimg))),
            fixedmask = fmask,
            inaff = regaladin[i].requestOutput('outfile') if regaladin[i] else None
          ) 
          if mmasks:
            regaladin.setInput(movingmask = mmasks[i])
          tp[i]  = regaladin[i].requestOutput('outfile')
        tpfiles.append(FancyList(tp))
        tpinfo.append({
         'paramfile': '[regaladin]',
         'fmask': fmask,
         'mmasks': mmasks
        })
    fancyLog(FancyList(tpfiles),'tpfiles')
    return FancyDict(
      tpfiles = FancyList(tpfiles),
      tpinfo = tpinfo,
      logfiles = FancyList(logfiles)
    )
#endclass

if __name__ == '__main__':
  print('This is not a command line utility.')
