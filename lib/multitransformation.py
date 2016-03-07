from fancypipe import *
import os.path as op
import util
import numpy

class RegResample(FancyExec):
  title="Call <reg_resample>."
  
  def main(self,ref_img,float_img,tpfile,outfile,interp):
    return FancyExec.main(
      self,
      'reg_resample',
      **{
        '-flo'   : str(float_img),
        '-ref'   : str(ref_img),
        '-trans' : tpfile,
        '-inter' : str(interp),
        '-res'   : FancyOutputFile(outfile)
      }
    )
#endclass
  

class Transformix(FancyExec):
  title= "Call <transformix>."

  def outputFile(self,key):
    files = dict(
      outfile = self.tempfile('result.nii.gz'),
      logfile = self.tempfile('transformix.log')
    )
    if key in files: return files[key]

  def main(self,infile,tpfile,dtype,interp):
    newtp = util.ElastixParamFile(tpfile)
    newtp['FinalBSplineInterpolationOrder'] = interp
    newtp['CompressResultImage'] = '"true"'
    newtp['ResultImageFormat'] = '"nii.gz"'
    pixeltypes = {'uint8':'unsigned char','int8':'char','uint16':'unsigned short','int16':'short','float32':'float','float64':'double'}
    newtp['ResultImagePixelType'] = '"{}"'.format(pixeltypes[str(dtype)])
    (basename,_) = op.splitext(op.basename(tpfile))
    newtpfile = self.tempfile('{}_{}_{}.txt'.format(basename,dtype,interp))
    newtp.saveAs(newtpfile, overwrite=True)
    FancyExec.main(
      self,
      'transformix',      
      **{
        '-in': infile,
        '-tp': newtpfile,
        '-out': self.tempdir(),
      }
    )

    return FancyDict(
      outfile = self.outputFile('outfile'),
      logfile = self.outputFile('logfile')
    )
#endclass


class MultiTransformation(FancyTask):
  title = "Transform an image using previously computed transformation parameters."
  inputs = {
    'infiles': {'type': assertList},
    'tpfiles': {'type': assertList},
    'dtype': {'type': str, 'help':'Datatype of the result file (Elastix only).'},
    'interp': {'type': int, 'help':'Interpolation order, 0 = nearest neighbor, 1 = linear, 3 = cubic.'},
    'prog': {'help':'Program: elastix or niftyreg.'},
    'ref_img': {'type': assertFile, 'help':'Reference image to infer sampling grid from (NiftyReg only).'}
  }

  def main(self,infiles,tpfiles,dtype,interp,prog,ref_img):
    outfiles = FancyList()
    logfiles = FancyList()
    for i,infile in enumerate(infiles):
      tpfile = tpfiles[i]
      if prog == 'elastix':
        transformix = Transformix().setInput(
          infile = infile,
          tpfile = tpfile,
          dtype  = dtype,
          interp = interp
        )
        outfiles.append(transformix.requestOutput('outfile'))
        logfiles.append(transformix.requestOutput('logfile'))
      elif prog == 'niftyreg':
        regresample = RegResample().setInput(
          float_img = infile,
          ref_img = ref_img,
          tpfile = tpfile,
          outfile = self.tempfile('regresample_{}_{}.nii.gz'.format(util.basename(infile),i)),
          interp = interp
        )
        outfiles.append(regresample.requestOutput('outfile'))
    return FancyDict(
      outfiles = outfiles,
      logfiles = logfiles
    )
#endclass


if __name__ == '__main__':
  MultiTransformation.fromCommandLine().run()
