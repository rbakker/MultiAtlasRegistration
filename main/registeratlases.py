import sys
sys.path.append('../fancypipe')
sys.path.append('../lib')
from fancypipe import *
import util
import os.path as op
import nibabel,numpy,json
import Image

from preprocess import PreProcess
from multiregistration import MultiRegistration
from multitransformation import MultiTransformation
from createprobabilitymaps import CreateProbabilityMaps
from calculate_volumes_mean import CalculateVolume, CalculateVolumesMean
from niftiinspect import NiftiInspect

def sortedAtlasesFromDir(atlasdir,
  scaninclude,scanexclude,
  labelinclude,labelexclude,
  headmaskinclude=None,headmaskexclude=None,
  brainmaskinclude=None,brainmaskexclude=None
):
  scanfiles = util.dir2files(atlasdir,include=scaninclude,exclude=scanexclude)
  basenames_scan = util.files2basenames(scanfiles, pattern=scaninclude)
  labelfiles = util.dir2files(atlasdir,include=labelinclude,exclude=labelexclude)
  basenames_label = util.files2basenames(labelfiles, pattern=labelinclude)
  diff = set(basenames_scan) ^ set(basenames_label)
  if diff:
    raise AssertionError('Some files in "{}" do not have a matching scan or label file:\n{}\nScan files: {}\nLabel files: {}\nDifference: {}'.format(atlasdir,diff,scanfiles,labelfiles,diff))
  # sort by basenames
  _,scanfiles = zip(*sorted(zip(basenames_scan,scanfiles)))
  _,labelfiles = zip(*sorted(zip(basenames_label,labelfiles)))
  if headmaskinclude:
    headmaskfiles = util.dir2files(atlasdir, include=headmaskinclude,exclude=headmaskexclude);
    basenames_headmask = util.files2basenames(headmaskfiles, pattern=headmaskinclude)
    diff = set(basenames_scan) ^ set(basenames_headmask)
    if diff:
      raise AssertionError('Some files in "{}" do not have a matching headmask file:\n{}\nScan files: {}\nMask files: {}'.format(atlasdir,diff,scanfiles,headmaskfiles))
    _,headmaskfiles = zip(*sorted(zip(basenames_headmask,headmaskfiles)))
  else:
    headmaskfiles = None
  if brainmaskinclude:
    brainmaskfiles = util.dir2files(atlasdir, include=brainmaskinclude,exclude=brainmaskexclude);
    basenames_brainmask = util.files2basenames(brainmaskfiles, pattern=brainmaskinclude)
    diff = set(basenames_scan) ^ set(basenames_brainmask)
    if diff:
      raise AssertionError('Some files in "{}" do not have a matching brainmask file:\n{}\nScan files: {}\nMask files: {}'.format(atlasdir,diff,scanfiles,brainmaskfiles))
    _,brainmaskfiles = zip(*sorted(zip(basenames_brainmask,brainmaskfiles)))
  else :
    brainmaskfiles = None
  return dict(
    scanfiles  = scanfiles,
    labelfiles = labelfiles,
    headmaskfiles = headmaskfiles,
    brainmaskfiles = brainmaskfiles
  )

class MultiRegistrationInspect(FancyTask):
  @staticmethod
  def basename(fname):
    return op.basename(fname) if fname else '-'
  
  def main(self, targetImage,targetMask,deformedAtlasImages,deformedAtlasLabels=None,deformationParams=None,atlasSelect=False):
    layers = [{
      'file':targetImage,
      'pctile':[0.01,99.99],
      'colormap':'alpha-#FA0'
    },{
      'file':targetMask,
      'colormap':'#000,#FF0'
    }]
    for s,step in enumerate(deformedAtlasImages):
      if isinstance(step,(tuple,list)):
        for i,atl in enumerate(step):
          if atlasSelect is not False and i != atlasSelect: continue
          layers.append({
            'file':atl,
            'pctile':[0.01,99.99],
            'colormap':'alpha-#0AF',
          })
          if deformationParams:
            layers[-1]['title'] = 'scan:{},tp:{},mask:{}'.format(i,self.basename(deformationParams[s]['paramfile']),self.basename(deformationParams[s]['fmask']))
          if deformedAtlasLabels:
            layers.append({
              'file':deformedAtlasLabels[s][i],
              'colormap':'#000,#0F0,#F00'
            })
            if deformationParams:
              layers[-1]['title'] = 'label:{},tp:{},mask:{}'.format(i,self.basename(deformationParams[s]['paramfile']),self.basename(deformationParams[s]['fmask']))
      else:
        if atlasSelect is not False: continue
        layers.append({
          'file':step,
          'pctile':[0.01,99.99],
          'colormap':'alpha-#0AF',
        })
        if deformationParams:
          layers[-1]['title'] = 'scan:{},tp:{},mask:{}'.format(s,self.basename(deformationParams[s]['paramfile']),self.basename(deformationParams[s]['fmask']))
        if deformedAtlasLabels:
          layers.append({
            'file':deformedAtlasLabels[s][i],
            'colormap':'#000,#0F0,#F00'
          })
          if deformationParams:
            layers[-1]['title'] = 'label:{},tp:{},mask:{}'.format(s,self.basename(deformationParams[s]['paramfile']),self.basename(deformationParams[s]['fmask']))
          
    niftiInspect = NiftiInspect().setInput(
      layers = layers,          
      htmlfile = self.tempfile('inspectregistration.html')
    )
    return FancyDict(
      htmlfile = niftiInspect.requestOutput('out')
    )
#endclass


class AverageHistogram(FancyTask):
  def main(self,targetImages,regionMask,numBins=64):
    if not isinstance(targetImages,(list,tuple)):
      targetImages = [targetImages]
    nii = nibabel.load(niiMask)
    mask = nii.get_data().astype(bool)
    allData = []
    for tgt in targetImages:
      nii = nibabel.load(niiScan)
      data = nii.get_data()
      all_data.append(data[mask])
    allData = numpy.array(allData)
    mn = allData.min()
    mx = allData.max()
    allData = ((allData-mn)/(mx-mn)*(numBins-1e-8)).astype(numpy.uint8)
    # fit kernel density estimate
    kde_mask = stats.gaussian_kde(mask)
    ind = numpy.linspace(0,numBins-1,numBins)
    pdf_mask = kde_mask.evaluate(ind)
    histogram_mask = [numpy.count_nonzero(allData==b) for b in range(0,numBins)]
    return FancyArgs(bar_histogram=histogram_mask,kde_histogram=pdf_mask,kde_fit=kde_mask,bin_edges=numpy.linspace(mn,mx,numBins+1))
#endtask


class JointHistograms(FancyTask):
  def main(self,targetImage,targetMask,deformedImages,numBins):
    def loaddata(niiFile,numBins,pctile=0.1):
      nii = nibabel.load(niiFile)
      data = nii.get_data()
      mn = numpy.percentile(data,pctile)
      data[data<mn] = 0
      mx = numpy.percentile(data,100.0-pctile)
      data[data>mx] = mx
      return ((data-mn)/(mx-mn)*(numBins-0.01)).astype(numpy.uint8)

    targetMask = nibabel.load(targetMask).get_data()
    histograms = []
    MI = []
    data1 = loaddata(targetImage,numBins)
    for s,step in enumerate(deformedImages):
      MI_s = []
      for i,atl in enumerate(step):        
        imgfile = self.tempfile('atl{}_step{}.png'.format(i,s))
        data2 = loaddata(atl,numBins)
        hh = numpy.zeros([numBins,numBins],numpy.float64)
        for b in range(0,numBins):
          idx = numpy.logical_and(targetMask,data1==b)
          if idx.size:
            for v in data2[idx]:
              hh[b,v] += 1
        hh = hh/hh.sum()
        hX = hh.sum(axis=0)
        hY = hh.sum(axis=1)
        hYX = numpy.outer(hY,hX)
        lt0 = hh>0
        mutInf = hh[lt0]*numpy.log2(hh[lt0]/hYX[lt0])
        MI_s.append(mutInf.sum())
        mx = hh.max()
        mn = hh.min()
        hh = (hh-mn)/(mx-mn)*(256-0.01)
        im = Image.fromarray(hh.astype(numpy.uint8))
        imgfile = self.tempfile('atl{}_step{}.png'.format(i,s))
        im.save(imgfile)
        histograms.append(imgfile)
      MI.append(MI_s)
    return FancyDict(histograms=histograms,MI=MI)
#endtask


class RegisterAtlases(FancyTask):
  title = 'Register multiple atlases to a target and create probability maps'
  description = None
  inputs = dict(
    scandir = {'type':assertDir,
      'help':'Directory that contains input scan files.'
    },
    scanselect = {'type':assertMultiSelect, 'default':None,
      'help':'Indices of sorted scans to read from <scandir> (None for all).'
    },
    atlasdir = {'type':assertDir,
      'help':'Directory that contains input atlas files (scan + label)'
    },
    atlasselect = {'type':assertMultiSelect, 'default':None,
      'help':'Indices of (sorted) atlases to read from <atlasdir> (None for all).'
    },
    paramfiles = {'type':assertList,
      'help':'List of Elastix registration parameter files.'
    },
    scaninclude = {'type':str, 'default':'(.*?)\.nii(\.gz)?$',
      'help':'Regular expression to include scan files from <scandir> and <atlasdir>. First (group) must contain basename.'
    },
    scanexclude = {'type':str, 'default':'-[a-z]+\.nii(\.gz)?$',
      'help':'Regular expression to exclude scan files from <scandir> and <atlasdir>.'
    },
    labelinclude = {'type':str, 'default':'(.*?)-([a-z]+)\.nii(\.gz)?$',
      'help':'Regular expression to include label files from <atlasdir>. First (group) must contain basename.'
    },
    labelexclude = {'type':str, 'default':None,
      'help':'Regular expression to exclude label files from <atlasdir>.'
    },
    mnit1 = {'type':assertFile,
      'help':'MN152 T1, used if no mask is specified.'
    },
    mnimask = {'type':assertFile,
      'help':'MNI52 mask, used if no mask is specified.'
    },
    prog = {'default':'elastix',
      'help':'Which program to use for the registration part.'
    },
    headmaskinclude = {'type':str, 'default':'(.*)\.mnioverlap\.nii(\.gz)?$',
      'help':'Regular expression to include headmask files from <scandir> and <atlasdir>.'
    },
    headmaskexclude = {'type':str, 'default':None,
      'help':'Regular expression to exclude headmask files from <scandir> and <atlasdir>.'
    },
    brainmaskinclude = {'type':str, 'default':None,
      'help':'Regular expression to include brainmask files from <scandir> and <atlasdir>.'
    },
    brainmaskexclude = {'type':str, 'default':None,
      'help':'Regular expression to exclude brainmask files from <scandir> and <atlasdir>.'
    },
    finalmask = {'type':int, 'default':-2,
      'help':'Mask to use at final stage, index into masks array. Default: -2'
    },
    bgfrommask = {'type':assertType(int,{'':None}), 'default':-1,
      'help':'Mask to use as background (everything outside mask) at final stage, index into masks array.'
    }
  )
  
  # if the target image is among the training scans, then leave it out of the training scans
  def atlasesWithoutTarget(self,targetImage,atlasImages,atlasLabels,atlasHeadmasks=None,atlasBrainmasks=None):
    try:
      i = atlasImages.index(targetImage)
      print('The target image is part of the atlas set.')
      atlasImages = list(atlasImages)
      atlasImages.pop(i)
      atlasLabels = list(atlasLabels)
      atlasLabels.pop(i)
      if atlasHeadmasks: 
        atlasHeadmasks = list(atlasHeadmasks)
        atlasHeadmasks.pop(i)
      if atlasBrainmasks: 
        atlasBrainmasks = list(atlasBrainmasks)
        atlasBrainmasks.pop(i)
    except ValueError:
      pass
    return (atlasImages,atlasLabels,atlasHeadmasks,atlasBrainmasks)
    
  def main(self,scandir,scanselect,atlasdir,atlasselect,paramfiles,
           scaninclude,scanexclude,labelinclude,labelexclude,mnit1,mnimask,prog,
           headmaskinclude,headmaskexclude,brainmaskinclude,brainmaskexclude,finalmask,bgfrommask):
    targets = sortedAtlasesFromDir(scandir,scaninclude,scanexclude,labelinclude,labelexclude,headmaskinclude,headmaskexclude,brainmaskinclude,brainmaskexclude)
    if scanselect is not None: 
      targets['scanfiles'] = [targets['scanfiles'][int(i)] for i in scanselect]
      targets['labelfiles'] = [targets['labelfiles'][int(i)] for i in scanselect]
      if targets['headmaskfiles']:
        targets['headmaskfiles'] = [targets['headmaskfiles'][int(i)] for i in scanselect]
      if targets['brainmaskfiles']:
        targets['brainmaskfiles'] = [targets['brainmaskfiles'][int(i)] for i in scanselect]
    atlases = sortedAtlasesFromDir(atlasdir,scaninclude,scanexclude,labelinclude,labelexclude,headmaskinclude,headmaskexclude,brainmaskinclude,brainmaskexclude)
    if atlasselect is not None: 
      atlases['scanfiles'] = [atlases['scanfiles'][int(i)] for i in atlasselect]
      atlases['labelfiles'] = [atlases['labelfiles'][int(i)] for i in atlasselect]
      if atlases['headmaskfiles']:
        atlases['headmaskfiles'] = [atlases['headmaskfiles'][int(i)] for i in atlasselect]
      if atlases['brainmaskfiles']:
        atlases['brainmaskfiles'] = [atlases['brainmaskfiles'][int(i)] for i in atlasselect]

    targetImages = targets['scanfiles']
    targetLabels = targets['labelfiles']
    targetHeadmasks = targets['headmaskfiles']
    targetBrainmasks = targets['brainmaskfiles']
    atlasImages = atlases['scanfiles']
    atlasLabels = atlases['labelfiles']
    atlasHeadmasks = atlases['headmaskfiles']
    atlasBrainmasks = atlases['brainmaskfiles']
  
    labelnames = ['background','hippoL','hippoR']
    tpfiles = []
    volumesMean_pre = []
    volumesMedian_pre = []
    for i,img in enumerate(targetImages):
      atlasImages_i,atlasLabels_i,atlasHeadmasks_i,atlasBrainmasks_i = \
        self.atlasesWithoutTarget(img,atlasImages,atlasLabels,atlasHeadmasks,atlasBrainmasks)
        
      # calculate volume of atlases before registration
      print('Calculating the volume of the labels')
      calcvolume_prereg = CalculateVolumesMean().setInput(
        infiles = atlasLabels_i,
        labelnames = labelnames
      )
      volumesMean_pre.append( calcvolume_prereg.requestOutput('meanVolume') )
      volumesMedian_pre.append( calcvolume_prereg.requestOutput('medianVolume') )

      print('Processing scan "{}"'.format(img))
      preproc = PreProcess().setInput(
        inp        = img,
        out        = self.tempfile('scan{}_prep.nii.gz'.format(i)),
        outmask    = self.tempfile('scan{}_mask.nii.gz'.format(i)),
        mnit1      = mnit1,
        mnimask    = mnimask
      )
      print('Registering scan "{}"'.format(img))
      print('WARNING: ASSUMING THAT SCAN IS ALREADY PREPROCESSED!')
      scanmasks = []
      atlasmasks = []
      if targetHeadmasks: 
        scanmasks.append(targetHeadmasks[i])
        atlasmasks.append(atlasHeadmasks_i)
      if targetBrainmasks: 
        scanmasks.append(targetBrainmasks[i])
        atlasmasks.append(atlasBrainmasks_i)
      # assume a preprocessed image
      multireg = MultiRegistration().setInput(
        fixedimg = img,
        movingimgs = atlasImages_i,
        paramfiles = paramfiles,
        fixedmasks = scanmasks,
        movingmasks = atlasmasks,
        prog = prog,
        finalmask = finalmask,
        bgfrommask = bgfrommask
      )
      tpfiles.append( multireg.requestOutput('tpfiles') )

    # transform images and labels
    deformedAtlasLabels = []
    deformedAtlasImages = []
    transformationLogsImages = []
    transformationLogsLabels = []
    for i,img in enumerate(targetImages):
      atlasImages_i,atlasLabels_i,atlasHeadmasks_i,atlasBrainmasks_i = \
        self.atlasesWithoutTarget(img,atlasImages,atlasLabels,atlasHeadmasks,atlasBrainmasks)
      deformedLabels,logs = MultiTransformation().setInput(
        infiles = atlasLabels_i,
        tpfiles = tpfiles[i][-1],
        dtype = numpy.dtype('int16'),
        interp = 0,
        prog = prog,
        ref_img = img
      ).linkOutput('outfiles','logfiles')
      deformedAtlasLabels.append( deformedLabels ) 
      transformationLogsLabels.append( logs )

      deformedImages,logs = MultiTransformation().setInput(
        infiles = atlasImages_i,
        tpfiles = tpfiles[i][-1],
        dtype = numpy.dtype('float'),
        interp = 3,
        prog = prog,
        ref_img = img
      ).linkOutput('outfiles','logfiles')
      deformedAtlasImages.append(deformedImages)
      transformationLogsImages.append(logs)
    
      multireginspect = MultiRegistrationInspect().setInput(
        targetImage = img,
        targetMask = targetHeadmasks[i],
        deformedAtlasImages = deformedImages
      )
    
    # create probability maps
    probabilityMaps = []
    histograms = []
    outliers = []
    volumes = []
    for i,img in enumerate(targetImages):
      probmap = CreateProbabilityMaps().setInput(
        infiles = deformedAtlasLabels[i],
        labelnames = labelnames,
        outfilebase = self.tempfile('probmap{}'.format(i)),
        groundtruth = targetLabels[i]
      )
      probabilityMaps.append( probmap.requestOutput('outfiles') )
      histograms.append( probmap.requestOutput('histograms') )
      #outliers.append( probmap.requestOutput('outliers') )
      volumes.append( probmap.requestOutput('volumes') )
      
    # calculate volumes of the labels after registration
    volumesRef = []
    volumesMean_post = []
    volumesMedian_post = []   
    for i,img in enumerate(targetImages):
      calcvolume_postreg = CalculateVolumesMean().setInput(
        infiles = deformedAtlasLabels[i],
        labelnames = labelnames
      )
      volumesMean_post.append( calcvolume_postreg.requestOutput('meanVolume') )
      volumesMedian_post.append( calcvolume_postreg.requestOutput('medianVolume') )
      
      # calculate the volume of the manual segmentation of the target image
      calcvolume_ref = CalculateVolume().setInput(
        infile = targetLabels[i],
        labelnames = labelnames
      )
      volumesRef.append(calcvolume_ref.requestOutput('volume'))
  
    return FancyDict(
      deformedAtlasLabels = FancyList(deformedAtlasLabels),
      probabilityMaps = FancyList(probabilityMaps),
      histograms = FancyList(histograms),
      outliers = FancyList(outliers),
      volumes = FancyList(volumes),
      volumesMean_pre = FancyList(volumesMean_pre),
      volumesMedian_pre = FancyList(volumesMedian_pre),
      volumesMean_post = FancyList(volumesMean_post),
      volumesMedian_post = FancyList(volumesMedian_post),
      volumesReference = FancyList(volumesRef),
      tpfiles = FancyList(tpfiles),
      transformationLogsLabels = FancyList(transformationLogsLabels),
      htmlfiles = multireginspect.requestOutput('htmlfile')
    )
#endclass


if __name__ == '__main__':
  mainTask = RegisterAtlases.runFromCommandLine()
  with open(mainTask.tempfile('result.json'),'w') as fp:
    json.dump(mainTask.myOutput.jsonify(),fp,indent=2)
