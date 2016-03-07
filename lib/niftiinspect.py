# Rembrandt Bakker, June 2014
from fancypipe import *
import sys,os,os.path as op
import numpy,scipy,re,json,math
import niitools

class NiftiInspect(FancyTask):
  title = 'Create html-page with selected slices from multiple nifti files.'
  description = """
    Colormap supports various formats:
    1. comma separated list of rgb-hex values: #RGB1,#RGB2,...
    2. range of rgb-hex values: #RGB1-#RGB2
    3. constant color with range of alpha values: alpha-#RGB
  """
  inputs = {
    'layers': {'type':assertList, 
      'help':"List of layers, each of which much have fields 'file','title','colormap','pctile'"
    },
    'htmlfile': {'type':str,
      'help':'Html output file'
    },
    'imagefolder': {'type':str,
      'default': None, 'help':'Image folder (None to use path/basename_files, with path and basename taken from the html file)'
    },
    'slices_x': {'type':int,
      'default': '0%:10%:100%', 'help':'Slices in the x-dimension, start%:step%:stop%, default 0%:10%:100%'
    },
    'slices_y': {'type':int,
      'default': '0%:10%:100%', 'help':'Slices in the y-dimension, start%:step%:stop%, default 0%:10%:100%'
    },
    'slices_z': {'type':int,
      'default': '0%:10%:100%', 'help':'Slices in the z-dimension, start%:step%:stop%, default 0%:10%:100%'
    }
  }
  
  def main(self,layers,htmlfile,imagefolder,slices_x,slices_y,slices_z):
    # Basic argument checking
    for lr in layers: assertFile(lr["file"])

    sliceRangePct = [[],[],[]]
    for d in [0,1,2]:
      dim = ['x','y','z'][d]
      s = locals()['slices_'+dim]
      if s:
        s = s.split(':');
        sliceRangePct[d] = [int(r.rstrip('%')) for r in s]
      else:
        sliceRangePct[d] = [0,10,100]

    print('Layer specification: {}'.format(layers))

    if op.isdir(htmlfile):
      htmlfile = op.join(htmlfile,'index.html')
    htmlfolder = op.dirname(htmlfile)
    if not imagefolder:
      htmlname,_ = op.splitext(op.basename(htmlfile))
      imagefolder = op.join(htmlfolder,htmlname+'_files')
    if not(op.exists(htmlfolder)):
      os.makedirs(htmlfolder)
      print('Created html output folder "{}".'.format(htmlfolder))
    if not(op.exists(imagefolder)):
      os.makedirs(imagefolder)
      print('Created image output folder "{}".'.format(imagefolder))
    imagefolder = op.realpath(imagefolder)
    scriptdir = op.realpath(op.dirname(__file__))
    parsedLayers = []
    for i,lr in enumerate(layers):
      nifti_src = lr["file"]
      if not nifti_src: continue
      baseName = 'lr{}_'.format(i)+re.sub('(\.nii|\.nii.gz)$','',op.basename(nifti_src))

      import nibabel
      nii = nibabel.load(nifti_src)
      nii = nibabel.as_closest_canonical(nii)
      img = numpy.squeeze(nii.get_data())
      hdr = nii.get_header()
      dims = img.shape
      print('Nifti image loaded, shape {}, data type "{}"'.format(dims,img.dtype))

      if len(dims)==4:
        raise Exception('NIFTI file with RGB color data not supported yet.')

      # apply colormap
      index2rgb = None
      rgbLen = 0
      if "colormap" in lr:
        index2rgb = niitools.parse_colormap(lr["colormap"])
        rgbLen = len(index2rgb[0])

      minmax = [None,None]
      rescale = ('pctile' in lr and lr["pctile"])
      if rescale:
        minmax = niitools.get_limits(img,lr["pctile"])

      fmt = 'png'
      if rescale and rgbLen<4:
        fmt = 'jpg'
                  
      sliceRange = [[],[],[]]
      for d in [0,1,2]:
        dim = ['x','y','z'][d]
        numSlices = dims[d];
        sliceStep = int(math.ceil(sliceRangePct[d][1]*numSlices/100.0))
        sliceStart = int(sliceRangePct[d][0]*(numSlices-1)/100)
        sliceEnd = int(sliceRangePct[d][2]*(numSlices-1)/100)
        sliceRange[d] = [sliceStart,sliceStep,sliceEnd]
        for i in range(sliceStart,sliceEnd+1,sliceStep):
          slice = niitools.get_slice(img,d,i)
                    
          pngFile = baseName+'_{}{:d}.{}'.format(dim,i,fmt)
          if index2rgb:
            slice = niitools.slice2rgb(slice,index2rgb,rescale,minmax[0],minmax[1])

          # Save image to PNG
          scipy.misc.toimage(slice,cmin=0,cmax=255).save(op.join(imagefolder,pngFile))
                        
          if i==sliceStart:
            print('image {}{} saved to png file "{}".'.format(dim,i,pngFile))

      pixdim = hdr['pixdim'][1:4]
      imgsize_mm = [
        round(pixdim[0]*dims[0],1),
        round(pixdim[1]*dims[1],1),
        round(pixdim[2]*dims[2],1)
      ]
      print('Image size in mm {}'.format(imgsize_mm))

      # update parsedLayers
      pl = {
        "name": baseName,
        "ext": fmt,
        "src": nifti_src,
        "imgsize_px": dims,
        "imgsize_mm": imgsize_mm
      }
      if "title" in lr:
        pl["title"] = lr["title"];
      parsedLayers.append(pl);

    inspectFile = '{}/../nifti-tools/nii_inspect.html'.format(scriptdir);
    with open(inspectFile, 'r') as fp:
      html = fp.read()
      html = html.replace(r"var defaultLayers = [];",
        r"var defaultLayers = {};".format(json.dumps(parsedLayers)))                            
      html = html.replace(r"var defaultSliceRange = [];",
        "var defaultSliceRange = {};".format(json.dumps(sliceRange)))
      html = html.replace(r"var imgDir = '';",
        "var imgDir = '{}/';".format(op.relpath(imagefolder,htmlfolder)))
  
    with open(htmlfile, 'w') as fp:
      fp.write(html)
                
    print('HTML viewer saved as "{}"'.format(htmlfile))

    return FancyDict(
      out = htmlfile
    )
#endclass

if __name__ == '__main__':
  NiftiInspect.fromCommandLine().run()
