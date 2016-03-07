# Rembrandt Bakker, June 2014

import os.path as op
import re,json,numpy
import nibabel

def get_limits(img,pctile=False):
    minmax = [ numpy.amin(img), numpy.amax(img) ]
    if pctile: 
        if isinstance(pctile,str):
            pctile = re.split('[,\- ]+',pctile)
            pctile = [float(p) for p in pctile]
        if len(pctile)<1:
            pctile = [0,100]
        elif len(pctile)<2:
            pctile = [0,pctile[0]]
        if pctile[1]<=pctile[0]:
            raise Exception('Max percentile must be larger than min percentile, not {},{}'.format(pctile[0],pctile[1]))
        elif pctile[0]<0:
            raise Exception('Min percentile must be >=0, not {}'.format(pctile[0]))
        elif pctile[1]>100:
            raise Exception('Max percentile must be <=100, not {}'.format(pctile[1]))
        if pctile[0]>0:
            minmax[0] = numpy.percentile(img,pctile[0])
        if pctile[1]>0:
            minmax[1] = numpy.percentile(img,pctile[1])
    return minmax

"""
In the code below, the transpose is needed because images are written with the 
first dimension as rows, second as columns. This must be flipped to align with
the common definition of x- and y axes .
The ::-1 part is a mirroring operation on the y-axis, which is needed because
images are written top to bottom, reversed from the common y-axis direction.
"""
def get_slice(img,dim,i):
    if dim == 0:
        slc = img[i,:,::-1].squeeze();
    elif dim == 1:
        slc = img[:,i,::-1].squeeze();
    elif dim == 2:
        slc = img[:,::-1,i].squeeze();
    else:
        raise Exception('Cannot return a slice for dimension "{}"'.format(dim))
    slc = slc.swapaxes(0,1)
    return slc

def rgb2hex(rgb,prefix='#'):
    return prefix+'{:02X}{:02X}{:02X}'.format(rgb[0],rgb[1],rgb[2])
    
def hex2rgb(v):
    v = v.lstrip('#')
    lv = len(v)
    if lv==3:
        rgb = [int(v[i]+v[i], 16) for i in range(0,lv)];
    elif lv==6:
        rgb = [int(v[i:i+2], 16) for i in range(0,lv,2)];
    else:
        raise Exception('Invalid hex color {}'.format(v))
    return rgb

def hex2rgba(v):
    rgba = hex2rgb(v)

    # transparent background
    if all(v==0 for v in rgba):
        rgba.append(0)
    else: 
        rgba.append(255)
    return rgba

def offsetcolormap(cmap0,offset):
  return {offset+i:rgb for (i,rgb) in enumerate(cmap0)}

def makecolormap(A,numColors,bg=[],bgColor=[0,0,0]):
  numBins = len(A)-1
  perBin = numColors/numBins
  if numBins*perBin<numColors:
      perBin += 1
  cmap = [];
  for j in range(0,perBin):
    for i in range(0,numBins):
      cmap.append([
        int(A[i][0]+(A[i+1][0]-A[i][0])*j*1.0/perBin),
        int(A[i][1]+(A[i+1][1]-A[i][1])*j*1.0/perBin),
        int(A[i][2]+(A[i+1][2]-A[i][2])*j*1.0/perBin),
      ])
  #if len(set(cmap))<len(cmap):
  #  raise Exception('Generated colormap is not invertible')
  for i in bg:
    cmap[bg[i]] = bgColor
  return cmap  

def contrastmap(numColors,bg=[],offset=0):
  A = [
    [  0,  0,255],
    [  0,255,255],
    [  0,255,127],
    [  0,255,  0],
    [127,255,  0],
    [255,255,  0],
    [255,127,  0],
    [255,  0,  0],
    [255,  0,127],
    [255,  0,255],
    [  0,  0,255]
  ]
  return makecolormap(A,numColors,bg)
  
def slice2rgb(slice,index2rgb,rescale,minLevel=None,maxLevel=None):
    if isinstance(index2rgb,dict):
        rgbLen = len(index2rgb[next(iter(index2rgb))])
    else:
        rgbLen = len(index2rgb[0])
        
    shape = slice.shape
    slice = slice.reshape(-1)
    if rescale:
        slice = 255.9999*(slice-minLevel)/(maxLevel-minLevel)
        slice[slice<0] = 0
        slice[slice>255] = 255
        slice = numpy.uint8(slice)
    rgbImg = numpy.zeros(shape=(slice.shape[0],rgbLen), dtype=numpy.uint8)
    for idx in numpy.unique(slice):
        mask = (slice == idx)
        try:
            val = index2rgb[numpy.uint16(idx)]
        except KeyError:
            try:
                val = index2rgb[str(idx)]
            except KeyError:
                raise Exception('Color index {} not in {},'.format(idx,index2rgb.keys()))
        except IndexError:
            raise Exception('Color index {} not in colormap of size {},'.format(idx,len(index2rgb)))

        rgbImg[mask] = val
    return rgbImg.reshape(shape[0],shape[1],rgbLen)

def slicedir_ras(dim):
    return ['saggital','coronal','axial'][dim]

def hexmap2rgbamap(cmap):
    index2rgb = {};
    if isinstance(cmap,list):
        for i,a in enumerate(cmap):
            index2rgb[i] = hex2rgba(a)
    else:
        for i,a in cmap.iteritems():
            index2rgb[i] = hex2rgba(a)
    return index2rgb
    
def parse_colormap(cmap,minmax=None,bg=[0]):
    index2rgb = None
    if cmap:
        matches = re.search('^(#[0-9a-fA-F]+)-(#[0-9a-fA-F]+)$',cmap)
        if matches:
            rgb1 = hex2rgb(matches.group(1))
            rgb2 = hex2rgb(matches.group(2))
            index2rgb = [
                numpy.array([
                    rgb1[0]+i/255.0*(rgb2[0]-rgb1[0]),
                    rgb1[1]+i/255.0*(rgb2[1]-rgb1[1]),
                    rgb1[2]+i/255.0*(rgb2[2]-rgb1[2])
                ],numpy.uint8) for i in range(256)
            ]
        elif cmap == 'auto':
            numColors = 256;
            if minmax:
                numColors = int(minmax[1]-minmax[0]+1)
            index2rgb = contrastmap(numColors,bg,offset=minmax[0])
            if minmax[0]>0: 
                index2rgb = offsetcolormap(index2rgb,minmax[0])
        elif cmap.startswith('alpha'):
            matches = re.search('^alpha-(#[0-9a-fA-F]+)$',cmap)
            if matches:
                rgb = hex2rgb(matches.group(1))
            else:
                rgb = list([255,255,255])
            index2rgb = [[rgb[0],rgb[1],rgb[2],i] for i in range(256)]
        elif cmap[0] == '#':
            index2rgb = hexmap2rgbamap(cmap.split(','))
        elif cmap[0] == '{' or cmap[0] == '[':
            index2rgb = hexmap2rgbamap(json.loads(cmap))
        elif op.exists(cmap):
            with open(cmap,'r') as fp:
                index2rgb = hexmap2rgbamap(json.load(fp))
        else:
            raise Exception('Do not know how to parse colormap "{}".'.format(cmap))
    return index2rgb

def rasLimits(hdr):
    import nibabel
    q = hdr.get_best_affine();
    ornt = nibabel.io_orientation(q)

    #DEBUG HERE
    #verbose('q {}'.format(q))
    #verbose('ornt {}'.format(ornt))
    #print 'hdr {}'.format(hdr)
    dims = hdr.get_data_shape()
    rasLimitsT = numpy.array([
        [-0.5,-0.5,-0.5],
        [dims[0]-0.5,dims[1]-0.5,dims[2]-0.5]
    ])
    rasLimits = nibabel.affines.apply_affine(q,rasLimitsT).T
    #verbose('rasLimits {}'.format(rasLimits))
    for lim in rasLimits:
        if lim[1]<lim[0]: 
            tmp = lim[0]
            lim[0] = lim[1]
            lim[1] = tmp
            
    #verbose('rasLimits {}'.format(rasLimits))
    #if args.out:
    #    with open(args.out, 'w') as fp:
    #        json.dump(rasLimits.tolist(),fp)
    
    return rasLimits.tolist()

def autoBackgroundColor(img):
    sz = img.shape
    rgbMode = (sz[-1]==3 or sz[-1]==4)
    if not rgbMode: raise Exception('Image must be in RGB color mode.')
    if len(sz)==4:
        a = [
            img[0,0,0,:],img[0,0,-1,:],img[0,-1,0,:],img[0,-1,-1,:],
            img[-1,0,0,:],img[-1,0,-1,:],img[-1,-1,0,:],img[-1,-1,-1,:],
        ]
    elif len(sz)==3:
        a = [
            img[0,0,:],
            img[0,-1,:],
            img[-1,0,:],
            img[-1,-1,:]
        ]
    med = numpy.median(a,0)
    bgColor = a[numpy.argmin(numpy.linalg.norm(a-med,1))]
    return bgColor

def imageMask(img,bgColors):
    sz = img.shape
    if hasattr(bgColors[0],"__len__") and sz[-1] == len(bgColors[0]):
        # assume bgColors contain rgb triplets or rbga quadruplets
        mask = numpy.zeros(sz[:-1],dtype=numpy.bool)
        for bg in bgColors:
            submask = numpy.all(img == bg,axis=len(sz)-1)
            submask.sum()   
            mask = numpy.bitwise_or(mask,submask)
    else:
        # assume bgColors contain bgIndices
        mask = numpy.zeros(sz,dtype=numpy.bool)
        for bg in bgColors:
            submask = (img == bg)
            mask = numpy.bitwise_or(mask,submask)
    return numpy.invert(mask)

# img is a 3D volume, q an affine transformation matrix, f is the integer downsampling factor
def downsample(img,q,f):
    if f>1:
        dims = img.shape
        dimg = img[0:dims[0]-f+1:f,:,:]

        # check if image type can support this operation
        try:
            dtype = img.dtype
            limit = numpy.iinfo(dtype).max
            if img.max() > limit/f/f/f:
                dimg = dimg.astype(numpy.float)
                print('Converting image to float, because of dtype limit "{}"'.format(limit))
        except:
            pass
            
        
        for di in range(1,f):
            dimg += img[di:dims[0]-f+1+di:f,:,:]

        img = dimg
        dimg = img[:,0:dims[1]-f+1:f,:]
        for di in range(1,f):
            dimg += img[:,di:dims[1]-f+1+di:f,:]

        img = dimg        
        dimg = img[:,:,0:dims[2]-f+1:f]
        for di in range(1,f):
            dimg += img[:,:,di:dims[2]-f+1+di:f]

        img = dimg/(f*f*f)
        img = img.astype(dtype)

        q[0:3,3] += q[0:3,0:3].dot([f-1,f-1,f-1])/2
        q[0:3,0:3] = q[0:3,0:3]*f
    return (img,q)

def reorient(nii,xixjxk):
    hdr = nii.get_header()
    dims = hdr.get_data_shape()
    dims = [x for x in dims if x > 1]
    print("Dims {}".format(dims))
    
    if len(xixjxk) is not 2*len(dims):
        raise Exception('Permutation string must be twice the length of the data, e.g. "+i+j+k"')
    
    Multiply = {'+':1,'-':-1}
    DirIndex = {'i':0,'j':1,'k':2}
    flip = []
    perm = []
    odd = True
    for ch in xixjxk:
        if odd: 
            flip.append(Multiply[ch])
        else: 
            perm.append(DirIndex[ch])
        odd = not odd
    perm.append(len(dims))
    flip.append(1)

    print('Permute: {}'.format(perm))
    print('Flip: {}'.format(flip))

    qf,qfc = hdr.get_qform(coded=True)
    sf,sfc = hdr.get_sform(coded=True)

    print('Qform before:\n{}'.format(qf))
    if (qf is None):
        qf = hdr.get_qform(coded=False)
        print('Qform reset to:\n{}'.format(qf))
    if (qf is not None):
        qf = qf[perm,:]
        print('qf.T {}'.format(qf.T))
        qf = (flip*qf.T).T
        hdr.set_qform(qf)
        print('Qform after:\n{}'.format(qf))

    print('Sform before:\n{}'.format(sf))
    if (sf is not None):
        sf = sf[perm,:]
        sf = (flip*sf.T).T
        hdr.set_sform(sf)
        print('Sform after:\n{}'.format(sf))
    
    return nibabel.nifti1.Nifti1Image(nii.get_data(),hdr.get_best_affine(),hdr)
    
