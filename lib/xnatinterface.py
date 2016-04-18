from fancypipe import *
import os,os.path as op
import pyxnat
import zipfile,re

def promptPassword(user):
  import getpass
  s = getpass.getpass('XNAT password for user "{}": '.format(user))
  return FancyPassword(s)
  
class ScansFromProjectToZip(FancyTask):
  title = 'Download (a subset of) scans from an XNAT project and save as zip files.'
  inputs = {
    'server'  : { 'help':'XNAT server URL.' },
    'project' : { 'help':'XNAT project' },
    'user'    : { 'help':'XNAT user with access rights to <project>' },
    'pwd'     : { 'type':assertPassword,'default':promptPassword,'help':'XNAT password for <user>' },
    'scantypes': { 'type':assertList, 'default':None, 
      'help':'Scan types to download (None for ALL)' }
  }
  
  def main(self, server,user,pwd,project,scantypes):
    cachedir = self.tempdir('XNAT')
    include = '{}(.*)\.zip$'.format(re.escape(project))
    zipfiles = { v:op.join(cachedir,v) for v in os.listdir(cachedir) if re.search(include,v) }

    if len(zipfiles) == 0:
      xnat = pyxnat.Interface(server=server, user=user, password=pwd.value, cachedir=cachedir)
      proj = xnat.select.project(project)
      if not proj.exists():
        raise LookupError('Project "{}" does not seem to exist on {} for user {}.'.format(project,server,user))
      subjects = proj.subjects()
      if not subjects:
        print('Project "{}" is empty!'.format(project,server,user))
        
      if isinstance(scantypes,(list,tuple)):
        scantypes = ','.join(scantypes)
      else:
        scantypes = 'ALL'

      for subj in subjects:
        print('{}'.format(subj))

        experiments = subj.experiments().get()
        for m in experiments:
          expmt = subj.experiment(m)
          scans = expmt.scans()
          sub_ses = subj.id()
          if len(experiments)>0: sub_ses += '_{}'.format(expmt.id())
          print("Downloading scantypes {} from project {}, subject {}, session {}".format(scantypes,proj.id(),subj.id(),expmt.id()))
          zipfiles[sub_ses] = scans.download(cachedir, type=scantypes, extract=False)

    return FancyDict(zipfiles)
#endclass


class ScansFromProjectNifti(ScansFromProjectToZip):
  inputs = ScansFromProjectToZip.inputs.copy()
  inputs['outdir'] = { 'help':'Output directory' }

  def main(self, server,user,pwd,project,outdir,scantypes):
    zipfiles = ScansFromProjectToZip.main(self, server,user,pwd,project,scantypes)
    niftifiles = []
    if not op.isdir(outdir):
      os.makedirs(outdir)
    for k,z in zipfiles.items():
      pattern = '{}\_([^_]+)\_([^_]+)\_([^_]+)\_([^_]+)\_scans\_([^_]+)\.zip'.format(re.escape(project))
      matches = re.match(pattern,op.basename(z))
      fzip = zipfile.ZipFile(z,'r')
      destdir = op.join(outdir,'{}_{}'.format(matches.group(2),matches.group(4)))
      if not op.isdir(destdir): os.mkdir(destdir)
      for info in fzip.infolist():
        fname = info.filename
        if fname.endswith(('.nii','.nii.gz')):
          dest = op.join(destdir,op.basename(fname))
          print('Extracting {} from {} to {}.'.format(fname,z,dest))
          data = fzip.read(fname)
          open(dest,'w+b').write(data)
          niftifiles.append(dest)
    return FancyList(niftifiles)
#endclass


if __name__ == '__main__':
  ScansFromProjectNifti.fromCommandLine().run()
