# MultiAtlasRegistration
Multi Atlas Registration pipeline for estimating brain areas in MRI images on the basis of a set of manually segmented examples. This pipeline is one of the outcomes of the 'Biomarker Boosting through eScience' project, funded by the Netherlands eScience Center under grant 027.011.304.

# Installation
The fancypipe folder in this repository is a symbolic link.
Clone https://github.com/rbakker/FancyPipe.git, and copy/link FancyPipe/fancypipe to MultiAtlasRegistration/fancypipe.

# Running the pipeline
The pipeline takes as input a folder that contains one or more target MRI scans in DICOM or NIFTI format, plus a folder that contains a set of atlases. These consist of a label file with manually segmented brain regions, and a corresponding MRI. Atlases are assumed to be preprocessed. The pipeline preprocesses the target MRI scans and nonlinearly registers each atlas to it in several steps, first without any mask to get a rough registration, and then with a brainmask that helps the registration to focus. Both the atlas and target MRI datasets are histogram equalized.
The output of the pipeline is a set of probabilistic atlases, one for each labelled region. They can be thresholded to generate segmented regions, or they can be used as input to further processing steps that use additional constraints to improve the segmentation.

Example:
First go to the folder ./MultiAtlasRegistration/main
Then invoke the pipeline by:
python registeratlases.py -config ../config/hippocampus.xml -atlasdir="/path/to/my/atlasdir" -paramfiles='["../config/Par0000affine.txt"]' -scandir="/path/to/my/scandir"

