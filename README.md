# Sensory Networks Functional Connectivity
This repository contains code used in the following publication: 
- Possidente, T.\*, Tripathi, V.\*, McGuire, J. T. & Somers, D. C. (in press) Interactions between sensory-biased and supramodal working memory networks in the human cerebral cortex. <i>Communications Biology</i>.

This repository contains code for preprocessing, analysing, and making figures for the sensory networks functional connectivity project. 

Running most of this code requires access to the data on the SCC.

Organization:
- ROI_search_spaces_fsaverage: contains all ROI search spaces used in the analyses as freesurfer label files on the fsaverage surface. Note that while many ROI search spaces overlap, no ROIs overlap at the individual-participant level.
- scripts: contains key scripts for file i/o and organization, file conversion, preprocessing, ROI creation, and analysis
- functions: key functions used by scripts
