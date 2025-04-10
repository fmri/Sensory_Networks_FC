# Sensory Networks Functional Connectivity
Tom Possidente 
tposside@bu.edu
Somers Lab - BU

This repository contains code for preprocessing, analysing, and making figures for the sensory networks functional connectivity project. 

Running most of this code requires access to the data on the SCC.

Organization:
- ROI_search_spaces_fsaverage: contains all ROI search spaces used in the analyses as freesurfer label files on the fsaverage surface. Note that while many ROI search spaces overlap, no ROIs overlap at the individual-participant level.
- scripts: contains key scripts for file i/o and organization, file conversion, preprocessing, ROI creation, and analysis
- archive_code: old scripts/functions no longer used/useful
- functions: key functions used by scripts
