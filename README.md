# CanadaJayMicrobiome
Repository for scripts associated with Alicia's Canada jay microbiome project.

## What can you find in this directory?

### scripts
There is a single QIIME2 script in this directory, labelled *gj-q2.sh*. This particular script was run as individual batch jobs, as some steps are long running. These individual jobs are stored in a subdirectory of *scripts* called *batchJobs*. The resulting QIIME2 visualizations and artifacts are stored in a directory above the repository on Graham (a Compute Canada cluster).

The remaining files in the scripts folder are a combination of shell and R scripts; corresponding scripts have the same naming but different extensions. All files in this directory, with the exception *gj-q2.sh* and scripts in *batchJobs*, are post-QIIME2 analysis in R. 

### plots
In addition to STDOUT, discussed below, some scripts output figures corresponding to different steps of the analysis. The labelling of these files is consistent with the hypothesis being investigated and the type of figure.

### outputs
Files labelled with the ".out" extension are text files containing the STDOUT from all SHARCNET batch jobs. These files are named according to the analyses the STDOUT is associated with, and for the SHARCNET job number.

### data
These files are neither plots nor STDOUT, but contain data/information relevant to the analyses. These are small data files only. Large data files are stored a directory above this GitHub repository, stored locally (backed up to OneDrive) and in a projects directory on Graham.
