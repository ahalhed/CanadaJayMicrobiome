# CanadaJayMicrobiome
Repository for scripts associated with Alicia's Canada jay microbiome project.

## What can you find in this directory?

### Script files
There is a single QIIME2 script in this directory, labelled *gj-q2.sh*. This particular script was run in chunks as individual batch jobs, as some steps were quite long running. These individual jobs are stored in a separate private GitHub repository. The resulting QIIME2 visualizations and artifacts are stored in a directory on the GRAHAM cluster of SHARCNET (ComputeCanada).

The remaining files in the scripts folder are a combination of shell and R scripts; corresponding scripts have the same naming but different extensions. All of these files, with the exception *gj-q2.sh*, are post-QIIME2 analysis in R. 

### Plots
In addition to STDOUT, discussed below, most scripts output a variety of figures corresponding to different steps of the analysis. The labelling of these files is consistent with the type of plot, microbial community (full (no added label), core, or rare), and grid/year combination.

### Output files
Files labelled with the ".out" extension are text files containing the STDOUT from all SHARCNET batch jobs. These files are named according to the grid/year combination whose spatial analysis STDOUT is contained within the file and for the SHARCNET job number.

### Data Files
These files are neither plots nor STDOUT, but contain data/information relevant to the analyses.
