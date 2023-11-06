# igNeurons

## overview 
These are scripts used for the paper: "A human-specific microRNA controls the timing of excitatory synaptogenesis."
Please contact me directly if you have any questions to one of the files or markdowns.

## requirements
R packages needed to run the scripts are listed in the "library" chunks. CellProfiler 4 (https://cellprofiler.org/)
was used for the image analysis. I am happy to help out if any of the scripts is not working with a specific package
version. 

## instructions
Generally the data and scripts are sorted by experiment. They are further divided in "data" folders containing raw data,
"R" folders with scripts used for data processing and statistical analysis and "Figure" folders to generate the figures 
used in the paper. The files in the "Figure Output" folders were generated with these scripts and it should be possible to 
reproduce them using the provided information. In case something is not working or there are questions you can contact me 
or create an issue.

## image analysis
The cell profiler pipelines can be found in the respective experimental folders. Example images for each experiment can 
be provided upon request since the file size of individual tif images is too big for a convenient handling in the git repository. A single
example image to illustrate the concept of the synaptogenesis cell profiler pipelines can be found in "Image_Analysis/pLNA_experiment/
Synapses/Example_Image". The cell profiler output files can be imported in R and should contain the relevant column titles used for further 
analysis or figure generation. 

## sequencing- and proteomics data
The general structure for the omics data is the same as described previously. Raw Summarized Experiment files can be found in the respective "data" folders,
the "..._R" folders contain relevant scripts for the analysis and the "Figure" folders the scripts to generate figures used in the paper. The igNeuron time course
data and scripts can be found in "TimeCourse_GRCh38.p10_Ens91" and the pLNA-sequencing data in "pLNA_GRCh38.p13_Ens107".

## apps
Custom developed Shiny Apps for Ca-Imaging analysis (https://ethz-ins.org/SpikeIt/), microRNA binding site analysis (https://ethz-ins.org/scanMiR/)
or microRNA enrichment analysis (https://ethz-ins.org/enrichMiR/) are freely available online.

## igNeuron time course data
The different igNeuron time course datasets can be accessed at https://ethz-ins.org/igNeuronsTimeCourse/.