# Computer code for `Gaussian process aided function comaprison using noisy scattered data`

## Description: 
This repository contains `MATLAB` code for reproducing the results in the paper: "Gaussian process aided function comaprison using noisy scattered data."

## Folders:
1. **Algorithm** - Contains the implementation for the funGP algorithm.
2. **ExperimentsSourceFiles** - Contains the source files for all the simulation experiments and comparison studies.
3. **WindApplicationSourceFiles** - Contains the source files for the turbine comparison case study as well as the turbine datasets under the sub-folder `Datasets`.

## Main Files:
The relevant source files for each of the figures and tables in the paper has been wrapped in the main files with their name corresponding to the tables and figures in the paper. For examples, `Table1.m` calls all the relevant source files to generate the results in Table 1 in the paper. 

## Execution:
The main files can be executed directly from the main folder to generate result for a particular table or figure in the `MATLAB` console. The results are displayed on the `MATLAB` console after the execution of the code. The relevant  `path` has been added to each main file. If one wishes to execute the files from some other folder, one must add the relevant `paths` using the `MATLAB` command: `addpath(<path to source code>)`.  

We do not recommend running the main files for the tables sequentially as each of the lines in the main files can take hours to run. We recommend running each line of the main files parallely on a computing cluster for efficient and timely computation.

## Output:
The output for each of the main files would be displayed on the `MATLAB` console.


### Log Files: 
Each of the source file generates a log file that stores the progress of the computation. Whenever a source file is called, the name of its corresponding log file is printed on the console for reference. 



