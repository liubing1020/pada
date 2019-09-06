============================================================================
==   Global Sensitivity Analysis (MPSA)                                   ==
==                                                                        ==
==   Author: Bing Liu (liubing@comp.nus.edu.sg)                           ==
==   Organization: National University of Singapore                       ==
==                                                                        ==
============================================================================

This module is for performing MPSA algorithm (Cho et al 2003) to compute
global sensitivity analysis. 

MPSA.c contains the code to compute intermediate results (frequencies). It 
should include the correct FF inference file (e.g. toy_ff0.c) for a model. 
To compile:
   $./make.sh
The input format is as follows:
   $./MPSA [par_id] [par_no] [dim] [threshold] [sampleNo]
 - par_id: the ID of parameter to be tested
 - par_no: the total number of paramaters
 - dim: the number of intervals a parameter can take
 - threshold: a value to distinguish good or bad samples (e.g. the average of 
   objective values
 - sampleNo: the number of samples
Thus, the following line compute the frequency results (stored in Nk2.txt) 
for parameter k2 of the toy model:
   $./MPSA 2 3 5 200.0 100
Note that generating Nkx.txt separately can facilitate the usage of PC cluster 
(e.g. allocate a processing core for each parameter), since it could be 
computational intensive for large sample size. 

After generating the frequency files, ksall.java can be used to compute the 
final global sensitivities:
   $javac ksall.java
   $java ksall [sampleNo] [parNo] [dim]
the global sensitivities will be display in the standard output
the generated mpsa-all.txt is used for plotting the frequency graphs

Example:
   $cd MPSA/toy
   $./make.sh
   $./run.sh
   $javac ksall.java
   $java ksall 100 3 5




