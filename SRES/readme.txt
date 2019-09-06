============================================================================
==   A Brief Introdcution of Using SRES                                 ==
==                                                                        ==
==   Author: Bing Liu (liubing@comp.nus.edu.sg)                           ==
==   Organization: National University of Singapore                       ==
==                                                                        ==
==   Note: the code is based on libSRES by Xinglai Ji (jix1@ornl.gov)     ==
============================================================================

FFSRES is a discretized version of libSRES which is used to perform SRES algorithm 
to do parameter estimation. The changes were made on ./source/ESES.c in order to
return solutions only in terms of integer vectors. 

To compile: 
  $make

To use FFSERS to do parameter estimation for the DBN approximation. We compute 
the fitness score using FF inference.
Example:
   $cd FFSRES/toy
   The sample code for estimating parameters is called toySRES.C, see comments 
   starting with '// [SPA]:'.

   Two required files are:
   - toy_ff0.c // FF inference code for toy model (see ../SPA/src/toy_ff0.c)
   - toy_data.txt // synthetic data

   To compile and execute: 
   $make
   $./toySRES

   The results look like:
   gen=10,dt=1,bestgen=10,bestfitness=133.211281,phi=0.000000,
bestindividual=	0.000000	3.000000	1.000000
      variance=	2.309401	2.309401	2.009375
   It means the current best solution is '0,3,1'

Note that for large models, a code template can be generated using the
'genSRESCcode' method in SPA/src/SPAmodel.java



