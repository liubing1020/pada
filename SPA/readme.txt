============================================================================
==   A Brief Introdcution of Using Signaling Pathway Approximator (SPA)   ==
==                                                                        ==
==   Author: Bing Liu (liubing@comp.nus.edu.sg)                           ==
==   Organization: National University of Singapore                       ==
==                                                                        ==
============================================================================


- Step 1: Generate the DBN construction code 
  SPAmodel.java in ./src/ is used for generating various code for an ODE model.
  The input files of SPAmodel.java are {ode.txt, var.csv, par.csv} specifying
  the differential equations, upper/lower bounds of variables and parameters.
  We can either manually write these input files or generate them from a
  COPASI model (http://www.copasi.org). 
  For example, we fisrt build a model for a toy pathway using COPASI.  
    [result] toy.cps in ./models/toy/
  Then export SBML file from the COPASI model
    [result] toy.xml in ./models/toy/
  Convert SBML file to the input files of SPAmodel.java
    $cd ./src/sbml2spa
	$java ToyGenInput
    [result] {ode.txt, par.csv, var.csv} in ./models/toy/
  Note that var.csv and par.csv are needed to be modified mannully to config the
  discretization setting (i.e. the upper/lower bounds, interval number/sizes,etc)
  
  Once we have the input files, code for DBN construction can be generated
    $cd ..
  Modify ToyGenBDM.java to specify the locations of the in/output files
	$javac ToyGenDBN.java 
	$java ToyGenDBN
    [result] {toy_bat.c, toy_mpi.c} in current directory
  Note that toy_bat.c is used for the execution on single PC while toy_mpi.c is 
  used for the execution on cluster. Next we compile it depending on the processor.
  For example:
    $gcc -O3 -finline-functions -fomit-frame-pointer -DNDEBUG -fno-strict-aliasing --param max-inline-insns-single=1800 -std=c99 -msse2 -DHAVE_SSE2 -DDSFMT_MEXP=521 -o toy_bat toy_bat.c
	[result] toy_bat	
	$/opt/mpich/bin/mpicc -O3 -finline-functions -fomit-frame-pointer -DNDEBUG -fno-strict-aliasing --param max-inline-insns-single=1800 -std=c99 -msse2 -DHAVE_SSE2 -DDSFMT_MEXP=521 -o toy_mpi toy_mpi.c
	[result] toy_mpi
	
- Step 2: Execute the DBN construction code and generate the CPTs
  Executing toy_bat or toy_mpi will generate counting files stored in 
  SPA/models/toy/batct or SPA/models/toy/mpict, and CPTs can be computed from them.
  For example:
    $./toy_bat 0  // 0 is an arbitary seed for the random number generator 
    [result] toyCTx[1-4]T[0-99]_0.txt in SPA/models/toy/batct/
  Note toyCTx1T0_0.txt record the CPT of x1 at T0. The format of this file is
    [index of the CPT entry] [counting] 
  For example, according the the ODEs, x1 has parents {k1, k2, x1, x2, x3}, then
  the entry:
    603 10
  means
    k1 k2 x1 x2 x3 | x1' | count
    0  0  4  4  0    3     10
    i.e. among the sample trajectories hitting {k1=0,k2=0,k3=4,x1=4,x2=0) at time t0, 
	10 of them hit x1=3 at time t1.
  
  Similarly, we can also run toy_mpi to generate countings, e.g.
    $/opt/mpich/bin/mpirun -np 5 ./toy_mpi
  	[result] toyCTx[1-4]T[0-99]_[0-4].txt in BDM/models/toy/mpict/

  Move those counting files to the input directory (ct0, ct1 ..) of SPAmodel.java:
    $mv ../models/toy/batct/* ../models/toy/ct0/
	$mv ../models/toy/mpict/* ../models/toy/ct1/
  Then countings from different runs can be simply merged, and CPTs can be computed.    
	$javac ToyGenCPD.java
    $java ToyGenCPD
	[result] toyCTx[1-4]T[0-99].txt in SPA/models/toy/ctn/ \\ merged countings
    [result] toyPx[1-4]T[0-99].txt in SPA/models/toy/ctn/tables \\ CPTs
  Here similarly, an entry on CPT files as:	
    603 0.31
  means
    k1 k2 x1 x2 x3 | x1' | p
    0  0  4  4  0    3     0.31
    i.e. P(x1'=3|k1=0,k2=0,k3=4,x1=4,x2=0)=0.31
    where x1' is the variable x1 on the next time point.
  Since a DBN is stored as CPTs of a dynamic bayesian network, we now finish the
  DBN construction for the toy pathway.  
	
- Step 3: Using DBN to perform model analysis 
  We can do probabilitic inference on DBN. For instance, given the prior distribution
  of the initial concentrations and parameters in terms of intervals. We can use
  FF algorithm to infer the marginal or mean of each species at each time point.
  To generate FF code:
    $javac ToyGenFF.java
    $java ToyGenFF
    [result] toy_ff0.c
  toy_ff0.c can be used for model analysis such as parameter estimation, global
  sensitivity anlaysis. Here testFF.c is an example of using toy_ff0.c.	
	$gcc -std=c99 -o testFF testFF.c
	$./testFF
	[result] toy_T.csv and toy_M.csv
  toy_T.csv records the mean of each catspecies at each time point, it can be used to
  plot the time profile of each species.
  toy_M.csv records the marginal probabilities of each species at each time points,
  format: 1,5,3,0.025 means Prob(X5 at T1 in Bin 3)=0.025)
  it can further be used for model checking.  
