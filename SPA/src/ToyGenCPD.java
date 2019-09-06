/***************************************************************
 *
 * Class Name:    ToyGenCPD
 *
 * Description:   Compute conditional probabilities for the toy model
 *               
 * Author: Liu Bing (lb@nus.edu.sg)
 *
 ***************************************************************/

public class ToyGenCPD{
    public static void main(String[] args){
	SPAmodel m=new SPAmodel("toy");
	m.loadVar("../models/toy/var_example.csv");
	m.loadPar("../models/toy/par_example.csv");
	m.loadEqn("../models/toy/ode_example.txt");
	/*
	  processN(dir_num,       // the number of input directories (i.e. ct0, ct1 ..)
	           node_num,      // the max number of CPT copies within 1 input dir
		   timepoint_num, // numer of DBN time points
		   input_dir,     // the location of input directories (i.e. ct0, ct1 ..)
		   output_dir,    // output location
		   option)        // set 0 for current implementation
	*/
	m.processN(2,5,100,
		   "../models/toy/",
		   "../models/toy/ctn",0);
    }
}
