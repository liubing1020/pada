/***************************************************************
 *
 * Class Name:    ToyGenDBN
 *
 * Description:   Generate sampling C code for the toy model
 *               
 * Author: Liu Bing (lb@nus.edu.sg)
 *
 ***************************************************************/

public class ToyGenDBN{
    public static void main(String[] args){
	SPAmodel m=new SPAmodel("toy");
	m.loadVar("../models/toy/var.csv");
	m.loadPar("../models/toy/par.csv");
	m.loadEqn("../models/toy/ode.txt");
	
	String rnd="../3rdparty/dSFMT-src-2.0/dSFMT.c";
	String out_bat="../models/toy/batct";
	String out_mpi="../models/toy/mpict";

        /* C code for execution on single PC 
	  genBATcode(sample_size,      // number of sampling trajectories
	             time_point,       // number of DBN time points
                     delta_t,          // delta t of internal ODE solver (RK4)
                     total_time,       // duration of simulation
                     random_generator, // file location of random generator
                     output_dir)       // location of output counting files for 
                                       // computing probability tables (CPTs) 
	*/
	m.genBATcode(1000,100,0.01,10,rnd,out_bat);

        /* MPI code for execution on cluster 
	  genMPIcode(sample_size,      // number of sampling trajectories
	             time_point,       // number of DBN time points
                     delta_t,          // delta t of internal ODE solver (RK4)
                     total_time,       // duration of simulation
                     random_generator, // file location of random generator
                     output_dir)       // location of output files
	*/
	m.genMPIcode(1000,100,0.01,10,rnd,out_mpi);

    }
}
