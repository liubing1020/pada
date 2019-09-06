/***************************************************************
 *
 * Class Name:    ToyGenMPI
 *
 * Description:   Generate sampling code for the toy model
 *               
 * Author: Liu Bing (lb@nus.edu.sg)
 *
 ***************************************************************/

public class ToyGenMPI{
    public static void main(String[] args){
	SPAmodel m=new SPAmodel("toy");
	m.loadVar("../models/toy/var.csv");
	m.loadPar("../models/toy/par.csv");
	m.loadEqn("../models/toy/ode.txt");

        /*
	  genMPIcode(sample_size,      // number of sampling trajectories
	             time_point,       // number of DBN time points
                     delta_t,          // delta t of internal ODE solver (RK4)
                     total_time,       // duration of simulation
                     random_generator, // file location of random generator
                     output_dir)       // location of output files
	*/
	m.genMPIcode(1000,100,0.01,10,
		     "/home/g06/g0601796/SPA/3rdparty/dSFMT-src-2.0/dSFMT.c",
		     "/home/g06/g0601796/SPA/models/toy/mpict");
    }
}
