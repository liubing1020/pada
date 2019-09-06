/***************************************************************
 *
 * Class Name:    ToyGenMPSA
 *
 * Description:   Generate MPSA code for toy model
 *               
 * Author: Liu Bing (lb@nus.edu.sg)
 *
 ***************************************************************/

public class ToyGenMPSA{
    public static void main(String[] args){
	SPAmodel m=new SPAmodel("toy");
	m.loadVar("../models/toy/var_example.csv");
	m.loadPar("../models/toy/par_example.csv");
	m.loadEqn("../models/toy/ode_example.txt");

	// unknown parameters
	int P[]={1,2,3};
	// data time points
	int T[]={10,20,30,40,50,60,70,80,90,100};
	// data variable
	int X[]={1,2,3,4};
	double weights[]={0.642113,1,0.776737,0.435402};

	// cvode version 
	m.genMPSAcode2(10000,"../../SPA/3rdparty/dSFMT-src-2.0/dSFMT.c",".",P,T,X,"../../SPA/models/toy/toydata.txt",200,3,0,0);

	//m.genMPSAcode(100,100,"../models/toy/ctn",0,P,T,X,weights,"../models/toy/toydata.txt",true);
	// FF version
	m.genMPSAcode(10,100,0.01,500,"../../SPA/3rdparty/dSFMT-src-2.0/dSFMT.c",".",P,T,X,weights,"../../SPA/models/toy/toydata.txt",3,1.0);
    }
}
