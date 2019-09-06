/***************************************************************
 *
 * Class Name:    ToyGenFF
 *
 * Description:   Generate FF code for toy model
 *               
 * Author: Liu Bing (lb@nus.edu.sg)
 *
 ***************************************************************/

public class ToyGenFF{
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

	// java version 
	//m.genFFcode(100,100,"../models/toy/ctn",0,P,T,X,weights,"../models/toy/toydata.txt",true);
	// C version
	m.genFFcodeC(100,100,"../models/toy/ctn",0,P,T,X,weights,"../models/toy/toydata.txt",true);
    }
}
