/***************************************************************
 *
 * Class Name:    BrownGenFF
 *
 * Description:   Generate FF c-code for Brown's model
 *               
 * Author: Liu Bing (lb@nus.edu.sg)
 *
 ***************************************************************/

public class BrownGenFF{
    public static void main(String[] args){
	SPAmodel m=new SPAmodel("brown");   
	m.loadVar("../models/brown/var.csv");
	m.loadPar("../models/brown/par.csv");
 	m.loadEqn("../models/brown/ode.txt");

	// unknown parameters
	int P[]={1,2,3,4,11,12,15,17,23,27,28,29,33,34,37,38,39,40,43,44};
	// data time points
	int T[]={2,5,10,20,30,40,50,60,80,100};
	// data variable
	int X[]={4,6,8,10,21,25,27};
	double weights[]={0.0336757,1,0.0455945,0.187954,0.00971087,0.184003,0.027279};

	// store all brownXBNT*.txt files in ""../models/brown/tables/"
	//m.genFFcode(100,100,"../models/brown/",0,P,T,X,weights,"../models/brown/brownNoisy.txt",true);
	m.genFFcodeC(100,100,"../models/brown/",0,P,T,X,weights,"../models/brown/brownNoisy.txt",true);

    }
}
