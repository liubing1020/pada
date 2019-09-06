/***************************************************************
 *
 * Class Name:    ToyGenInput
 *
 * Description:   Generate Input {ode.txt,var.csv,par,csv} for the toy model
 *               
 * Author: Liu Bing (lb@nus.edu.sg)
 *
 ***************************************************************/

public class ToyGenInput{
    public static void main(String[] args){
	SBML2SPA ss=new SBML2SPA("../../models/toy/toy.xml"); 
	ss.printODEs("../../models/toy/");
	ss.printPar("../../models/toy/");
	ss.printVar("../../models/toy/");
    }
}
