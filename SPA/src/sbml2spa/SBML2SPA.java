/**
 * @file    SBML2SPA.java
 * @brief   Convert SBML Document to SPA input files
 * @author  Liu Bing (translated from libSBML Java examples)
 *
 */


import org.sbml.libsbml.*;
import java.io.*;

public class SBML2SPA{
    Model model;

    public SBML2SPA(String filename){
	SBMLReader reader = new SBMLReader();
	SBMLDocument document;
	int level, version;
	document = reader.readSBML(filename);
	if (document.getNumErrors() > 0){
	    document.printErrors();
	    System.out.println("Please correct the above problems first.");
	    System.exit(1);
	}
	this.model = document.getModel();
    }

    String getRateEqns (String species){
	String result="",formula,sid;
	boolean flag;
	Reaction r;
	for (int i = 0; i < model.getNumReactions(); i++){
	    r=model.getReaction(i);
	    formula="";
	    flag=false;
	    // scan reactants
	    for(int j=0;j<r.getNumReactants();j++){
		sid=r.getReactant(j).getSpecies();
		if(sid.equals(species)){
		    formula="-"+getReactionMath(r);
		    flag=true;break;
		}
	    }
	    
	    // scan products
	    for(int j=0;j<r.getNumProducts();j++){
		sid=r.getProduct(j).getSpecies();
		if(sid.equals(species)){
		    if(flag){
			formula="";flag=false;break;
		    } else{
			formula="+"+getReactionMath(r);flag=true;break;
		    }
		}
	    }
	    if(flag) result+=formula;
	}
	return result;
    }

    String getReactionMath (Reaction r){
	String formula="";
	String fun="function_1";
	if (r.isSetKineticLaw()){
	    KineticLaw kl = r.getKineticLaw();
	    //formula = kl.getFormula();
	    
	    if (kl.isSetMath()){
		ASTNode math = kl.getMath();
		if(!fun.equals(math.getRightChild().getName())){ // mass law
		    formula=libsbml.formulaToString(math);
		} else { // pre-defined function
		    math=math.getRightChild();
		    String k=math.getChild(0).getName();
		    String E=math.getChild(1).getName();
		    String S=math.getChild(2).getName();
		    String Km=math.getChild(3).getName();
		    formula=k+"*"+E+"*"+S+"/("+Km+"+"+S+")";
		}
	    }
		//formula = libsbml.formulaToString(kl.getMath());
		//}
	}
	return formula;
    }
    
    
    void printODEs(String outdir){
	try{
	    FileOutputStream outfile=new FileOutputStream(outdir+"./ode.txt");
	    PrintWriter out=new PrintWriter(outfile);
	    //out.println(model.getNumSpecies());
	    for (int i=0;i<model.getNumSpecies();i++){
		out.println(convert(model.getSpecies(i).getId()+
					   "="+getRateEqns(model.getSpecies(i).getId())));
	    }
	    out.flush();
	} catch(IOException e){
	    e.printStackTrace();
	}
    }
    
    void printPar(String outdir){
	try{
	    FileOutputStream outfile=new FileOutputStream(outdir+"./par.csv");
	    PrintWriter out=new PrintWriter(outfile);
	    // header
	    out.println("NAME,INIT,LowerBound,UpperBound,BoundNum,BoundSize");
	    for (int i=0;i<model.getNumParameters();i++){
		Parameter p=model.getParameter(i);
		//out.println(convert(p.getId())+","+p.getValue()+",,,");
		double x=p.getValue();
		if(x<1)
		    out.println(convert(p.getId())+","+p.getValue()+",0,1,5");
		else if(x<100)
		    out.println(convert(p.getId())+","+p.getValue()+",0,100,5");
		else
		    out.println(convert(p.getId())+","+p.getValue()+",0,10000,5");
	    }
	    
	    out.flush();
	} catch(IOException e){
	    e.printStackTrace();
	}
    }

    void printVar(String outdir){
	try{
	    FileOutputStream outfile=new FileOutputStream(outdir+"./var.csv");
	    PrintWriter out=new PrintWriter(outfile);
	    // header
	    out.println("NAME,INIT,LowerBound,UpperBound,BoundNum,BoundSize");
	    for (int i=0;i<model.getNumSpecies();i++){
		Species p=model.getSpecies(i);
		out.println(convert(p.getId())+","
			    +p.getInitialConcentration()+",0,15,5");
	    }
	    
	    out.flush();
	} catch(IOException e){
	    e.printStackTrace();
	}
    }

    
    
    public static void main (String[] args){

	if (args.length != 1){
	    System.out.println("Usage: java SBML2SPA filename");
	    System.exit(1);
	}
	
	String filename = args[0];
	SBML2SPA ss=new SBML2SPA(filename);
	//ss.printODEs();
	//ss.printPar();
	//ss.printVar();
	/*
	for (int n = 0; n < ss.model.getNumFunctionDefinitions(); ++n){
	    ss.printFunctionDefinition(n + 1, ss.model.getFunctionDefinition(n));
	    System.out.println("");
	}
	//for (int n = 0; n < ss.model.getNumReactions(); ++n)
	//  System.out.println(n+ " "+ss.getReactionMath(ss.model.getReaction(n)));
	int n=3;
	System.out.println(n+ " "+ss.getReactionMath(ss.model.getReaction(n)));
	*/
    }
    
    
    void printFunctionDefinition (int n, FunctionDefinition fd){
	if (fd.isSetMath()){
	    System.out.print("FunctionDefinition " + n + ", " + fd.getId() + "(");
	    
	    ASTNode math = fd.getMath();
	    
	    /* Print function arguments. */
	    if (math.getNumChildren() > 1){
		System.out.print(" " + math.getLeftChild().getName());
                
		for (int i = 1; i < math.getNumChildren() - 1; ++i){
		    System.out.print(",  " + math.getChild(i).getName());
		}
	    }

	    System.out.print(") := ");
	    
	    /* Print function body. */
	    if (math.getNumChildren() == 0){
		System.out.println("(no body defined)");
	    } else {
		math = math.getChild(math.getNumChildren() - 1);
		System.out.println(libsbml.formulaToString(math));
	    }
	}
    }
    

    void printRuleMath (int n, Rule r){
	if (r.isSetMath()){
	    String formula = libsbml.formulaToString(r.getMath());
	    String var     = r.getVariable();
	    
	    if (r.getVariable().length() > 0)
		System.out.println("Rule " + n + ", formula: " + var + " = " + formula);
	    else
		System.out.println("Rule " + n + ", formula: " + formula + " = 0");
	}
    }
    
    
    void printReactionMath (int n, Reaction r){
	if (r.isSetKineticLaw()){
	    KineticLaw kl = r.getKineticLaw();
            
	    if ( kl.isSetMath() ){
		
		//ASTNode ast = kl.getMath();
		//long formula = ast.getNumChildren();
		String formula = libsbml.formulaToString( kl.getMath() );
		System.out.println("Reaction " + n + ", formula: " + convert(formula));
	    }
	}
    }
    
    
    void printEventAssignmentMath (int n, EventAssignment ea){
	if (ea.isSetMath()){
	    String variable = ea.getVariable();
	    String formula  = libsbml.formulaToString( ea.getMath() );
	    
	    System.out.println("  EventAssignment " + n + ", trigger: " + variable + " = " +
		    formula);
	}
    }
    

    void printEventMath (int n, Event e){
	String formula;
	
	
	if (e.isSetDelay()){
	    formula = libsbml.formulaToString(e.getDelay().getMath());
	    System.out.println("Event " + n + " delay: " + formula);
	}

	if (e.isSetTrigger()){
	    formula = libsbml.formulaToString(e.getTrigger().getMath());
	    System.out.println("Event " + n + " trigger: " + formula);
	}
	
	for (int i = 0; i < e.getNumEventAssignments(); ++i){
	    printEventAssignmentMath(i + 1, e.getEventAssignment(i));
	}
	
	System.out.println("\n");
    }
    
    
    void printMath (){
	for (int n = 0; n < model.getNumFunctionDefinitions(); ++n){
	    printFunctionDefinition(n + 1, model.getFunctionDefinition(n));
	    System.out.println("");
	}
	
	for (int n = 0; n < model.getNumRules(); ++n){
	    printRuleMath(n + 1, model.getRule(n));
	    System.out.println("");
	}

	for (int n = 0; n < model.getNumReactions(); ++n){
	    printReactionMath(n + 1, model.getReaction(n));
	    System.out.println("");
	}

	for (int n = 0; n < model.getNumEvents(); ++n){
	    printEventMath(n + 1, model.getEvent(n));
	}

	


    }
    
    String convert (String formula){
	// step1: remove compartment_1
	formula=formula.replace("compartment_1 * ","");
	formula=formula.replace("species_","x");
	formula=formula.replace("parameter_","k");
	formula=formula.replace(" ","");
	formula=formula.replace("=+","=");
	// step2: replace function_1
// 	String result="";
// 	StringTokenizer st=new StringTokenizer(line, " =()+-*/");
// 	int tNum=st.countTokens();
// 	for(int j=0;j<tNum;j++){
// 	    String element=st.nextToken();
// 	    if(element.charAt(0)==''	      
// 	       result+=model.getSpecies();
	       
	return formula;
	       
    }
  /**
   * Loads the SWIG-generated libSBML Java module when this class is
   * loaded, or reports a sensible diagnostic message about why it failed.
   */
  static
  {
    String varname;

    if (System.getProperty("mrj.version") != null)
      varname = "DYLD_LIBRARY_PATH";	// We're on a Mac.
    else
      varname = "LD_LIBRARY_PATH";	// We're not on a Mac.

    try
    {
      System.loadLibrary("sbmlj");
      // For extra safety, check that the jar file is in the classpath.
      Class.forName("org.sbml.libsbml.libsbml");
    }
    catch (UnsatisfiedLinkError e)
    {
      System.err.println("Error: could not link with the libSBML library."+
			 "  It is likely\nyour " + varname +
			 " environment variable does not include\nthe"+
			 " directory containing the libsbml library file.");
      System.exit(1);
    }
    catch (ClassNotFoundException e)
    {
      System.err.println("Error: unable to load the file libsbmlj.jar."+
			 "  It is likely\nyour " + varname + " environment"+
			 " variable or CLASSPATH variable\ndoes not include"+
			 " the directory containing the libsbmlj.jar file.");
      System.exit(1);
    }
    catch (SecurityException e)
    {
      System.err.println("Could not load the libSBML library files due to a"+
			 " security exception.");
    }
  }
}

