import java.io.*;
import java.util.*;

public class SPAmodel{
    String modelName;
    int varNum; // N+1
    int parNum; // M+1
    Variable [] vars;
    Parameter [] pars;
	
    boolean hasCPT;

    // option 1: compact representation
    Hashtable cpt;

    public SPAmodel(String name){
	modelName=name;
	hasCPT=false;
	cpt=new Hashtable();
    }

    
    //------------- loading input data ----------------------

    public void loadVar(String filename){
	try{
	    int MAX=1000;
	    String []lines=new String[MAX];
	    FileInputStream fis =new FileInputStream(filename);
	    BufferedReader in=new BufferedReader(new InputStreamReader(fis));
	    in.readLine(); // skip table header
	    int ctr=1;
	    
	    while((lines[ctr++]=in.readLine())!=null){} // read data
	    
	    // Variable initialization
	    varNum=ctr-1;
	    vars=new Variable[varNum];
	    for(int i=1;i<varNum;i++) 
		vars[i]=new Variable(lines[i],i);
	} catch(IOException e){
	    e.printStackTrace();
	}
    
    }

    public void loadPar(String filename){
	try{
	    int MAX=1000;
	    String []lines=new String[MAX];
	    FileInputStream fis =new FileInputStream(filename);
	    BufferedReader in=new BufferedReader(new InputStreamReader(fis));
	    in.readLine(); // skip table header
	    int ctr=1;
	    
	    while((lines[ctr++]=in.readLine())!=null){} // read data
	    
	    // Variable initialization
	    parNum=ctr-1;
	    pars=new Parameter[parNum];
	    for(int i=1;i<parNum;i++) 
		pars[i]=new Parameter(lines[i],i);
	} catch(IOException e){
	    e.printStackTrace();
	}
    }

    public void loadEqn(String filename){
	try{
	    if(varNum==0||parNum==0){
		System.out.println("load variable and parameter first.");
		return;
	    }
	    
	    updateEqn(filename);

	    FileInputStream fis =new FileInputStream(filename);
	    BufferedReader in=new BufferedReader(new InputStreamReader(fis));
	    String line;
	    
	    for(int i=1;i<varNum;i++){
		int []vMap=new int[varNum];
		int []pMap=new int[parNum];

		line=in.readLine(); System.out.println(line);
		StringTokenizer st=new StringTokenizer(line, "=()+-*/ \t");
		int tNum=st.countTokens();
		for(int j=0;j<tNum;j++){
		    String element=st.nextToken();
		    boolean find=false;

		    for(int v=1;v<varNum;v++){
			// only add dynamic variables
			if((vars[v].name).equals(element)&&
			   !((vars[v].eqn).equals(vars[v].name+"=0"))){
			    vMap[v]=1;
			    find=true;break;
			}
		    }
		    if(!find){
			for(int p=1;p<parNum;p++){
			    // only add unknown parameters
			    if((pars[p].name).equals(element)&&(pars[p].binNum>1)){
				pMap[p]=1;break;
			    }
			}
		    }
		    
		}
		
		if(vars[i].init<0) vMap[i]=0; else vMap[i]=1;

		// [patch]: for the wrong CPT of brown's model
		if(modelName.equals("brown"))
		    if(i==20||i==21) vMap[31]=1;

		int ctr,ctr2;
		
		ctr=0;
		for(int v=1;v<varNum;v++) if(vMap[v]==1) ctr++;
		int varIDs[]=new int[ctr];
		ctr2=0;
		for(int v=1;v<varNum;v++) if(vMap[v]==1) varIDs[ctr2++]=v;
		
		ctr=0;
		for(int p=1;p<parNum;p++) if(pMap[p]==1) ctr++;
		int parIDs[]=new int[ctr];
		ctr2=0;
		for(int p=1;p<parNum;p++) if(pMap[p]==1) parIDs[ctr2++]=p;
		
		vars[i].setEquation(line,varIDs,parIDs);
	    }
	} catch(IOException e){
	    e.printStackTrace();
	}
    }

    public void updateEqn(String filename){
	try{
	    if(varNum==0||parNum==0){
		System.out.println("load variable and parameter first.");
		return;
	    }
	    
	    FileInputStream fis =new FileInputStream(filename);
	    BufferedReader in=new BufferedReader(new InputStreamReader(fis));
	    String line;
	    
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0){ // only update real var 
		    line=in.readLine();
		    vars[i].updateEquation(line);
		} else if(vars[i].eqn==null){
		    vars[i].updateEquation("empty");
		}
	    }
	} catch(IOException e){
	    e.printStackTrace();
	}
    }

    
    
    //------------- generating MPI code -------------------
    // sampling code for cluster
    
    
    public void genMPIcode(int sampleSize,int timeptNum,double dt,double simTime,String dSFMTPath,String outputDir){
	try{
	    FileOutputStream outfile=new FileOutputStream(modelName+"_mpi.c");
	    PrintWriter out=new PrintWriter(outfile);

	    // header
	    out.println("#include \"mpi.h\"");
	    out.println("#include <stdio.h>");
	    out.println("#include <math.h>");
	    out.println("#include \""+dSFMTPath+"\"");
	    out.println("");

	    // global variable
	    for(int i=1;i<varNum;i++) 
		out.println("double "+vars[i].name+","+vars[i].name+"p,"+vars[i].name+"preN;");
	    for(int i=1;i<parNum;i++)
		out.println("double "+pars[i].name+";");
	    for(int i=1;i<varNum;i++){
		int n=(vars[i].varIDs).length+(vars[i].parIDs).length+1;
		out.print("int "+vars[i].name+"ctr["+timeptNum+"]");
		//for(int j=0;j<n;j++) out.print("["+vars[i].binNum+"]"); // bug july13
		int tmp_v=(vars[i].varIDs).length;
		int tmp_p=(vars[i].parIDs).length;
		for(int j=0;j<tmp_p;j++)
		    out.print("["+pars[(vars[i].parIDs)[j]].binNum+"]");
		for(int j=0;j<tmp_v;j++)
		    out.print("["+vars[(vars[i].varIDs)[j]].binNum+"]");
		out.print("["+vars[i].binNum+"]");
		out.println(";");
	    }
	    out.println("");

	    // functions
	    for(int i=1;i<varNum;i++){
		StringTokenizer st = new StringTokenizer(vars[i].eqn, "=");
		String name=st.nextToken(); //name
		String body=st.nextToken(); // body
		out.println("double f"+name+"(double "+name+"){");
		out.println("return "+body+";");
		out.println("}");	
	    }
	    out.println("");

	    // discretization function
	    out.println("");
	    out.println("int discretize(double v, double xi[], int length){");
	    out.println("for(int j=1;j<length-1;j++) if(v<xi[j]) return j-1;");
	    out.println("return length-2;");
	    out.println("}");
	    out.println("");

	    // main
	    out.println("int main(int argc, char *argv[]){");
	    out.println("double dt="+dt+",halfdt=dt/2.0;");
	    out.println("int tps=(int)("+simTime+"/dt),t,i;");
	    out.println("int block=tps/"+timeptNum+",tb;");
	    out.println("double halfF1,halfF2,F3,F4;");

	    // --------- variable 
	    for(int i=1;i<varNum;i++){
		out.print("double "+vars[i].name+"i[]={");
		double tmp_r=vars[i].upperBound-vars[i].lowerBound;
		double tmp_b=vars[i].lowerBound;
		for(int j=0;j<vars[i].binNum;j++){
		    out.print(tmp_b+",");
		    tmp_b+=vars[i].binSizes[j]*tmp_r;
		}
		out.println(vars[i].upperBound+"};");
	    }
	    out.println("");  
	    for(int i=1;i<varNum;i++) 
		out.println("int "+vars[i].name+"pre,"+vars[i].name+
			    "post,"+vars[i].name+"init="+vars[i].binInit+";");
	    // --------- parameter
	    for(int i=1;i<parNum;i++){
		out.print("double "+pars[i].name+"i[]={");
		double tmp_r=pars[i].upperBound-pars[i].lowerBound;
		double tmp_b=pars[i].lowerBound;
		for(int j=0;j<pars[i].binNum;j++){
		    out.print(tmp_b+",");
		    tmp_b+=pars[i].binSizes[j]*tmp_r;
		}
		out.println(pars[i].upperBound+"};");
	    }
	    for(int i=1;i<parNum;i++) 
		out.println("int "+pars[i].name+"bin;");
	    
	    Random seed=new Random();
	    out.println("");
	    out.println("");
	    out.println("// MPI ----------------");
	    out.println("double startwtime = 0.0, endwtime;");
	    out.println("int  namelen,numprocs,myid;");
	    out.println("char processor_name[MPI_MAX_PROCESSOR_NAME];");
	    out.println("MPI_Init(&argc,&argv);");
	    out.println("MPI_Comm_size(MPI_COMM_WORLD,&numprocs);");
	    out.println("MPI_Comm_rank(MPI_COMM_WORLD,&myid);");
	    out.println("MPI_Get_processor_name(processor_name,&namelen);");
	    out.println("fprintf(stderr,\"Process %d on %s\\n\",myid, processor_name);");
	    out.println("if(myid==0){");
	    out.println("startwtime = MPI_Wtime();");
	    out.println("}");
	    out.println("// -------------------");
	    out.println("");
	    out.println("");
	    out.println("int sampleNo="+sampleSize+";");
	    out.println("dsfmt_t dsfmt;");
	    out.println("int seed="+seed.nextInt(12345)+"+myid;");
	    out.println("dsfmt_init_gen_rand(&dsfmt,seed);");
	    out.println("");
	    out.println("");
	    out.println("for(int i=0;i<sampleNo;i++){");
	    for(int i=1;i<parNum;i++){
		out.println(pars[i].name+"="+pars[i].lowerBound+
			    "+dsfmt_genrand_close_open(&dsfmt)*"+
			    (pars[i].upperBound-pars[i].lowerBound)+";");
	    }
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0)
		    out.println(vars[i].name+"="+vars[i].name+"i["+vars[i].name+
				"init]+dsfmt_genrand_close_open(&dsfmt)*("+
				vars[i].name+"i["+vars[i].name+"init+1]-"+
				vars[i].name+"i["+vars[i].name+"init]);");
		else
		    out.println(vars[i].name+"=f"+vars[i].name+"(0);");
	    }
	    for(int i=1;i<varNum;i++){
		//if(vars[i].init>=0)
		out.println(vars[i].name+"preN="+vars[i].name+";");
		    //else
		    //out.println(vars[i].name+"preN=f"+vars[i].name+"(0);");
	    }
	    for(int i=1;i<parNum;i++) 
		out.println(pars[i].name+"bin=discretize("+pars[i].name+","+pars[i].name+"i,"+(pars[i].binNum+1)+");");
	    
	    out.println("");

	    // ODE solver
	    out.println("for(int t=1;t<=tps;t++){");
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0){
		    out.println("// "+vars[i].name);
		    out.println("halfF1=halfdt*f"+vars[i].name+"("+vars[i].name+");");
		    out.println("halfF2=halfdt*f"+vars[i].name+"("+vars[i].name+"+halfF1);");
		    out.println("F3=dt*f"+vars[i].name+"("+vars[i].name+"+halfF2);");
		    out.println("F4=dt*f"+vars[i].name+"("+vars[i].name+"+F3);");
		    out.println(vars[i].name+"p="+vars[i].name+"+(2*halfF1+4*halfF2+2*F3+F4)/6.0;");
		}
	    }
	    
	    out.println("if(t%block==0){");
	    out.println("tb=t/block-1;");
	    
	    for(int i=1;i<varNum;i++) 
		out.println(vars[i].name+"pre=discretize("+
			    vars[i].name+"preN,"+vars[i].name+"i,"+(vars[i].binNum+1)+");");
	    out.println("");

	    for(int i=1;i<varNum;i++) 
		if(vars[i].init>=0)
		    out.println(vars[i].name+"post=discretize("+
				vars[i].name+","+vars[i].name+"i,"+(vars[i].binNum+1)+");");
	    out.println("");

	    out.println("");
	    for(int i=1;i<varNum;i++){
		out.print(vars[i].name+"ctr[tb]");
		for(int j=0;j<(vars[i].parIDs).length;j++) 
		    out.print("["+pars[vars[i].parIDs[j]].name+"bin]");
		for(int j=0;j<(vars[i].varIDs).length;j++) 
		    out.print("["+vars[vars[i].varIDs[j]].name+"pre]");
		if(vars[i].init>=0)
		    out.println("["+vars[i].name+"post]++;");
		else
		    out.println("["+vars[i].name+"pre]++;");
	    }
	    for(int i=1;i<varNum;i++)
		//if(vars[i].init>=0)
		out.println(vars[i].name+"preN="+vars[i].name+";");
	    //else
	    //    out.println(vars[i].name+"preN=f"+vars[i].name+"(0);");
	    out.println("}");
	    
	    for(int i=1;i<varNum;i++) 
		if(vars[i].init>=0)
		    out.println(vars[i].name+"="+vars[i].name+"p;");
		else
		    out.println(vars[i].name+"=f"+vars[i].name+"(0);");
	    
	    out.println("}");
	    out.println("}");
	    
	    out.println("");
	    
	    // output
	    out.println("");
	    out.println("// output");
	    out.println("FILE *out;");
	    out.println("char buffer[256];");
	    out.println("snprintf(buffer, sizeof(buffer), \"dummy.txt\");");
	    out.println("int idx=0;");
	    out.println("");
	    out.println("for(tb=0;tb<"+timeptNum+";tb++){");
	    for(int i=1;i<varNum;i++){
		out.println("snprintf(buffer, sizeof(buffer), \""+
			    outputDir+"/"+modelName+"CTx"+i+
			    "T%d_%d.txt\", tb,myid);");
		out.println("out=fopen(buffer, \"w\");");
		out.println("idx=0;");
		out.println("");
		out.println("");
		int tmp_p=(vars[i].parIDs).length;
		int tmp_v=+(vars[i].varIDs).length;
		for(int j=0;j<tmp_p;j++) 
		    out.println("for(int ki"+j+"=0;ki"+j+"<"
				+pars[vars[i].parIDs[j]].binNum+";ki"+j+"++)");
		for(int j=0;j<tmp_v;j++) 
		    out.println("for(int vi"+j+"=0;vi"+j+"<"
				+vars[vars[i].varIDs[j]].binNum+";vi"+j+"++)");
		out.println("for(int vi=0;vi<"+vars[i].binNum+";vi++)");
		out.println("{");
		
		out.print("int ctrtmp=("+vars[i].name+"ctr[tb]");
		for(int j=0;j<tmp_p;j++) out.print("[ki"+j+"]");
		for(int j=0;j<tmp_v;j++) out.print("[vi"+j+"]");
		out.println("[vi]);");
		out.println("if(ctrtmp>0){");
		out.println("fprintf(out,\"%d %d\\n\",idx,ctrtmp);");
		out.println("}");
		out.println("idx++;");
		out.println("}");
		out.println("fclose(out);");	
	    }

	    out.println("}");
	    out.println("// MPI ---------------------------------");
	    out.println("");

	    out.println("if (myid == 0){");
	    out.println("endwtime = MPI_Wtime();");
	    out.println("printf(\"time = %f\\n\",endwtime-startwtime);");
	    out.println("}");
	    out.println("");
	    out.println("MPI_Finalize();");
	    out.println("//--------------------------------------");  


	    out.println("");
	    out.println("return 0;");

	    out.println("}");

	    out.flush();
	    
	    
	} catch(IOException e){
	    e.printStackTrace();
	}
    }



    //------------- generating MPI code -------------------
    // new sampling method
    
    public void genMPIcode3(int sampleSize,int timeptNum,double dt,double simTime,String dSFMTPath,String outputDir){
	try{
	    FileOutputStream outfile=new FileOutputStream(modelName+"_mpi.c");
	    PrintWriter out=new PrintWriter(outfile);

	    // header
	    out.println("#include \"mpi.h\"");
	    out.println("#include <stdio.h>");
	    out.println("#include <math.h>");
	    out.println("#include \""+dSFMTPath+"\"");
	    out.println("");

	    // initilization
	    out.println("");
	    out.println("");
	    out.println("");
	    out.println("int sampleNo="+sampleSize+";");
	    out.println("dsfmt_t dsfmt;");
	    out.println("");
	    out.println("");

	    // global variable
	    out.println("int tps,t,i,block,tb;");
	    out.println("double dt,halfdt,halfF1,halfF2,F3,F4;");

	    // --------- variable 
	    for(int i=1;i<varNum;i++){
		out.print("double "+vars[i].name+"i[]={");
		double tmp_r=vars[i].upperBound-vars[i].lowerBound;
		double tmp_b=vars[i].lowerBound;
		for(int j=0;j<vars[i].binNum;j++){
		    out.print(tmp_b+",");
		    tmp_b+=vars[i].binSizes[j]*tmp_r;
		}
		out.println(vars[i].upperBound+"};");
	    }
	    out.println("");  
	    for(int i=1;i<varNum;i++) 
		out.println("int "+vars[i].name+"pre,"+vars[i].name+
			    "post,"+vars[i].name+"init="+vars[i].binInit+";");
	    // --------- parameter
	    for(int i=1;i<parNum;i++){
		out.print("double "+pars[i].name+"i[]={");
		double tmp_r=pars[i].upperBound-pars[i].lowerBound;
		double tmp_b=pars[i].lowerBound;
		for(int j=0;j<pars[i].binNum;j++){
		    out.print(tmp_b+",");
		    tmp_b+=pars[i].binSizes[j]*tmp_r;
		}
		out.println(pars[i].upperBound+"};");
	    }
	    for(int i=1;i<parNum;i++) 
		out.println("int "+pars[i].name+"bin;");


	    for(int i=1;i<varNum;i++) 
		out.println("double "+vars[i].name+","+vars[i].name+"p,"+vars[i].name+"preN;");
	    for(int i=1;i<parNum;i++)
		out.println("double "+pars[i].name+";");
	    for(int i=1;i<varNum;i++){
		int n=(vars[i].varIDs).length+(vars[i].parIDs).length+1;
		out.print("int "+vars[i].name+"ctr["+timeptNum+"]");
		int tmp_v=(vars[i].varIDs).length;
		int tmp_p=(vars[i].parIDs).length;
		for(int j=0;j<tmp_p;j++)
		    out.print("["+pars[(vars[i].parIDs)[j]].binNum+"]");
		for(int j=0;j<tmp_v;j++)
		    out.print("["+vars[(vars[i].varIDs)[j]].binNum+"]");
		out.print("["+vars[i].binNum+"]");
		out.println(";");
	    }
	    out.println("");


	    // functions
	    for(int i=1;i<varNum;i++){
		StringTokenizer st = new StringTokenizer(vars[i].eqn, "=");
		String name=st.nextToken(); //name
		String body=st.nextToken(); // body
		out.println("double f"+name+"(double "+name+"){");
		out.println("return "+body+";");
		out.println("}");	
	    }
	    out.println("");

	    // discretization function
	    out.println("");
	    out.println("int discretize(double v, double xi[], int length){");
	    out.println("for(int j=1;j<length-1;j++) if(v<xi[j]) return j-1;");
	    out.println("return length-2;");
	    out.println("}");
	    out.println("");


	    // sampling function
	    out.println("");
	    out.println("void sampling(){");
	    for(int i=1;i<parNum;i++){
		out.println(pars[i].name+"="+pars[i].lowerBound+
			    "+dsfmt_genrand_close_open(&dsfmt)*"+
			    (pars[i].upperBound-pars[i].lowerBound)+";");
	    }
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0)
		    out.println(vars[i].name+"="+vars[i].name+"i["+vars[i].name+
				"init]+dsfmt_genrand_close_open(&dsfmt)*("+
				vars[i].name+"i["+vars[i].name+"init+1]-"+
				vars[i].name+"i["+vars[i].name+"init]);");
		else
		    out.println(vars[i].name+"=f"+vars[i].name+"(0);");
	    }
	    out.println("}");
	    out.println("");


	    // solve function
	    out.println("");
	    out.println("void solve(){");
	    for(int i=1;i<varNum;i++){
		//if(vars[i].init>=0)
		out.println(vars[i].name+"preN="+vars[i].name+";");
		    //else
		    //out.println(vars[i].name+"preN=f"+vars[i].name+"(0);");
	    }
	    for(int i=1;i<parNum;i++) 
		out.println(pars[i].name+"bin=discretize("+pars[i].name+","+pars[i].name+"i,"+(pars[i].binNum+1)+");");
	    
	    out.println("");

	    // ODE solver
	    out.println("for(int t=1;t<=tps;t++){");
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0){
		    out.println("// "+vars[i].name);
		    out.println("halfF1=halfdt*f"+vars[i].name+"("+vars[i].name+");");
		    out.println("halfF2=halfdt*f"+vars[i].name+"("+vars[i].name+"+halfF1);");
		    out.println("F3=dt*f"+vars[i].name+"("+vars[i].name+"+halfF2);");
		    out.println("F4=dt*f"+vars[i].name+"("+vars[i].name+"+F3);");
		    out.println(vars[i].name+"p="+vars[i].name+"+(2*halfF1+4*halfF2+2*F3+F4)/6.0;");
		}
	    }
	    
	    out.println("if(t%block==0){");
	    out.println("tb=t/block-1;");
	    
	    for(int i=1;i<varNum;i++) 
		out.println(vars[i].name+"pre=discretize("+
			    vars[i].name+"preN,"+vars[i].name+"i,"+(vars[i].binNum+1)+");");
	    out.println("");

	    for(int i=1;i<varNum;i++) 
		if(vars[i].init>=0)
		    out.println(vars[i].name+"post=discretize("+
				vars[i].name+","+vars[i].name+"i,"+(vars[i].binNum+1)+");");
	    out.println("");

	    out.println("");
	    for(int i=1;i<varNum;i++){
		out.print(vars[i].name+"ctr[tb]");
		for(int j=0;j<(vars[i].parIDs).length;j++) 
		    out.print("["+pars[vars[i].parIDs[j]].name+"bin]");
		for(int j=0;j<(vars[i].varIDs).length;j++) 
		    out.print("["+vars[vars[i].varIDs[j]].name+"pre]");
		if(vars[i].init>=0)
		    out.println("["+vars[i].name+"post]++;");
		else
		    out.println("["+vars[i].name+"pre]++;");
	    }
	    for(int i=1;i<varNum;i++)
		//if(vars[i].init>=0)
		out.println(vars[i].name+"preN="+vars[i].name+";");
	    //else
	    //    out.println(vars[i].name+"preN=f"+vars[i].name+"(0);");
	    out.println("}");
	    
	    for(int i=1;i<varNum;i++) 
		if(vars[i].init>=0)
		    out.println(vars[i].name+"="+vars[i].name+"p;");
		else
		    out.println(vars[i].name+"=f"+vars[i].name+"(0);");
	    
	    out.println("}");
	    out.println("}");
	    //out.println("}");
	    out.println("");


	    // main
	    out.println("int main(int argc, char *argv[]){");
	    //out.println("int myid=atoi(argv[1]);");
	    out.println("dt="+dt+",halfdt=dt/2.0;");
	    out.println("tps=(int)("+simTime+"/dt);");
	    out.println("block=tps/"+timeptNum+";");
	    
 	    Random seed=new Random();
	    out.println("");
	    out.println("");
	    out.println("");
	    out.println("");
	    out.println("// MPI ----------------");
	    out.println("double startwtime = 0.0, endwtime;");
	    out.println("int  namelen,numprocs,myid;");
	    out.println("char processor_name[MPI_MAX_PROCESSOR_NAME];");
	    out.println("MPI_Init(&argc,&argv);");
	    out.println("MPI_Comm_size(MPI_COMM_WORLD,&numprocs);");
	    out.println("MPI_Comm_rank(MPI_COMM_WORLD,&myid);");
	    out.println("MPI_Get_processor_name(processor_name,&namelen);");
	    out.println("fprintf(stderr,\"Process %d on %s\\n\",myid, processor_name);");
	    out.println("if(myid==0){");
	    out.println("startwtime = MPI_Wtime();");
	    out.println("}");
	    out.println("// -------------------");
	    out.println("");
 	    out.println("int seed="+seed.nextInt(12345)+"+myid;");
	    out.println("dsfmt_init_gen_rand(&dsfmt,seed);");
	    
	    // sampling
 	    out.println("for(int i=0;i<sampleNo;i++){");
	    for(int i=1;i<varNum;i++){
		out.println("// "+vars[i].name+" samples");
		int tmp_p=(vars[i].parIDs).length;
		for(int j=0;j<tmp_p;j++) 
		    out.println("for(int ki"+j+"=0;ki"+j+"<"
				+pars[vars[i].parIDs[j]].binNum+";ki"+j+"++)");
		out.println("{");
		out.println("sampling();");
		out.println(""); 
		for(int j=0;j<tmp_p;j++)
		    out.println(pars[vars[i].parIDs[j]].name+"="+pars[vars[i].parIDs[j]].name+"i[ki"+j+
				"]+dsfmt_genrand_close_open(&dsfmt)*("+
				pars[vars[i].parIDs[j]].name+"i[ki"+j+"+1]-"+
				pars[vars[i].parIDs[j]].name+"i[ki"+j+"]);");
		
		out.println(""); 
		out.println("solve();");
		out.println("}");
		out.println("");
	    }
	    out.println("}");
	    
	    out.println("");
	    
	    // output
	    out.println("");
	    out.println("// output");
	    out.println("FILE *out;");
	    out.println("char buffer[256];");
	    out.println("snprintf(buffer, sizeof(buffer), \"dummy.txt\");");
	    out.println("int idx=0;");
	    out.println("");
	    out.println("for(tb=0;tb<"+timeptNum+";tb++){");
	    for(int i=1;i<varNum;i++){
		out.println("snprintf(buffer, sizeof(buffer), \""+
			    outputDir+"/"+modelName+"CTx"+i+
			    "T%d_%d.txt\", tb,myid);");
		out.println("out=fopen(buffer, \"w\");");
		out.println("idx=0;");
		out.println("");
		out.println("");
		int tmp_p=(vars[i].parIDs).length;
		int tmp_v=+(vars[i].varIDs).length;
		for(int j=0;j<tmp_p;j++) 
		    out.println("for(int ki"+j+"=0;ki"+j+"<"
				+pars[vars[i].parIDs[j]].binNum+";ki"+j+"++)");
		for(int j=0;j<tmp_v;j++) 
		    out.println("for(int vi"+j+"=0;vi"+j+"<"
				+vars[vars[i].varIDs[j]].binNum+";vi"+j+"++)");
		out.println("for(int vi=0;vi<"+vars[i].binNum+";vi++)");
		out.println("{");
		
		out.print("int ctrtmp=("+vars[i].name+"ctr[tb]");
		for(int j=0;j<tmp_p;j++) out.print("[ki"+j+"]");
		for(int j=0;j<tmp_v;j++) out.print("[vi"+j+"]");
		out.println("[vi]);");
		out.println("if(ctrtmp>0){");
		out.println("fprintf(out,\"%d %d\\n\",idx,ctrtmp);");
		out.println("}");
		out.println("idx++;");
		out.println("}");
		out.println("fclose(out);");	
	    }

	    out.println("}");
	    out.println("// MPI ---------------------------------");
	    out.println("");

	    out.println("if (myid == 0){");
	    out.println("endwtime = MPI_Wtime();");
	    out.println("printf(\"time = %f\\n\",endwtime-startwtime);");
	    out.println("}");
	    out.println("");
	    out.println("MPI_Finalize();");
	    out.println("//--------------------------------------");  

	    out.println("");
	    out.println("return 0;");

	    out.println("}");

	    out.flush();
	    
    
	} catch(IOException e){
	    e.printStackTrace();
	}
    }
    





    //------------- generating Batch code -------------------
    // sampling code for single PC
    
    
    public void genBATcode(int sampleSize,int timeptNum,double dt,double simTime,String dSFMTPath,String outputDir){
	try{
	    FileOutputStream outfile=new FileOutputStream(modelName+"_bat.c");
	    PrintWriter out=new PrintWriter(outfile);

	    // header
	    out.println("#include <stdio.h>");
	    out.println("#include <math.h>");
	    out.println("#include \""+dSFMTPath+"\"");
	    out.println("");

	    // global variable
	    for(int i=1;i<varNum;i++) 
		out.println("double "+vars[i].name+","+vars[i].name+"p,"+vars[i].name+"preN;");
	    for(int i=1;i<parNum;i++)
		out.println("double "+pars[i].name+";");
	    for(int i=1;i<varNum;i++){
		int n=(vars[i].varIDs).length+(vars[i].parIDs).length+1;
		out.print("int "+vars[i].name+"ctr["+timeptNum+"]");
		int tmp_v=(vars[i].varIDs).length;
		int tmp_p=(vars[i].parIDs).length;
		for(int j=0;j<tmp_p;j++)
		    out.print("["+pars[(vars[i].parIDs)[j]].binNum+"]");
		for(int j=0;j<tmp_v;j++)
		    out.print("["+vars[(vars[i].varIDs)[j]].binNum+"]");
		out.print("["+vars[i].binNum+"]");
		out.println(";");
	    }
	    out.println("");

	    // functions
	    for(int i=1;i<varNum;i++){
		StringTokenizer st = new StringTokenizer(vars[i].eqn, "=");
		String name=st.nextToken(); //name
		String body=st.nextToken(); // body
		out.println("double f"+name+"(double "+name+"){");
		out.println("return "+body+";");
		out.println("}");	
	    }
	    out.println("");

	    // discretization function
	    out.println("");
	    out.println("int discretize(double v, double xi[], int length){");
	    out.println("for(int j=1;j<length-1;j++) if(v<xi[j]) return j-1;");
	    out.println("return length-2;");
	    out.println("}");
	    out.println("");

	    // main
	    out.println("int main(int argc, char *argv[]){");
	    out.println("int myid=atoi(argv[1]);");
	    out.println("double dt="+dt+",halfdt=dt/2.0;");
	    out.println("int tps=(int)("+simTime+"/dt),t,i;");
	    out.println("int block=tps/"+timeptNum+",tb;");
	    out.println("double halfF1,halfF2,F3,F4;");

	    // --------- variable 
	    for(int i=1;i<varNum;i++){
		out.print("double "+vars[i].name+"i[]={");
		double tmp_r=vars[i].upperBound-vars[i].lowerBound;
		double tmp_b=vars[i].lowerBound;
		for(int j=0;j<vars[i].binNum;j++){
		    out.print(tmp_b+",");
		    tmp_b+=vars[i].binSizes[j]*tmp_r;
		}
		out.println(vars[i].upperBound+"};");
	    }
	    out.println("");  
	    for(int i=1;i<varNum;i++) 
		out.println("int "+vars[i].name+"pre,"+vars[i].name+
			    "post,"+vars[i].name+"init="+vars[i].binInit+";");
	    // --------- parameter
	    for(int i=1;i<parNum;i++){
		out.print("double "+pars[i].name+"i[]={");
		double tmp_r=pars[i].upperBound-pars[i].lowerBound;
		double tmp_b=pars[i].lowerBound;
		for(int j=0;j<pars[i].binNum;j++){
		    out.print(tmp_b+",");
		    tmp_b+=pars[i].binSizes[j]*tmp_r;
		}
		out.println(pars[i].upperBound+"};");
	    }
	    for(int i=1;i<parNum;i++) 
		out.println("int "+pars[i].name+"bin;");
	    
	    Random seed=new Random();
	    out.println("");
	    out.println("");
	    out.println("");
	    out.println("int sampleNo="+sampleSize+";");
	    out.println("dsfmt_t dsfmt;");
	    out.println("int seed="+seed.nextInt(12345)+"+myid;");
	    out.println("dsfmt_init_gen_rand(&dsfmt,seed);");
	    out.println("");
	    out.println("");
	    out.println("for(int i=0;i<sampleNo;i++){");
	    for(int i=1;i<parNum;i++){
		out.println(pars[i].name+"="+pars[i].lowerBound+
			    "+dsfmt_genrand_close_open(&dsfmt)*"+
			    (pars[i].upperBound-pars[i].lowerBound)+";");
	    }
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0)
		    out.println(vars[i].name+"="+vars[i].name+"i["+vars[i].name+
				"init]+dsfmt_genrand_close_open(&dsfmt)*("+
				vars[i].name+"i["+vars[i].name+"init+1]-"+
				vars[i].name+"i["+vars[i].name+"init]);");
		else
		    out.println(vars[i].name+"=f"+vars[i].name+"(0);");
	    }
	    for(int i=1;i<varNum;i++){
		//if(vars[i].init>=0)
		out.println(vars[i].name+"preN="+vars[i].name+";");
		    //else
		    //out.println(vars[i].name+"preN=f"+vars[i].name+"(0);");
	    }
	    for(int i=1;i<parNum;i++) 
		out.println(pars[i].name+"bin=discretize("+pars[i].name+","+pars[i].name+"i,"+(pars[i].binNum+1)+");");
	    
	    out.println("");

	    // ODE solver
	    out.println("for(int t=1;t<=tps;t++){");
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0){
		    out.println("// "+vars[i].name);
		    out.println("halfF1=halfdt*f"+vars[i].name+"("+vars[i].name+");");
		    out.println("halfF2=halfdt*f"+vars[i].name+"("+vars[i].name+"+halfF1);");
		    out.println("F3=dt*f"+vars[i].name+"("+vars[i].name+"+halfF2);");
		    out.println("F4=dt*f"+vars[i].name+"("+vars[i].name+"+F3);");
		    out.println(vars[i].name+"p="+vars[i].name+"+(2*halfF1+4*halfF2+2*F3+F4)/6.0;");
		}
	    }
	    
	    out.println("if(t%block==0){");
	    out.println("tb=t/block-1;");
	    
	    for(int i=1;i<varNum;i++) 
		out.println(vars[i].name+"pre=discretize("+
			    vars[i].name+"preN,"+vars[i].name+"i,"+(vars[i].binNum+1)+");");
	    out.println("");

	    for(int i=1;i<varNum;i++) 
		if(vars[i].init>=0)
		    out.println(vars[i].name+"post=discretize("+
				vars[i].name+","+vars[i].name+"i,"+(vars[i].binNum+1)+");");
	    out.println("");

	    out.println("");
	    for(int i=1;i<varNum;i++){
		out.print(vars[i].name+"ctr[tb]");
		for(int j=0;j<(vars[i].parIDs).length;j++) 
		    out.print("["+pars[vars[i].parIDs[j]].name+"bin]");
		for(int j=0;j<(vars[i].varIDs).length;j++) 
		    out.print("["+vars[vars[i].varIDs[j]].name+"pre]");
		if(vars[i].init>=0)
		    out.println("["+vars[i].name+"post]++;");
		else
		    out.println("["+vars[i].name+"pre]++;");
	    }
	    for(int i=1;i<varNum;i++)
		//if(vars[i].init>=0)
		out.println(vars[i].name+"preN="+vars[i].name+";");
	    //else
	    //    out.println(vars[i].name+"preN=f"+vars[i].name+"(0);");
	    out.println("}");
	    
	    for(int i=1;i<varNum;i++) 
		if(vars[i].init>=0)
		    out.println(vars[i].name+"="+vars[i].name+"p;");
		else
		    out.println(vars[i].name+"=f"+vars[i].name+"(0);");
	    
	    out.println("}");
	    out.println("}");
	    
	    out.println("");
	    
	    // output
	    out.println("");
	    out.println("// output");
	    out.println("FILE *out;");
	    out.println("char buffer[256];");
	    out.println("snprintf(buffer, sizeof(buffer), \"dummy.txt\");");
	    out.println("int idx=0;");
	    out.println("");
	    out.println("for(tb=0;tb<"+timeptNum+";tb++){");
	    for(int i=1;i<varNum;i++){
		out.println("snprintf(buffer, sizeof(buffer), \""+
			    outputDir+"/"+modelName+"CTx"+i+
			    "T%d_%d.txt\", tb,myid);");
		out.println("out=fopen(buffer, \"w\");");
		out.println("idx=0;");
		out.println("");
		out.println("");
		int tmp_p=(vars[i].parIDs).length;
		int tmp_v=+(vars[i].varIDs).length;
		for(int j=0;j<tmp_p;j++) 
		    out.println("for(int ki"+j+"=0;ki"+j+"<"
				+pars[vars[i].parIDs[j]].binNum+";ki"+j+"++)");
		for(int j=0;j<tmp_v;j++) 
		    out.println("for(int vi"+j+"=0;vi"+j+"<"
				+vars[vars[i].varIDs[j]].binNum+";vi"+j+"++)");
		out.println("for(int vi=0;vi<"+vars[i].binNum+";vi++)");
		out.println("{");
		
		out.print("int ctrtmp=("+vars[i].name+"ctr[tb]");
		for(int j=0;j<tmp_p;j++) out.print("[ki"+j+"]");
		for(int j=0;j<tmp_v;j++) out.print("[vi"+j+"]");
		out.println("[vi]);");
		out.println("if(ctrtmp>0){");
		out.println("fprintf(out,\"%d %d\\n\",idx,ctrtmp);");
		out.println("}");
		out.println("idx++;");
		out.println("}");
		out.println("fclose(out);");	
	    }

	    out.println("}");
	    out.println("");
	    out.println("return 0;");

	    out.println("}");

	    out.flush();
	    
	    
	} catch(IOException e){
	    e.printStackTrace();
	}
    }


    //------------- generating Batch code using Cvode -------------------
    // todo: intermedia nodes
    
    public void genBATcvode(int sampleSize,int timeptNum,double dt,double simTime,String dSFMTPath,String outputDir){
	try{
	    FileOutputStream outfile=new FileOutputStream(modelName+"bat.c");
	    PrintWriter out=new PrintWriter(outfile);

	    // header
	    out.println("#include <stdio.h>");
	    out.println("#include <math.h>");
	    out.println("#include \""+dSFMTPath+"\"");
	    out.println(" ");
	    //suchee 11Nov
	    out.println("#include <cvode/cvode.h>");
	    out.println("#include <nvector/nvector_serial.h>");
	    out.println("#include <cvode/cvode_dense.h>");
	    out.println("#include <sundials/sundials_dense.h>");
	    out.println("#include <sundials/sundials_types.h>");
	    out.println(" ");
	    out.println("#define Ith(v,i)    NV_Ith_S(v,i-1)");
	    out.println("#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)");
	    out.println(" ");
	    out.println("#define NEQ   "+(varNum-1));
        out.println("#define RTOL  1.0e-4");
            
        int ppp=0;
        for(int i=1;i<varNum;i++)
        {
        if(vars[i].init>=0)
        {
			ppp+=1;
			out.println("#define ATOL"+ppp+"  1.0e-26");
        }        
        
        }
        
        out.println("#define T0    0.0");
        out.println("#define T1    1.0");
        out.println("#define TMULT 1.0");
        out.println("#define NOUT  "+timeptNum);
        out.println(" ");
        out.println(" ");
        out.println("static int f(realtype t, N_Vector x,  N_Vector xdot, void *user_data);");
        out.println(" ");
        out.println("static int check_flag(void *flagvalue, char *funcname, int opt);") ; 
        out.println(" "); 
        
        // End add Suchee 11 Nov    
    
        // global variable
	    for(int i=1;i<varNum;i++) 
		out.println("double p"+vars[i].name+","+vars[i].name+","+vars[i].name+"p,"+vars[i].name+"preN;");
	    for(int i=1;i<parNum;i++)
		out.println("double "+pars[i].name+";");
	    for(int i=1;i<varNum;i++){
		int n=(vars[i].varIDs).length+(vars[i].parIDs).length+1;
		out.print("int "+vars[i].name+"ctr["+timeptNum+"]");
		//for(int j=0;j<n;j++) out.print("["+vars[i].binNum+"]"); // bug july13
		int tmp_v=(vars[i].varIDs).length;
		int tmp_p=(vars[i].parIDs).length;
		for(int j=0;j<tmp_p;j++)
		    out.print("["+pars[(vars[i].parIDs)[j]].binNum+"]");
		for(int j=0;j<tmp_v;j++)
		    out.print("["+vars[(vars[i].varIDs)[j]].binNum+"]");
		out.print("["+vars[i].binNum+"]");
		out.println(";");
	    }
	    out.println("");

	    // functions
	    for(int i=1;i<varNum;i++){
	    
	    if(!(vars[i].init>=0))
	    {
		StringTokenizer st = new StringTokenizer(vars[i].eqn, "=");
		String name=st.nextToken(); //name
		String body=st.nextToken(); // body
		out.println("double f"+name+"(double "+name+"){");
		out.println("return "+body+";");
		out.println("}");
		}	
	    }
	    out.println("");

	    // discretization function
	    out.println("");
	    out.println("int discretize(double v, double xi[], int length){");
	    out.println("for(int j=1;j<length-1;j++) if(v<xi[j]) return j-1;");
	    out.println("return length-2;");
	    out.println("}");
	    out.println("");

	    // main
	    out.println("int main(int argc, char *argv[]){");
	    out.println(" ");
	    out.println("int myid=atoi(argv[1]);");
	    out.println("realtype reltol, t, tout;");//Suchee 11 Nov
		out.println("N_Vector x, abstol;");//Suchee 11 Nov
		out.println("void *cvode_mem;");//Suchee 11 Nov
		out.println("int flag,iout;");//Suchee 11 Nov
		out.println("x = abstol = NULL;");//Suchee 11 Nov
		out.println("cvode_mem = NULL;");//Suchee 11 Nov
	    
	    //out.println("double dt="+dt+",halfdt=dt/2.0;"); Suchee 11 Nov
	    //out.println("int tps=(int)("+simTime+"/dt),t,i;"); Suchee 11 Nov
	    //out.println("int block=tps/"+timeptNum+",tb;"); Suchee 11 Nov
	    //out.println("double halfF1,halfF2,F3,F4;"); Suchee 11 Nov

	    // --------- variable 
	    for(int i=1;i<varNum;i++){
		out.print("double "+vars[i].name+"i[]={");
		double tmp_r=vars[i].upperBound-vars[i].lowerBound;
		double tmp_b=vars[i].lowerBound;
		for(int j=0;j<vars[i].binNum;j++){
		    out.print(tmp_b+",");
		    tmp_b+=vars[i].binSizes[j]*tmp_r;
		}
		out.println(vars[i].upperBound+"};");
	    }
	    out.println("");  
	    for(int i=1;i<varNum;i++) 
		out.println("int "+vars[i].name+"pre,"+vars[i].name+
			    "post,"+vars[i].name+"init="+vars[i].binInit+";");
	    // --------- parameter
	    for(int i=1;i<parNum;i++){
		out.print("double "+pars[i].name+"i[]={");
		double tmp_r=pars[i].upperBound-pars[i].lowerBound;
		double tmp_b=pars[i].lowerBound;
		for(int j=0;j<pars[i].binNum;j++){
		    out.print(tmp_b+",");
		    tmp_b+=pars[i].binSizes[j]*tmp_r;
		}
		out.println(pars[i].upperBound+"};");
	    }
	    for(int i=1;i<parNum;i++) 
		out.println("int "+pars[i].name+"bin;");
	    
	    Random seed=new Random();
	    out.println("");
	    out.println("");
	    out.println("");
	    out.println("int sampleNo="+sampleSize+";");
	    out.println("dsfmt_t dsfmt;");
	    out.println("int seed="+seed.nextInt(12345)+"+myid;");
	    out.println("dsfmt_init_gen_rand(&dsfmt,seed);");
	    out.println("");
	    out.println("");
	    out.println("for(int i=0;i<sampleNo;i++){");
	    out.println("");	    
	    out.println("x = N_VNew_Serial(NEQ);");//Suchee 11-11
		out.println("if (check_flag((void *)x, \"N_VNew_Serial\", 0)) return(1);");// Suchee 11 Nov
		out.println("abstol = N_VNew_Serial(NEQ);");// Suchee 11 Nov
		out.println("if (check_flag((void *)abstol, \"N_VNew_Serial\", 0)) return(1);");// Suchee 11 Nov 
		out.println("");
		out.println("");
		
	    
	    
	    for(int i=1;i<parNum;i++){
		out.println(pars[i].name+"="+pars[i].lowerBound+
			    "+dsfmt_genrand_close_open(&dsfmt)*"+
			    (pars[i].upperBound-pars[i].lowerBound)+";");
	    }
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0)
		{
		    out.println("p"+vars[i].name+"="+vars[i].name+"i["+vars[i].name+
				"init]+dsfmt_genrand_close_open(&dsfmt)*("+
				vars[i].name+"i["+vars[i].name+"init+1]-"+
				vars[i].name+"i["+vars[i].name+"init]);");
			out.println("Ith(x,"+i+")=p"+vars[i].name+";");
			//out.println(vars[i].name+"=p"+vars[i].name+";");			
			out.println("Ith(abstol,"+i+") = ATOL"+i+";");
			}
		//else
		    //out.println(vars[i].name+"=f"+vars[i].name+"(0);");
	    }
	    
	  	    
	    out.println("");
	    out.println("");
	    out.println("");
	    out.println("reltol = RTOL;");
	    
    
	    out.println("cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);");
		out.println("if (check_flag((void *)cvode_mem, \"CVodeCreate\", 0)) return(1);");
	    out.println("");
	    out.println("");
	    out.println("");
	    
	    out.println("flag = CVodeInit(cvode_mem, f, T0, x);");
		out.println("if (check_flag(&flag, \"CVodeInit\", 1)) return(1);");
		
		out.println("");
		out.println("");
		out.println("");


		out.println("flag = CVodeSVtolerances(cvode_mem, reltol, abstol);");
		out.println("if (check_flag(&flag, \"CVodeSVtolerances\", 1)) return(1);");
		
		out.println("");
		out.println("");
		out.println("");
		

		out.println("flag = CVDense(cvode_mem, NEQ);");
		out.println("if (check_flag(&flag, \"CVDense\", 1)) return(1);");
		
		out.println("");
		out.println("");
		out.println("");    
	    
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0)
		{
		out.println(vars[i].name+"preN=p"+vars[i].name+";");
		out.println(vars[i].name+"=p"+vars[i].name+";");
		}
		
		else
		out.println(vars[i].name+"preN=f"+vars[i].name+"(0);");
	    }
	    for(int i=1;i<parNum;i++) 
		out.println(pars[i].name+"bin=discretize("+pars[i].name+","+pars[i].name+"i,"+(pars[i].binNum+1)+");");
	    
	    out.println("");
	    
	    out.println("iout = 0;  tout = T1;"); // Suchee 11 Nov
	    

	    // ODE solver
	    out.println("for(int tt=0;tt<NOUT;tt++){");
	    out.println("");
	    out.println("");
	    
	     for(int i=1;i<varNum;i++) 
		out.println(vars[i].name+"pre=discretize("+
			    vars[i].name+"preN,"+vars[i].name+"i,"+(vars[i].binNum+1)+");");
	    out.println("");
	    
	    
	    
	    
	    
	    //for(int i=1;i<varNum;i++){
					
		//if(vars[i].init>=0){
		    //out.println("// "+vars[i].name);
		    //out.println("halfF1=halfdt*f"+vars[i].name+"("+vars[i].name+");");
		    //out.println("halfF2=halfdt*f"+vars[i].name+"("+vars[i].name+"+halfF1);");
		    //out.println("F3=dt*f"+vars[i].name+"("+vars[i].name+"+halfF2);");
		    //out.println("F4=dt*f"+vars[i].name+"("+vars[i].name+"+F3);");
		    //out.println(vars[i].name+"p="+vars[i].name+"+(2*halfF1+4*halfF2+2*F3+F4)/6.0;");
		    
		    
		    
		//}
		//}
		
		out.println("flag = CVode(cvode_mem, tout, x, &t, CV_NORMAL);");
		out.println("");
	  
	    
	    //out.println("if(t%block==0){");
	    //out.println("tb=t/block-1;");
	    
	    //for(int i=1;i<varNum;i++) 
		//out.println(vars[i].name+"pre=discretize("+
			    //vars[i].name+"preN,"+vars[i].name+"i,"+(vars[i].binNum+1)+");");
	    //out.println("");

	    //for(int i=1;i<varNum;i++) 
		//if(vars[i].init>=0)
		    //out.println(vars[i].name+"post=discretize("+
				//vars[i].name+","+vars[i].name+"i,"+(vars[i].binNum+1)+");");
	    //out.println("");
	    
	    // suchee 11-Nov
	    int p=0;
	    int pp=0;
	    for(int i=1;i<varNum;i++)
	    {
	    if(vars[i].init>=0)
			{
			pp+=1;
			out.println(vars[i].name+"=Ith(x,"+pp+");");
			}
	    }
	    
	    for(int i=1;i<varNum;i++) 
	    {
			if(vars[i].init>=0)
			 out.println(vars[i].name+"post=discretize("+vars[i].name+","+vars[i].name+"i,"+(vars[i].binNum+1)+");");
			 else
			 out.println(vars[i].name+"post=discretize(f"+vars[i].name+"(0),"+vars[i].name+"i,"+(vars[i].binNum+1)+");");
			    
		  }
		 out.println("");
	    
	     out.println("if (check_flag(&flag, \"CVode\", 1)) break;");

		 out.println("if (flag == CV_SUCCESS) {");
		 out.println("tout += TMULT;");
		 out.println("	}");
	        

	    out.println("");
	    for(int i=1;i<varNum;i++){
		out.print(vars[i].name+"ctr[tt]");
		for(int j=0;j<(vars[i].parIDs).length;j++) 
		    out.print("["+pars[vars[i].parIDs[j]].name+"bin]");
		for(int j=0;j<(vars[i].varIDs).length;j++) 
		    out.print("["+vars[vars[i].varIDs[j]].name+"pre]");
		if(vars[i].init>=0)
		    out.println("["+vars[i].name+"post]++;");
		else
		    out.println("["+vars[i].name+"pre]++;");
	    }
	    for(int i=1;i<varNum;i++)
	    {
		if(vars[i].init>=0)
		//out.println(vars[i].name+"preN=Ith(x,"+i+");");
		out.println(vars[i].name+"preN="+vars[i].name+";");
	    else
	    out.println(vars[i].name+"preN=f"+vars[i].name+"(0);");
	    }
	    out.println("}");
	    
	    //for(int i=1;i<varNum;i++) 
		//if(vars[i].init>=0)
		    //out.println(vars[i].name+"="+vars[i].name+"p;");
		//else
		    //out.println(vars[i].name+"=f"+vars[i].name+"(0);");
	    
	    out.println("}");
	  	    
	    out.println("");
	    
	    // output
	    out.println("");
	    out.println("// output");
	    out.println("FILE *out;");
	    out.println("char buffer[256];");
	    out.println("snprintf(buffer, sizeof(buffer), \"dummy.txt\");");
	    out.println("int idx=0;");
	    out.println("int tb=0;");
	    
	    out.println("");
	    out.println("for(tb=0;tb<"+timeptNum+";tb++){");
	    for(int i=1;i<varNum;i++){
		out.println("snprintf(buffer, sizeof(buffer), \""+
			    outputDir+"/"+modelName+"CTx"+i+
			    "T%d_%d.txt\", tb,myid);");
		out.println("out=fopen(buffer, \"w\");");
		out.println("idx=0;");
		out.println("");
		out.println("");
		int tmp_p=(vars[i].parIDs).length;
		int tmp_v=+(vars[i].varIDs).length;
		for(int j=0;j<tmp_p;j++) 
		    out.println("for(int ki"+j+"=0;ki"+j+"<"
				+pars[vars[i].parIDs[j]].binNum+";ki"+j+"++)");
		for(int j=0;j<tmp_v;j++) 
		    out.println("for(int vi"+j+"=0;vi"+j+"<"
				+vars[vars[i].varIDs[j]].binNum+";vi"+j+"++)");
		out.println("for(int vi=0;vi<"+vars[i].binNum+";vi++)");
		out.println("{");
		
		out.print("int ctrtmp=("+vars[i].name+"ctr[tb]");
		for(int j=0;j<tmp_p;j++) out.print("[ki"+j+"]");
		for(int j=0;j<tmp_v;j++) out.print("[vi"+j+"]");
		out.println("[vi]);");
		out.println("if(ctrtmp>0){");
		out.println("fprintf(out,\"%d %d\\n\",idx,ctrtmp);");
		out.println("}");
		out.println("idx++;");
		out.println("}");
		out.println("fclose(out);");	
	    }

	    out.println("}");
	    out.println("");
	    out.println("return 0;");

	    out.println("}");
	    
	    out.println("static int check_flag(void *flagvalue, char *funcname, int opt)");
	    out.println("	{");
	    out.println("		int *errflag;");
	    out.println("");
	    out.println("		/* Check if SUNDIALS function returned NULL pointer - no memory allocated */");
	    out.println("		if (opt == 0 && flagvalue == NULL) {");
	    out.println("			fprintf(stderr, \"\\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\\n\\n\",");
	    out.println("				funcname);");
	    out.println("			return(1); }");
	    out.println("");
	    out.println("		/* Check if flag < 0 */");
	    out.println("		else if (opt == 1) {");
	    out.println("			errflag = (int *) flagvalue;");
	    out.println("			if (*errflag < 0) {");
	    out.println("				fprintf(stderr, \"\\nSUNDIALS_ERROR: %s() failed with flag = %d\\n\\n\",");
	    out.println("					funcname, *errflag);");
	    out.println("				return(1); }}");
	    out.println("");
	    out.println("		/* Check if function returned NULL pointer - no memory allocated */");
	    out.println("		else if (opt == 2 && flagvalue == NULL) {");
	    out.println("			fprintf(stderr,\"\\nMEMORY_ERROR: %s() failed - returned NULL pointer\\n\\n\",");
	    out.println("				funcname);");
	    out.println("			return(1); }");
	    out.println("");
	    out.println("		return(0);");
	    out.println("	}");
	    
	    out.println("");
	    out.println("");
		
	    
	    out.println("static int f(realtype t, N_Vector x, N_Vector xdot, void *user_data)");
	    out.println("	{");
	    for(int i=1;i<varNum;i++) 
		{
		    if(vars[i].init>=0)
			{
			    out.println("realtype "+vars[i].name+";");
			}
		}
	    
	    int pppp=0;
	    for(int i=1;i<varNum;i++){
		
		StringTokenizer st = new StringTokenizer(vars[i].eqn, "=");
		String name=st.nextToken(); //name
		String body=st.nextToken(); // body
		if(vars[i].init>=0)
		    {
			pppp+=1;
			out.println(name+"= Ith(x,"+pppp+");");
		    }
		
	    }
	    for(int i=1;i<varNum;i++){
		StringTokenizer st = new StringTokenizer(vars[i].eqn, "=");
		String name=st.nextToken(); //name
		String body=st.nextToken(); // body
		int ppppp=0;
		if(vars[i].init>=0)
		    {
			ppppp+=1;
			out.println("		Ith(xdot,"+i+") = "+body+";");
		    }
		
		
	    }
	    out.println("		return(0);");
	    out.println("	}");
	    
	    
	    
	    out.flush();
	    
	    
	} catch(IOException e){
	    e.printStackTrace();
	}
    }

    //------------- generating Batch code -------------------
    // fix par bin
    
    public void genBATcode2(int sampleSize,int timeptNum,double dt,double simTime,String dSFMTPath,String outputDir){
	try{
	    FileOutputStream outfile=new FileOutputStream(modelName+"bat2.c");
	    PrintWriter out=new PrintWriter(outfile);

	    // header
	    out.println("#include <stdio.h>");
	    out.println("#include <math.h>");
	    out.println("#include \""+dSFMTPath+"\"");
	    out.println("");

	    // global variable
	    for(int i=1;i<varNum;i++) 
		out.println("double "+vars[i].name+","+vars[i].name+"p,"+vars[i].name+"preN;");
	    for(int i=1;i<parNum;i++)
		out.println("double "+pars[i].name+";");
	    for(int i=1;i<varNum;i++){
		int n=(vars[i].varIDs).length+(vars[i].parIDs).length+1;
		out.print("int "+vars[i].name+"ctr["+timeptNum+"]");
		//for(int j=0;j<n;j++) out.print("["+vars[i].binNum+"]"); //bug july13
		int tmp_v=(vars[i].varIDs).length;
		int tmp_p=(vars[i].parIDs).length;
		for(int j=0;j<tmp_p;j++)
		    out.print("["+pars[(vars[i].parIDs)[j]].binNum+"]");
		for(int j=0;j<tmp_v;j++)
		    out.print("["+vars[(vars[i].varIDs)[j]].binNum+"]");
		out.print("["+vars[i].binNum+"]");
		out.println(";");
	    }
	    out.println("");

	    // functions
	    for(int i=1;i<varNum;i++){
		StringTokenizer st = new StringTokenizer(vars[i].eqn, "=");
		String name=st.nextToken(); //name
		String body=st.nextToken(); // body
		out.println("double f"+name+"(double "+name+"){");
		out.println("return "+body+";");
		out.println("}");	
	    }
	    out.println("");

	    // discretization function
	    out.println("");
	    out.println("int discretize(double v, double xi[], int length){");
	    out.println("for(int j=1;j<length-1;j++) if(v<xi[j]) return j-1;");
	    out.println("return length-2;");
	    out.println("}");
	    out.println("");

	    // main
	    out.println("int main(int argc, char *argv[]){");
	    out.println("int myid=atoi(argv[1]);");
	    out.println("double dt="+dt+",halfdt=dt/2.0;");
	    out.println("int tps=(int)("+simTime+"/dt),t,i;");
	    out.println("int block=tps/"+timeptNum+",tb;");
	    out.println("double halfF1,halfF2,F3,F4;");

	    // --------- variable 
	    for(int i=1;i<varNum;i++){
		out.print("double "+vars[i].name+"i[]={");
		double tmp_r=vars[i].upperBound-vars[i].lowerBound;
		double tmp_b=vars[i].lowerBound;
		for(int j=0;j<vars[i].binNum;j++){
		    out.print(tmp_b+",");
		    tmp_b+=vars[i].binSizes[j]*tmp_r;
		}
		out.println(vars[i].upperBound+"};");
	    }
	    out.println("");  
	    for(int i=1;i<varNum;i++) 
		out.println("int "+vars[i].name+"pre,"+vars[i].name+
			    "post,"+vars[i].name+"init="+vars[i].binInit+";");
	    // --------- parameter
	    for(int i=1;i<parNum;i++){
		out.print("double "+pars[i].name+"i[]={");
		double tmp_r=pars[i].upperBound-pars[i].lowerBound;
		double tmp_b=pars[i].lowerBound;
		for(int j=0;j<pars[i].binNum;j++){
		    out.print(tmp_b+",");
		    tmp_b+=pars[i].binSizes[j]*tmp_r;
		}
		out.println(pars[i].upperBound+"};");
	    }
	    for(int i=1;i<parNum;i++) 
		out.println("int "+pars[i].name+"bin;");
	    for(int i=1;i<parNum;i++) 
		out.println("int "+pars[i].name+"init="+pars[i].binInit+";");
	    
	    Random seed=new Random();
	    out.println("");
	    out.println("");
	    out.println("");
	    out.println("int sampleNo="+sampleSize+";");
	    out.println("dsfmt_t dsfmt;");
	    out.println("int seed="+seed.nextInt(12345)+"+myid;");
	    out.println("dsfmt_init_gen_rand(&dsfmt,seed);");
	    out.println("");
	    out.println("");
	    out.println("for(int i=0;i<sampleNo;i++){");
	    for(int i=1;i<parNum;i++){
		out.println(pars[i].name+"="+pars[i].name+"i["+pars[i].name+
			    "init]+dsfmt_genrand_close_open(&dsfmt)*("+
			    pars[i].name+"i["+pars[i].name+"init+1]-"+
			    pars[i].name+"i["+pars[i].name+"init]);");
	    }
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0)
		    out.println(vars[i].name+"="+vars[i].name+"i["+vars[i].name+
				"init]+dsfmt_genrand_close_open(&dsfmt)*("+
				vars[i].name+"i["+vars[i].name+"init+1]-"+
				vars[i].name+"i["+vars[i].name+"init]);");
		else
		    out.println(vars[i].name+"=f"+vars[i].name+"(0);");
	    }
	    for(int i=1;i<varNum;i++){
		//if(vars[i].init>=0)
		out.println(vars[i].name+"preN="+vars[i].name+";");
		    //else
		    //out.println(vars[i].name+"preN=f"+vars[i].name+"(0);");
	    }
	    for(int i=1;i<parNum;i++) 
		out.println(pars[i].name+"bin=discretize("+pars[i].name+","+pars[i].name+"i,"+(pars[i].binNum+1)+");");
	    
	    out.println("");

	    // ODE solver
	    out.println("for(int t=1;t<=tps;t++){");
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0){
		    out.println("// "+vars[i].name);
		    out.println("halfF1=halfdt*f"+vars[i].name+"("+vars[i].name+");");
		    out.println("halfF2=halfdt*f"+vars[i].name+"("+vars[i].name+"+halfF1);");
		    out.println("F3=dt*f"+vars[i].name+"("+vars[i].name+"+halfF2);");
		    out.println("F4=dt*f"+vars[i].name+"("+vars[i].name+"+F3);");
		    out.println(vars[i].name+"p="+vars[i].name+"+(2*halfF1+4*halfF2+2*F3+F4)/6.0;");
		}
	    }
	    
	    out.println("if(t%block==0){");
	    out.println("tb=t/block-1;");
	    
	    for(int i=1;i<varNum;i++) 
		out.println(vars[i].name+"pre=discretize("+
			    vars[i].name+"preN,"+vars[i].name+"i,"+(vars[i].binNum+1)+");");
	    out.println("");

	    for(int i=1;i<varNum;i++) 
		if(vars[i].init>=0)
		    out.println(vars[i].name+"post=discretize("+
				vars[i].name+","+vars[i].name+"i,"+(vars[i].binNum+1)+");");
	    out.println("");

	    out.println("");
	    for(int i=1;i<varNum;i++){
		out.print(vars[i].name+"ctr[tb]");
		for(int j=0;j<(vars[i].parIDs).length;j++) 
		    out.print("["+pars[vars[i].parIDs[j]].name+"bin]");
		for(int j=0;j<(vars[i].varIDs).length;j++) 
		    out.print("["+vars[vars[i].varIDs[j]].name+"pre]");
		if(vars[i].init>=0)
		    out.println("["+vars[i].name+"post]++;");
		else
		    out.println("["+vars[i].name+"pre]++;");
	    }
	    for(int i=1;i<varNum;i++)
		//if(vars[i].init>=0)
		out.println(vars[i].name+"preN="+vars[i].name+";");
	    //else
	    //    out.println(vars[i].name+"preN=f"+vars[i].name+"(0);");
	    out.println("}");
	    
	    for(int i=1;i<varNum;i++) 
		if(vars[i].init>=0)
		    out.println(vars[i].name+"="+vars[i].name+"p;");
		else
		    out.println(vars[i].name+"=f"+vars[i].name+"(0);");
	    
	    out.println("}");
	    out.println("}");
	    
	    out.println("");
	    
	    // output
	    out.println("");
	    out.println("// output");
	    out.println("FILE *out;");
	    out.println("char buffer[256];");
	    out.println("snprintf(buffer, sizeof(buffer), \"dummy.txt\");");
	    out.println("int idx=0;");
	    out.println("");
	    out.println("for(tb=0;tb<"+timeptNum+";tb++){");
	    for(int i=1;i<varNum;i++){
		out.println("snprintf(buffer, sizeof(buffer), \""+
			    outputDir+"/"+modelName+"CTx"+i+
			    "T%d_%d.txt\", tb,myid);");
		out.println("out=fopen(buffer, \"w\");");
		out.println("idx=0;");
		out.println("");
		out.println("");
		int tmp_p=(vars[i].parIDs).length;
		int tmp_v=+(vars[i].varIDs).length;
		for(int j=0;j<tmp_p;j++) 
		    out.println("for(int ki"+j+"=0;ki"+j+"<"
				+pars[vars[i].parIDs[j]].binNum+";ki"+j+"++)");
		for(int j=0;j<tmp_v;j++) 
		    out.println("for(int vi"+j+"=0;vi"+j+"<"
				+vars[vars[i].varIDs[j]].binNum+";vi"+j+"++)");
		out.println("for(int vi=0;vi<"+vars[i].binNum+";vi++)");
		out.println("{");
		
		out.print("int ctrtmp=("+vars[i].name+"ctr[tb]");
		for(int j=0;j<tmp_p;j++) out.print("[ki"+j+"]");
		for(int j=0;j<tmp_v;j++) out.print("[vi"+j+"]");
		out.println("[vi]);");
		out.println("if(ctrtmp>0){");
		out.println("fprintf(out,\"%d %d\\n\",idx,ctrtmp);");
		out.println("}");
		out.println("idx++;");
		out.println("}");
		out.println("fclose(out);");	
	    }

	    out.println("}");
	    out.println("");
	    out.println("return 0;");

	    out.println("}");

	    out.flush();
	    
	    
	} catch(IOException e){
	    e.printStackTrace();
	}
    }




    //------------- generating Batch code -------------------
    // new sampling
    
    public void genBATcode3(int sampleSize,int timeptNum,double dt,double simTime,String dSFMTPath,String outputDir){
	try{
	    FileOutputStream outfile=new FileOutputStream(modelName+"bat.c");
	    PrintWriter out=new PrintWriter(outfile);

	    // header
	    out.println("#include <stdio.h>");
	    out.println("#include <math.h>");
	    out.println("#include \""+dSFMTPath+"\"");
	    out.println("");

	    // initilization
	    out.println("");
	    out.println("");
	    out.println("");
	    out.println("int sampleNo="+sampleSize+";");
	    out.println("dsfmt_t dsfmt;");
	    out.println("");
	    out.println("");

	    // global variable
	    out.println("int tps,t,i,block,tb;");
	    out.println("double dt,halfdt,halfF1,halfF2,F3,F4;");

	    // --------- variable 
	    for(int i=1;i<varNum;i++){
		out.print("double "+vars[i].name+"i[]={");
		double tmp_r=vars[i].upperBound-vars[i].lowerBound;
		double tmp_b=vars[i].lowerBound;
		for(int j=0;j<vars[i].binNum;j++){
		    out.print(tmp_b+",");
		    tmp_b+=vars[i].binSizes[j]*tmp_r;
		}
		out.println(vars[i].upperBound+"};");
	    }
	    out.println("");  
	    for(int i=1;i<varNum;i++) 
		out.println("int "+vars[i].name+"pre,"+vars[i].name+
			    "post,"+vars[i].name+"init="+vars[i].binInit+";");
	    // --------- parameter
	    for(int i=1;i<parNum;i++){
		out.print("double "+pars[i].name+"i[]={");
		double tmp_r=pars[i].upperBound-pars[i].lowerBound;
		double tmp_b=pars[i].lowerBound;
		for(int j=0;j<pars[i].binNum;j++){
		    out.print(tmp_b+",");
		    tmp_b+=pars[i].binSizes[j]*tmp_r;
		}
		out.println(pars[i].upperBound+"};");
	    }
	    for(int i=1;i<parNum;i++) 
		out.println("int "+pars[i].name+"bin;");


	    for(int i=1;i<varNum;i++) 
		out.println("double "+vars[i].name+","+vars[i].name+"p,"+vars[i].name+"preN;");
	    for(int i=1;i<parNum;i++)
		out.println("double "+pars[i].name+";");
	    for(int i=1;i<varNum;i++){
		int n=(vars[i].varIDs).length+(vars[i].parIDs).length+1;
		out.print("int "+vars[i].name+"ctr["+timeptNum+"]");
		//for(int j=0;j<n;j++) out.print("["+vars[i].binNum+"]"); // bug july13
		int tmp_v=(vars[i].varIDs).length;
		int tmp_p=(vars[i].parIDs).length;
		for(int j=0;j<tmp_p;j++)
		    out.print("["+pars[(vars[i].parIDs)[j]].binNum+"]");
		for(int j=0;j<tmp_v;j++)
		    out.print("["+vars[(vars[i].varIDs)[j]].binNum+"]");
		out.print("["+vars[i].binNum+"]");
		out.println(";");
	    }
	    out.println("");

	    // functions
	    for(int i=1;i<varNum;i++){
		StringTokenizer st = new StringTokenizer(vars[i].eqn, "=");
		String name=st.nextToken(); //name
		String body=st.nextToken(); // body
		out.println("double f"+name+"(double "+name+"){");
		out.println("return "+body+";");
		out.println("}");	
	    }
	    out.println("");

	    // discretization function
	    out.println("");
	    out.println("int discretize(double v, double xi[], int length){");
	    out.println("for(int j=1;j<length-1;j++) if(v<xi[j]) return j-1;");
	    out.println("return length-2;");
	    out.println("}");
	    out.println("");


	    // sampling function
	    out.println("");
	    out.println("void sampling(){");
	    for(int i=1;i<parNum;i++){
		out.println(pars[i].name+"="+pars[i].lowerBound+
			    "+dsfmt_genrand_close_open(&dsfmt)*"+
			    (pars[i].upperBound-pars[i].lowerBound)+";");
	    }
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0)
		    out.println(vars[i].name+"="+vars[i].name+"i["+vars[i].name+
				"init]+dsfmt_genrand_close_open(&dsfmt)*("+
				vars[i].name+"i["+vars[i].name+"init+1]-"+
				vars[i].name+"i["+vars[i].name+"init]);");
		else
		    out.println(vars[i].name+"=f"+vars[i].name+"(0);");
	    }
	    out.println("}");
	    out.println("");


	    // solve function
	    out.println("");
	    out.println("void solve(){");
	    for(int i=1;i<varNum;i++){
		//if(vars[i].init>=0)
		out.println(vars[i].name+"preN="+vars[i].name+";");
		    //else
		    //out.println(vars[i].name+"preN=f"+vars[i].name+"(0);");
	    }
	    for(int i=1;i<parNum;i++) 
		out.println(pars[i].name+"bin=discretize("+pars[i].name+","+pars[i].name+"i,"+(pars[i].binNum+1)+");");
	    
	    out.println("");

	    // ODE solver
	    out.println("for(int t=1;t<=tps;t++){");
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0){
		    out.println("// "+vars[i].name);
		    out.println("halfF1=halfdt*f"+vars[i].name+"("+vars[i].name+");");
		    out.println("halfF2=halfdt*f"+vars[i].name+"("+vars[i].name+"+halfF1);");
		    out.println("F3=dt*f"+vars[i].name+"("+vars[i].name+"+halfF2);");
		    out.println("F4=dt*f"+vars[i].name+"("+vars[i].name+"+F3);");
		    out.println(vars[i].name+"p="+vars[i].name+"+(2*halfF1+4*halfF2+2*F3+F4)/6.0;");
		}
	    }
	    
	    out.println("if(t%block==0){");
	    out.println("tb=t/block-1;");
	    
	    for(int i=1;i<varNum;i++) 
		out.println(vars[i].name+"pre=discretize("+
			    vars[i].name+"preN,"+vars[i].name+"i,"+(vars[i].binNum+1)+");");
	    out.println("");

	    for(int i=1;i<varNum;i++) 
		if(vars[i].init>=0)
		    out.println(vars[i].name+"post=discretize("+
				vars[i].name+","+vars[i].name+"i,"+(vars[i].binNum+1)+");");
	    out.println("");

	    out.println("");
	    for(int i=1;i<varNum;i++){
		out.print(vars[i].name+"ctr[tb]");
		for(int j=0;j<(vars[i].parIDs).length;j++) 
		    out.print("["+pars[vars[i].parIDs[j]].name+"bin]");
		for(int j=0;j<(vars[i].varIDs).length;j++) 
		    out.print("["+vars[vars[i].varIDs[j]].name+"pre]");
		if(vars[i].init>=0)
		    out.println("["+vars[i].name+"post]++;");
		else
		    out.println("["+vars[i].name+"pre]++;");
	    }
	    for(int i=1;i<varNum;i++)
		//if(vars[i].init>=0)
		out.println(vars[i].name+"preN="+vars[i].name+";");
	    //else
	    //    out.println(vars[i].name+"preN=f"+vars[i].name+"(0);");
	    out.println("}");
	    
	    for(int i=1;i<varNum;i++) 
		if(vars[i].init>=0)
		    out.println(vars[i].name+"="+vars[i].name+"p;");
		else
		    out.println(vars[i].name+"=f"+vars[i].name+"(0);");
	    
	    out.println("}");
	    out.println("}");
	    //out.println("}");
	    out.println("");


	    // main
	    out.println("int main(int argc, char *argv[]){");
	    out.println("int myid=atoi(argv[1]);");
	    out.println("dt="+dt+",halfdt=dt/2.0;");
	    out.println("tps=(int)("+simTime+"/dt);");
	    out.println("block=tps/"+timeptNum+";");
	    
 	    Random seed=new Random();
	    out.println("");
	    out.println("");
	    out.println("");
 	    out.println("int seed="+seed.nextInt(12345)+"+myid;");
	    out.println("dsfmt_init_gen_rand(&dsfmt,seed);");
	    
	    // sampling
 	    out.println("for(int i=0;i<sampleNo;i++){");
	    for(int i=1;i<varNum;i++){
		out.println("// "+vars[i].name+" samples");
		int tmp_p=(vars[i].parIDs).length;
		for(int j=0;j<tmp_p;j++) 
		    out.println("for(int ki"+j+"=0;ki"+j+"<"
				+pars[vars[i].parIDs[j]].binNum+";ki"+j+"++)");
		out.println("{");
		out.println("sampling();");
		out.println(""); 
		for(int j=0;j<tmp_p;j++)
		    out.println(pars[vars[i].parIDs[j]].name+"="+pars[vars[i].parIDs[j]].name+"i[ki"+j+
				"]+dsfmt_genrand_close_open(&dsfmt)*("+
				pars[vars[i].parIDs[j]].name+"i[ki"+j+"+1]-"+
				pars[vars[i].parIDs[j]].name+"i[ki"+j+"]);");
		
		out.println(""); 
		out.println("solve();");
		out.println("}");
		out.println("");
	    }
	    out.println("}");
	    
	    out.println("");
	    
	    // output
	    out.println("");
	    out.println("// output");
	    out.println("FILE *out;");
	    out.println("char buffer[256];");
	    out.println("snprintf(buffer, sizeof(buffer), \"dummy.txt\");");
	    out.println("int idx=0;");
	    out.println("");
	    out.println("for(tb=0;tb<"+timeptNum+";tb++){");
	    for(int i=1;i<varNum;i++){
		out.println("snprintf(buffer, sizeof(buffer), \""+
			    outputDir+"/"+modelName+"CTx"+i+
			    "T%d_%d.txt\", tb,myid);");
		out.println("out=fopen(buffer, \"w\");");
		out.println("idx=0;");
		out.println("");
		out.println("");
		int tmp_p=(vars[i].parIDs).length;
		int tmp_v=+(vars[i].varIDs).length;
		for(int j=0;j<tmp_p;j++) 
		    out.println("for(int ki"+j+"=0;ki"+j+"<"
				+pars[vars[i].parIDs[j]].binNum+";ki"+j+"++)");
		for(int j=0;j<tmp_v;j++) 
		    out.println("for(int vi"+j+"=0;vi"+j+"<"
				+vars[vars[i].varIDs[j]].binNum+";vi"+j+"++)");
		out.println("for(int vi=0;vi<"+vars[i].binNum+";vi++)");
		out.println("{");
		
		out.print("int ctrtmp=("+vars[i].name+"ctr[tb]");
		for(int j=0;j<tmp_p;j++) out.print("[ki"+j+"]");
		for(int j=0;j<tmp_v;j++) out.print("[vi"+j+"]");
		out.println("[vi]);");
		out.println("if(ctrtmp>0){");
		out.println("fprintf(out,\"%d %d\\n\",idx,ctrtmp);");
		out.println("}");
		out.println("idx++;");
		out.println("}");
		out.println("fclose(out);");	
	    }

	    out.println("}");
	    out.println("");
	    out.println("return 0;");

	    out.println("}");

	    out.flush();
	    
	    
	} catch(IOException e){
	    e.printStackTrace();
	}
    }




    




    //------------- processing results on cluster and output probabilities -------------------

    public void process(int nodeNum,int timeptNum,String outputDir,int option){
	try{
	    int ctr;
	    String line;
	    FileInputStream fis;
	    BufferedReader stdin;
	    FileOutputStream outfile;
	    PrintWriter out;
	    StringTokenizer st;
	    for(int tb=0;tb<timeptNum;tb++){
		for(int x=1;x<varNum;x++){
		    // combine counts on N nodes
		    int row=vars[x].binNum;
		    for(int i=0;i<(vars[x].parIDs).length;i++)
			row*=pars[(vars[x].parIDs)[i]].binNum;
		    for(int i=0;i<(vars[x].varIDs).length;i++)
			row*=vars[(vars[x].varIDs)[i]].binNum;
		    int []count=new int[row];
		    outfile=new FileOutputStream(outputDir+"/"+modelName+"CTx"+x+"T"+tb+".txt");
		    out = new PrintWriter(outfile);
		    for(int i=0;i<nodeNum;i++){
			fis = new FileInputStream(outputDir+"/"+modelName+"CTx"+x+"T"+tb+"_"+i+".txt");
			stdin= new BufferedReader(new InputStreamReader(fis));
			while((line=stdin.readLine())!=null){
			    st = new StringTokenizer(line, " ");
			    count[Integer.parseInt(st.nextToken())]+=Integer.parseInt(st.nextToken());
			}
			stdin.close();fis.close();
		    }
		    // output total counts
		    for(int j=0;j<row;j++)
			if(count[j]>0) out.println(j+" "+count[j]);
		    out.flush();out.close();
		    
		    
		    // compute and output CPT 
		    switch(option){
		    case 0: { // option 0 --- use high-dim arrays to store CPT
			outfile=new FileOutputStream(outputDir+"/tables/"+modelName+"Px"+x+"T"+tb+".txt");
			out = new PrintWriter(outfile);
			int B=vars[x].binNum;
			int pa=row/B;
			for(int p=0;p<pa;p++){
			    int sum=0;
			    for(int j=0;j<B;j++)
				sum+=count[B*p+j];
			    if(sum>0){
				for(int j=0;j<B;j++)
				    out.println((p*B+j)+" "+count[p*B+j]*1.0/sum);	
			    }
			}
			out.flush();out.close();break;
		    }
		    case 1: { // option 1 --- use hashtable to store CPT
			int B=vars[x].binNum;
			int pa=row/B;
			String key="x"+x+"t"+tb+"i";
			for(int p=0;p<pa;p++){
			    int sum=0;
			    for(int j=0;j<B;j++)
				sum+=count[B*p+j];
			    if(sum>0){
				for(int j=0;j<B;j++){
				    key=key+(p*B+j);
				    cpt.put(key,count[p*B+j]*1.0/sum);
				}	
			    }
			}
			break;}
		    default: System.out.println("No tables."); break;
		    }
		}

	    }
	    
	    if(option==1){
		FileOutputStream fos = new FileOutputStream(outputDir+"/tables/"+modelName+"_cpt.dat");
		ObjectOutputStream oos = new ObjectOutputStream(fos);
		oos.writeObject(cpt);
		oos.close();
	    }
	    
	    
	    
	} catch(IOException e){
	    e.printStackTrace();
	}
    }

    //------------- processing results on cluster and output probabilities -------------------
    /* e.g. inputDir='/home/course/cs1101s/m201/'  // then search .../m201/ct0, ct1, ct2...
            outputDir='/home/course/cs1101s/m201/ctn' // store results in .../ctn/tables
     */

    public void processN(int dirNum,int nodeNum,int timeptNum,String inputDir,String outputDir,int option){
	try{
	    int ctr;
	    String line;
	    FileInputStream fis;
	    BufferedReader stdin;
	    FileOutputStream outfile;
	    PrintWriter out;
	    StringTokenizer st;
	    String fileOutPrefix=outputDir+"/"+modelName+"CT";
	    String tablesDir=outputDir+"/tables";
	    new File(tablesDir).mkdirs();
	    int [][]hitmap=new int[dirNum][nodeNum+1];
	    for(int dir=0;dir<dirNum;dir++){
		if(new File(inputDir+"/ct"+dir+"/"+modelName+"CTx1T0.txt").exists())
		    hitmap[dir][nodeNum]=1;
		else{
		    for(int node=0;node<nodeNum;node++){
			if(new File(inputDir+"/ct"+dir+"/"+modelName+"CTx1T0_"+node+".txt").exists())
			    hitmap[dir][node]=1;
			//System.out.print(hitmap[dir][node]+" ");
		    }
		    //System.out.println();
		}
	    }

	    String tag="P";
	    if(modelName.equals("brown")) tag="XBNT"; // [patch] brown's model use XBNT lable

	    
	    for(int tb=0;tb<timeptNum;tb++){
		for(int x=1;x<varNum;x++){
		    // combine counts on N nodes
		    int row=vars[x].binNum;
		    for(int i=0;i<(vars[x].parIDs).length;i++)
			row*=pars[(vars[x].parIDs)[i]].binNum;
		    for(int i=0;i<(vars[x].varIDs).length;i++)
			row*=vars[(vars[x].varIDs)[i]].binNum;
		    int []count=new int[row];
		    outfile=new FileOutputStream(fileOutPrefix+"x"+x+"T"+tb+".txt");
		    out = new PrintWriter(outfile);
		    for(int dir=0;dir<dirNum;dir++){
			for(int i=0;i<=nodeNum;i++){
			    if(hitmap[dir][i]==1){
				if(i==nodeNum) fis = new FileInputStream(inputDir+"/ct"+dir+"/"+modelName+"CTx"+x+"T"+tb+".txt");
				else fis = new FileInputStream(inputDir+"/ct"+dir+"/"+modelName+"CTx"+x+"T"+tb+"_"+i+".txt");
				stdin= new BufferedReader(new InputStreamReader(fis));
				while((line=stdin.readLine())!=null){
				    st = new StringTokenizer(line, " ");
				    count[Integer.parseInt(st.nextToken())]+=Integer.parseInt(st.nextToken());
				}
				stdin.close();fis.close();
			    }
			}
		    }
		    // output total counts
		    for(int j=0;j<row;j++)
			if(count[j]>0) out.println(j+" "+count[j]);
		    out.flush();out.close();
		    
		    
		    // compute and output CPT 
		    switch(option){
		    case 0: { // option 0 --- use high-dim arrays to store CPT
			outfile=new FileOutputStream(outputDir+"/tables/"+modelName+""+tag+"x"+x+"T"+tb+".txt");
			out = new PrintWriter(outfile);
			int B=vars[x].binNum;
			int pa=row/B;
			for(int p=0;p<pa;p++){
			    int sum=0;
			    for(int j=0;j<B;j++)
				sum+=count[B*p+j];
			    if(sum>0){
				for(int j=0;j<B;j++)
				    if(count[p*B+j]>0) // [update] do not output 0; date 072209 
					out.println((p*B+j)+" "+count[p*B+j]*1.0/sum);	
			    }
			}
			out.flush();out.close();break;
		    }
		    case 1: { // option 1 --- use hashtable to store CPT
			int B=vars[x].binNum;
			int pa=row/B;
			String key="x"+x+"t"+tb+"i";
			for(int p=0;p<pa;p++){
			    int sum=0;
			    for(int j=0;j<B;j++)
				sum+=count[B*p+j];
			    if(sum>0){
				for(int j=0;j<B;j++){
				    key=key+(p*B+j);
				    cpt.put(key,count[p*B+j]*1.0/sum);
				}	
			    }
			}
			break;}
		    default: System.out.println("Invalid option."); break;
		    }
		}

	    }
	    
	    if(option==1){
		FileOutputStream fos = new FileOutputStream(outputDir+"/tables/"+modelName+"_cpt.dat");
		ObjectOutputStream oos = new ObjectOutputStream(fos);
		oos.writeObject(cpt);
		oos.close();
	    }
	    
	    
	    
	} catch(IOException e){
	    e.printStackTrace();
	}
    }



    //------------- population-based ODE simulation (Runge-Kutta) --------------    
    
    public void genNORMcode(int sampleSize,int timeptNum,double dt,double simTime,String dSFMTPath,String outputDir){
	try{
	    FileOutputStream outfile=new FileOutputStream(modelName+"Norm.c");
	    PrintWriter out=new PrintWriter(outfile);
	    
	    // header
	    out.println("#include <stdio.h>");
	    out.println("#include <math.h>");
	    out.println("#include \""+dSFMTPath+"\"");
	    out.println("");
	    
	    // global variable
	    for(int i=1;i<varNum;i++) 
		out.println("double "+vars[i].name+","+vars[i].name+"p,"+vars[i].name+"preN;");
	    for(int i=1;i<parNum;i++)
		out.println("double "+pars[i].name+";");
	    out.println("double resultMean["+(timeptNum+1)+"]["+varNum+"];");
	    out.println("double resultMin["+(timeptNum+1)+"]["+varNum+"];");
	    out.println("double resultMax["+(timeptNum+1)+"]["+varNum+"];");
	    
	    out.println("");
	    
	    // functions
	    for(int i=1;i<varNum;i++){
		StringTokenizer st = new StringTokenizer(vars[i].eqn, "=");
		String name=st.nextToken(); //name
		String body=st.nextToken(); // body
		out.println("double f"+name+"(double "+name+"){");
		out.println("return "+body+";");
		out.println("}");	
	    }
	    out.println("");

	    // main
	    out.println("int main(int argc, char *argv[]){");
	    out.println("double dt="+dt+",halfdt=dt/2.0;");
	    out.println("int tps=(int)("+simTime+"/dt),t,i;");
	    out.println("int block=tps/"+timeptNum+",tb;");
	    out.println("double halfF1,halfF2,F3,F4;");

	    // -------- results initialization
	    out.println("for(int tb=0;tb<"+(timeptNum+1)+";tb++)");
	    out.println("for(int i=0;i<"+varNum+";i++){");
	    out.println("resultMean[tb][i]=0;");
	    out.println("resultMin[tb][i]=10000;");
	    out.println("resultMax[tb][i]=0;");
	    out.println("}");
	    out.println("");

	    // --------- variable 
	    for(int i=1;i<varNum;i++){
		out.print("double "+vars[i].name+"i[]={");
		double tmp_r=vars[i].upperBound-vars[i].lowerBound;
		double tmp_b=vars[i].lowerBound;
		for(int j=0;j<vars[i].binNum;j++){
		    out.print(tmp_b+",");
		    tmp_b+=vars[i].binSizes[j]*tmp_r;
		}
		out.println(vars[i].upperBound+"};");
	    }
	    out.println("");  
	    for(int i=1;i<varNum;i++) 
		out.println("int "+vars[i].name+"init="+vars[i].binInit+";");
	    // --------- parameter
	    for(int i=1;i<parNum;i++){
		out.print("double "+pars[i].name+"i[]={");
		double tmp_r=pars[i].upperBound-pars[i].lowerBound;
		double tmp_b=pars[i].lowerBound;
		for(int j=0;j<pars[i].binNum;j++){
		    out.print(tmp_b+",");
		    tmp_b+=pars[i].binSizes[j]*tmp_r;
		}
		out.println(pars[i].upperBound+"};");
	    }
	    for(int i=1;i<parNum;i++) 
		out.println("int "+pars[i].name+"init="+pars[i].binInit+";");
	    Random seed=new Random();
	    out.println("");
	    out.println("");
	    out.println("");
	    out.println("");
	    out.println("int sampleNo="+sampleSize+";");
	    out.println("dsfmt_t dsfmt;");
	    out.println("int seed="+seed.nextInt(12345)+";");
	    out.println("dsfmt_init_gen_rand(&dsfmt,seed);");
	    out.println("");
	    out.println("");
	    out.println("for(int i=0;i<sampleNo;i++){");
	    for(int i=1;i<parNum;i++){ // updated Apr 15, 2009
		out.println(pars[i].name+"="+pars[i].name+"i["+pars[i].name+
			    "init]+dsfmt_genrand_close_open(&dsfmt)*("+
			    pars[i].name+"i["+pars[i].name+"init+1]-"+
			    pars[i].name+"i["+pars[i].name+"init]);");
	    }
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0)
		    ///*
		    out.println(vars[i].name+"="+vars[i].name+"i["+vars[i].name+
				"init]+dsfmt_genrand_close_open(&dsfmt)*("+
				vars[i].name+"i["+vars[i].name+"init+1]-"+
				vars[i].name+"i["+vars[i].name+"init]);");
		    //*/
		    //out.println(vars[i].name+"="+vars[i].init+";");
		else
		    out.println(vars[i].name+"=f"+vars[i].name+"(0);");
	    }
	    out.println("");

	    // ODE solver
	    out.println("for(int t=1;t<=tps+1;t++){");
	    out.println("if(t%block==1){ // collect results");
	    out.println("tb=t/block;");
	    for(int i=1;i<varNum;i++){ //collect results
		out.println("resultMean[tb]["+i+"]+="+vars[i].name+";");
		out.println("if("+vars[i].name+"<resultMin[tb]["+i+"]) resultMin[tb]["+i+"]="+vars[i].name+";");
		out.println("if("+vars[i].name+">resultMax[tb]["+i+"]) resultMax[tb]["+i+"]="+vars[i].name+";");
	    }
	    out.println("}");

	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0){
		    out.println("// "+vars[i].name);
		    out.println("halfF1=halfdt*f"+vars[i].name+"("+vars[i].name+");");
		    out.println("halfF2=halfdt*f"+vars[i].name+"("+vars[i].name+"+halfF1);");
		    out.println("F3=dt*f"+vars[i].name+"("+vars[i].name+"+halfF2);");
		    out.println("F4=dt*f"+vars[i].name+"("+vars[i].name+"+F3);");
		    out.println(vars[i].name+"p="+vars[i].name+"+(2*halfF1+4*halfF2+2*F3+F4)/6.0;");
		}
	    }
	    
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0)
		    out.println(vars[i].name+"="+vars[i].name+"p;");
		else
		    out.println(vars[i].name+"=f"+vars[i].name+"(0);");
	    }
	    out.println("}");
	    out.println("}");
	    
	    out.println("");
	    
	    // output
	    out.println("");
	    out.println("// output");
	    out.println("FILE *out;");
	    out.println("char buffer[256];");
	    out.println("snprintf(buffer, sizeof(buffer), \"dummy.txt\");");
	    out.println("int idx=0;");
	    out.println("");
	    
	    out.println("snprintf(buffer, sizeof(buffer), \""+outputDir+"/"+modelName+"Mean.txt\");");
	    out.println("out=fopen(buffer, \"w\");");
	    out.println("fprintf(out,\"Time,\");"); // header
	    out.println("for(int i=1;i<"+varNum+";i++)");
	    out.println("fprintf(out,\"x%d,\",i);");
	    out.println("fprintf(out,\"\\n\");");
	    out.println("for(int tb=0;tb<"+(timeptNum+1)+";tb++){");
	    out.println("fprintf(out,\"%d,\",tb);");
	    out.println("for(int i=1;i<"+varNum+";i++)");
	    out.println("fprintf(out,\"%f,\",resultMean[tb][i]/sampleNo);");
	    out.println("fprintf(out,\"\\n\");");
	    out.println("}");
	    out.println("fclose(out);");

	    
	    out.println("snprintf(buffer, sizeof(buffer), \""+outputDir+"/"+modelName+"Min.txt\");");
	    out.println("out=fopen(buffer, \"w\");");
	    out.println("fprintf(out,\"Time,\");"); // header
	    out.println("for(int i=1;i<"+varNum+";i++)");
	    out.println("fprintf(out,\"x%d,\",i);");
	    out.println("fprintf(out,\"\\n\");");
	    out.println("for(int tb=0;tb<"+(timeptNum+1)+";tb++){");
	    out.println("fprintf(out,\"%d,\",tb);");
	    out.println("for(int i=1;i<"+varNum+";i++)");
	    out.println("fprintf(out,\"%f,\",resultMin[tb][i]);");
	    out.println("fprintf(out,\"\\n\");");
	    out.println("}");
	    out.println("fclose(out);");

	    
	    out.println("snprintf(buffer, sizeof(buffer), \""+outputDir+"/"+modelName+"Max.txt\");");
	    out.println("out=fopen(buffer, \"w\");");
	    out.println("fprintf(out,\"Time,\");"); // header
	    out.println("for(int i=1;i<"+varNum+";i++)");
	    out.println("fprintf(out,\"x%d,\",i);");
	    out.println("fprintf(out,\"\\n\");");
	    out.println("for(int tb=0;tb<"+(timeptNum+1)+";tb++){");
	    out.println("fprintf(out,\"%d,\",tb);");
	    out.println("for(int i=1;i<"+varNum+";i++)");
	    out.println("fprintf(out,\"%f,\",resultMax[tb][i]);");
	    out.println("fprintf(out,\"\\n\");");
	    out.println("}");
	    out.println("fclose(out);");

	    
	    out.println("return 0;");
	    
	    out.println("}");
	    
	    out.flush();
	    
	    
	} catch(IOException e){
	    e.printStackTrace();
	}
    }
    

    //------------- population-based ODE simulation (Runge-Kutta) --------------    
    // parameter-less

    public void genNORMcode2(int sampleSize,int timeptNum,double dt,double simTime,String dSFMTPath,String outputDir){
	try{
	    FileOutputStream outfile=new FileOutputStream(modelName+"Norm2.c");
	    PrintWriter out=new PrintWriter(outfile);
	    
	    // header
	    out.println("#include <stdio.h>");
	    out.println("#include <math.h>");
	    out.println("#include \""+dSFMTPath+"\"");
	    out.println("");
	    
	    // global variable
	    for(int i=1;i<varNum;i++) 
		out.println("double "+vars[i].name+","+vars[i].name+"p,"+vars[i].name+"preN;");
	    for(int i=1;i<parNum;i++)
		out.println("double "+pars[i].name+";");
	    out.println("double resultMean["+(timeptNum+1)+"]["+varNum+"];");
	    out.println("double resultMin["+(timeptNum+1)+"]["+varNum+"];");
	    out.println("double resultMax["+(timeptNum+1)+"]["+varNum+"];");
	    
	    out.println("");
	    
	    // functions
	    for(int i=1;i<varNum;i++){
		StringTokenizer st = new StringTokenizer(vars[i].eqn, "=");
		String name=st.nextToken(); //name
		String body=st.nextToken(); // body
		out.println("double f"+name+"(double "+name+"){");
		out.println("return "+body+";");
		out.println("}");	
	    }
	    out.println("");

	    // main
	    out.println("int main(int argc, char *argv[]){");
	    out.println("double dt="+dt+",halfdt=dt/2.0;");
	    out.println("int tps=(int)("+simTime+"/dt),t,i;");
	    out.println("int block=tps/"+timeptNum+",tb;");
	    out.println("double halfF1,halfF2,F3,F4;");

	    // -------- results initialization
	    out.println("for(int tb=0;tb<"+(timeptNum+1)+";tb++)");
	    out.println("for(int i=0;i<"+varNum+";i++){");
	    out.println("resultMean[tb][i]=0;");
	    out.println("resultMin[tb][i]=10000;");
	    out.println("resultMax[tb][i]=0;");
	    out.println("}");
	    out.println("");

	    // --------- variable 
	    for(int i=1;i<varNum;i++){
		out.print("double "+vars[i].name+"i[]={");
		double tmp_r=vars[i].upperBound-vars[i].lowerBound;
		double tmp_b=vars[i].lowerBound;
		for(int j=0;j<vars[i].binNum;j++){
		    out.print(tmp_b+",");
		    tmp_b+=vars[i].binSizes[j]*tmp_r;
		}
		out.println(vars[i].upperBound+"};");
	    }
	    out.println("");  
	    for(int i=1;i<varNum;i++) 
		out.println("int "+vars[i].name+"init="+vars[i].binInit+";");
	    // --------- parameter
	    for(int i=1;i<parNum;i++){
		out.print("double "+pars[i].name+"i[]={");
		double tmp_r=pars[i].upperBound-pars[i].lowerBound;
		double tmp_b=pars[i].lowerBound;
		for(int j=0;j<pars[i].binNum;j++){
		    out.print(tmp_b+",");
		    tmp_b+=pars[i].binSizes[j]*tmp_r;
		}
		out.println(pars[i].upperBound+"};");
	    }
	    Random seed=new Random();
	    out.println("");
	    out.println("");
	    out.println("");
	    out.println("");
	    out.println("int sampleNo="+sampleSize+";");
	    out.println("dsfmt_t dsfmt;");
	    out.println("int seed="+seed.nextInt(12345)+";");
	    out.println("dsfmt_init_gen_rand(&dsfmt,seed);");
	    out.println("");
	    out.println("");
	    out.println("for(int i=0;i<sampleNo;i++){");
	    for(int i=1;i<parNum;i++){
		out.println(pars[i].name+"="+pars[i].lowerBound+
			    "+dsfmt_genrand_close_open(&dsfmt)*"+
			    (pars[i].upperBound-pars[i].lowerBound)+";");
	    }
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0)
		    ///*
		    out.println(vars[i].name+"="+vars[i].name+"i["+vars[i].name+
				"init]+dsfmt_genrand_close_open(&dsfmt)*("+
				vars[i].name+"i["+vars[i].name+"init+1]-"+
				vars[i].name+"i["+vars[i].name+"init]);");
		    //*/
		    //out.println(vars[i].name+"="+vars[i].init+";");
		else
		    out.println(vars[i].name+"=f"+vars[i].name+"(0);");
	    }
	    out.println("");

	    // ODE solver
	    out.println("for(int t=1;t<=tps+1;t++){");
	    out.println("if(t%block==1){ // collect results");
	    out.println("tb=t/block;");
	    for(int i=1;i<varNum;i++){ //collect results
		out.println("resultMean[tb]["+i+"]+="+vars[i].name+";");
		out.println("if("+vars[i].name+"<resultMin[tb]["+i+"]) resultMin[tb]["+i+"]="+vars[i].name+";");
		out.println("if("+vars[i].name+">resultMax[tb]["+i+"]) resultMax[tb]["+i+"]="+vars[i].name+";");
	    }
	    out.println("}");

	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0){
		    out.println("// "+vars[i].name);
		    out.println("halfF1=halfdt*f"+vars[i].name+"("+vars[i].name+");");
		    out.println("halfF2=halfdt*f"+vars[i].name+"("+vars[i].name+"+halfF1);");
		    out.println("F3=dt*f"+vars[i].name+"("+vars[i].name+"+halfF2);");
		    out.println("F4=dt*f"+vars[i].name+"("+vars[i].name+"+F3);");
		    out.println(vars[i].name+"p="+vars[i].name+"+(2*halfF1+4*halfF2+2*F3+F4)/6.0;");
		}
	    }
	    
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0)
		    out.println(vars[i].name+"="+vars[i].name+"p;");
		else
		    out.println(vars[i].name+"=f"+vars[i].name+"(0);");
	    }
	    out.println("}");
	    out.println("}");
	    
	    out.println("");
	    
	    // output
	    out.println("");
	    out.println("// output");
	    out.println("FILE *out;");
	    out.println("char buffer[256];");
	    out.println("snprintf(buffer, sizeof(buffer), \"dummy.txt\");");
	    out.println("int idx=0;");
	    out.println("");
	    
	    out.println("snprintf(buffer, sizeof(buffer), \""+outputDir+"/"+modelName+"Mean2.txt\");");
	    out.println("out=fopen(buffer, \"w\");");
	    out.println("fprintf(out,\"Time,\");"); // header
	    out.println("for(int i=1;i<"+varNum+";i++)");
	    out.println("fprintf(out,\"x%d,\",i);");
	    out.println("fprintf(out,\"\\n\");");
	    out.println("for(int tb=0;tb<"+(timeptNum+1)+";tb++){");
	    out.println("fprintf(out,\"%d,\",tb);");
	    out.println("for(int i=1;i<"+varNum+";i++)");
	    out.println("fprintf(out,\"%f,\",resultMean[tb][i]/sampleNo);");
	    out.println("fprintf(out,\"\\n\");");
	    out.println("}");
	    out.println("fclose(out);");

	    
	    out.println("snprintf(buffer, sizeof(buffer), \""+outputDir+"/"+modelName+"Min2.txt\");");
	    out.println("out=fopen(buffer, \"w\");");
	    out.println("fprintf(out,\"Time,\");"); // header
	    out.println("for(int i=1;i<"+varNum+";i++)");
	    out.println("fprintf(out,\"x%d,\",i);");
	    out.println("fprintf(out,\"\\n\");");
	    out.println("for(int tb=0;tb<"+(timeptNum+1)+";tb++){");
	    out.println("fprintf(out,\"%d,\",tb);");
	    out.println("for(int i=1;i<"+varNum+";i++)");
	    out.println("fprintf(out,\"%f,\",resultMin[tb][i]);");
	    out.println("fprintf(out,\"\\n\");");
	    out.println("}");
	    out.println("fclose(out);");

	    
	    out.println("snprintf(buffer, sizeof(buffer), \""+outputDir+"/"+modelName+"Max2.txt\");");
	    out.println("out=fopen(buffer, \"w\");");
	    out.println("fprintf(out,\"Time,\");"); // header
	    out.println("for(int i=1;i<"+varNum+";i++)");
	    out.println("fprintf(out,\"x%d,\",i);");
	    out.println("fprintf(out,\"\\n\");");
	    out.println("for(int tb=0;tb<"+(timeptNum+1)+";tb++){");
	    out.println("fprintf(out,\"%d,\",tb);");
	    out.println("for(int i=1;i<"+varNum+";i++)");
	    out.println("fprintf(out,\"%f,\",resultMax[tb][i]);");
	    out.println("fprintf(out,\"\\n\");");
	    out.println("}");
	    out.println("fclose(out);");

	    
	    out.println("return 0;");
	    
	    out.println("}");
	    
	    out.flush();
	    
	    
	} catch(IOException e){
	    e.printStackTrace();
	}
    }

   //------------- population-based ODE simulation (Runge-Kutta) --------------    
    // MPSA (c++ code)

    public void genMPSAcode(int sampleSize,int timeptNum,double dt,double simTime,String dSFMTPath,String outputDir,int[]freeParIDs, int[]dataTPs,int[]dataVarIDs, double[]weights, String dataFile, int ktest, double threshold){
	try{
	    FileOutputStream outfile=new FileOutputStream(modelName+"MPSA.c");
	    PrintWriter out=new PrintWriter(outfile);
	    
	    // header
	    out.println("#include <stdio.h>");
	    out.println("#include <math.h>");
	    out.println("#include \""+dSFMTPath+"\"");
	    out.println("");
	    
	    // global variable
	    for(int i=1;i<varNum;i++) 
		out.println("double "+vars[i].name+","+vars[i].name+"p,"+vars[i].name+"preN;");
	    for(int i=1;i<parNum;i++)
		out.println("double "+pars[i].name+";");
	    out.println("double result["+(timeptNum+1)+"]["+varNum+"];");
	    //out.println("double resultMean["+(timeptNum+1)+"]["+varNum+"];");
	    //out.println("double resultMin["+(timeptNum+1)+"]["+varNum+"];");
	    //out.println("double resultMax["+(timeptNum+1)+"]["+varNum+"];");
	    
	    out.println("");

	    // parameter estimation data
	    out.println("");
	    out.print("int P[]={");
	    for(int i=0;i<freeParIDs.length;i++)
		out.print(freeParIDs[i]+",");
	    out.println("};");

	    out.print("int T[]={");
	    for(int i=0;i<dataTPs.length;i++)
		out.print(dataTPs[i]+",");
	    out.println("};");
	    
	    out.print("int X[]={");
	    for(int i=0;i<dataVarIDs.length;i++)
		out.print(dataVarIDs[i]+",");
	    out.println("};");

	    out.print("double weights[]={");
	    for(int i=0;i<weights.length;i++)
		out.print(weights[i]+",");
	    out.println("};");

	    out.print("int k[]={-1,");
	    for(int i=1;i<parNum;i++)
		out.print(pars[i].binInit+",");	    
	    out.println("};");

	    out.println("double data["+dataTPs.length+"]["+dataVarIDs.length+"];");
	    out.println("");
	    // -------------------


	    
	    // functions
	    for(int i=1;i<varNum;i++){
		StringTokenizer st = new StringTokenizer(vars[i].eqn, "=");
		String name=st.nextToken(); //name
		String body=st.nextToken(); // body
		out.println("double f"+name+"(double "+name+"){");
		out.println("return "+body+";");
		out.println("}");	
	    }
	    out.println("");

	    // main
	    out.println("int main(int argc, char *argv[]){");
	    out.println("double dt="+dt+",halfdt=dt/2.0;");
	    out.println("int tps=(int)("+simTime+"/dt),t,i;");
	    out.println("int block=tps/"+timeptNum+",tb;");
	    out.println("double halfF1,halfF2,F3,F4;");

	    // --------- variable 
	    for(int i=1;i<varNum;i++){
		out.print("double "+vars[i].name+"i[]={");
		double tmp_r=vars[i].upperBound-vars[i].lowerBound;
		double tmp_b=vars[i].lowerBound;
		for(int j=0;j<vars[i].binNum;j++){
		    out.print(tmp_b+",");
		    tmp_b+=vars[i].binSizes[j]*tmp_r;
		}
		out.println(vars[i].upperBound+"};");
	    }
	    out.println("");  
	    for(int i=1;i<varNum;i++) 
		out.println("int "+vars[i].name+"init="+vars[i].binInit+";");

	    // output
	    out.println("");
	    out.println("// output setting");
	    out.println("FILE *out;");
	    out.println("char buffer[256];");
	    out.println("snprintf(buffer, sizeof(buffer), \"Nk"+ktest+".txt\");");
	    out.println("out=fopen(buffer, \"w\");");


	    //loadExp
	    out.println("char *cp;");
	    out.println("FILE *in;");
	    out.println("snprintf(buffer, sizeof(buffer), \""+dataFile+"\");");
	    out.println("in = fopen(buffer, \"r\");");
	    out.println("fgets(buffer, sizeof(buffer), in); //header");
	    out.println("fgets(buffer, sizeof(buffer), in); //init");
	    out.println("for(int j=0;j<"+dataTPs.length+";j++){");
	    out.println("int ctr2=0;");
	    out.println("fgets(buffer, sizeof(buffer), in);");
	    out.println("cp = (char *)strtok(buffer, \",\"); // time");
	    out.println("for(int k=1;k<"+varNum+";k++){");
	    out.println("cp = (char *)strtok(NULL, \",\");");
	    out.println("double tmp = strtod(cp, NULL);");
	    out.println("if(k==X[ctr2]){");
	    out.println("data[j][ctr2]=tmp;");
	    out.println("ctr2++;");
	    out.println("}");
	    out.println("}");
	    out.println("}");
	    out.println("");
	    out.println("");

	    Random seed=new Random();
	    out.println("");
	    out.println("");
	    out.println("");
	    out.println("");
	    out.println("int sampleNo="+sampleSize+";");
	    out.println("dsfmt_t dsfmt;");
	    out.println("int seed="+seed.nextInt(12345)+";");
	    out.println("dsfmt_init_gen_rand(&dsfmt,seed);");
	    out.println("");
	    out.println("double sumsum=0,threshold="+threshold+",freq=0;");
	    out.println("int colNo=40;");
	    out.println("int ctr[colNo];");
	    out.println("for(int i=0;i<colNo;i++) ctr[i]=0;");
	    out.println("");
	    out.println("");
	    out.println("");

	    out.println("for(int i=0;i<sampleNo;i++){");
	    out.println("for(int b=0;b<colNo;b++){");
	    out.println("double value=0;");


	    for(int i=1;i<parNum;i++){
		out.println(pars[i].name+"="+pars[i].lowerBound+
			    "+dsfmt_genrand_close_open(&dsfmt)*"+
			    (pars[i].upperBound-pars[i].lowerBound)+";");
	    }
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0)
		    ///*
		    out.println(vars[i].name+"="+vars[i].name+"i["+vars[i].name+
				"init]+dsfmt_genrand_close_open(&dsfmt)*("+
				vars[i].name+"i["+vars[i].name+"init+1]-"+
				vars[i].name+"i["+vars[i].name+"init]);");
		    //*/
		    //out.println(vars[i].name+"="+vars[i].init+";");
		else
		    out.println(vars[i].name+"=f"+vars[i].name+"(0);");
	    }
	    out.println("");


	    out.println(pars[ktest].name+"="+pars[ktest].lowerBound+"+b*"+(pars[ktest].upperBound-pars[ktest].lowerBound)+"/(colNo*1.0)+0.00001;");
	    out.println("");
	    out.println("");


	    // ODE solver
	    out.println("for(int t=1;t<=tps+1;t++){");
	    out.println("if(t%block==1){ // collect results");
	    out.println("tb=t/block;");
	    for(int i=1;i<varNum;i++){ //collect results
		out.println("result[tb]["+i+"]="+vars[i].name+";");
	    }
	    out.println("}");

	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0){
		    out.println("// "+vars[i].name);
		    out.println("halfF1=halfdt*f"+vars[i].name+"("+vars[i].name+");");
		    out.println("halfF2=halfdt*f"+vars[i].name+"("+vars[i].name+"+halfF1);");
		    out.println("F3=dt*f"+vars[i].name+"("+vars[i].name+"+halfF2);");
		    out.println("F4=dt*f"+vars[i].name+"("+vars[i].name+"+F3);");
		    out.println(vars[i].name+"p="+vars[i].name+"+(2*halfF1+4*halfF2+2*F3+F4)/6.0;");
		}
	    }
	    
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0)
		    out.println(vars[i].name+"="+vars[i].name+"p;");
		else
		    out.println(vars[i].name+"=f"+vars[i].name+"(0);");
	    }
	    out.println("}");

	    out.println("// evaluate");
	    out.println("double err=0;");
	    out.println("for(int zz=0;zz<"+10+";zz++){");
	    out.println("for(int j=0;j<"+7+";j++){");
	    out.println("err=result[T[zz]][X[j]]-data[zz][j];");
	    out.println("value+=weights[j]*err*err;");
	    out.println("}");
	    out.println("}");
	    out.println("if(value<threshold){");
	    out.println("ctr[b]++;");
	    out.println("}");
	    out.println("freq=ctr[b]/(i*1.0);");
	    out.println("fprintf(out,\"%f\\t\",freq);");
	    out.println("}");
	    out.println("fprintf(out,\"\\n\");");
	    out.println("");
	    out.println("");
	    out.println("");
	    out.println("");
	    out.println("");

	    out.println("}");
	    
	    out.println("");
	    
	    out.println("fclose(out);");
	    out.println("return 0;");
	    
	    out.println("}");
	    
	    out.flush();
	    
	    
	} catch(IOException e){
	    e.printStackTrace();
	}
    }

    //------------- population-based ODE simulation (CVODE 1.0 solver) --------------    
    // MPSA (c++ code)

    public void genMPSAcode2(int sampleSize,String dSFMTPath,String outputDir,int[]freeParIDs, int[]dataTPs,int[]dataVarIDs, String dataFile, double threshold, int colNo, int normalize, int scale){
	try{
	    FileOutputStream outfile=new FileOutputStream(modelName+"MPSA2.C");
	    PrintWriter out=new PrintWriter(outfile);
	    
	    // header
	    out.println("#include <stdio.h>");
	    out.println("#include <math.h>");
	    out.println("#include \""+dSFMTPath+"\"");
	    out.println("/* header files and macros for CVODE */");
	    out.println("#include \"llnltyps.h\"");
	    out.println("#include \"cvode.h\"");
	    out.println("#include \"cvdense.h\"");
	    out.println("#include \"nvector.h\"");
	    out.println("#include \"dense.h\"");
	    out.println("#define Ith(v,i) N_VIth(v,i-1)");
	    out.println("#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)");
	    out.println("");
	    
	    // global variable
	    out.println("int NEQ;");
	    out.println("double RTOL, ATOL;");
	    out.println("double T0, T1, Tm;");

	    for(int i=1;i<varNum;i++) 
		out.println("double "+vars[i].name+";");
	    for(int i=1;i<parNum;i++)
		out.println("double "+pars[i].name+";");
	    out.println("double result["+dataVarIDs.length+"]["+dataVarIDs.length+"];");
	    out.println("double data["+dataVarIDs.length+"]["+dataTPs.length+"];");
	    out.println("int tn,xn;");
	    out.println("");

	    // parameter estimation data
	    out.println("");
	    out.print("double *P[]={");
	    for(int i=0;i<freeParIDs.length;i++)
		out.print("&"+pars[freeParIDs[i]].name+",");
	    out.println("};");

	    out.print("int T[]={");
	    for(int i=0;i<dataTPs.length;i++)
		out.print(dataTPs[i]+",");
	    out.println("};");
	    
	    out.print("int X[]={");
	    for(int i=0;i<dataVarIDs.length;i++)
		out.print(dataVarIDs[i]+",");
	    out.println("};");

	    out.println("");
	    out.println("");
	    // -------------------

	    out.println("static void difeq(integer N, real t, N_Vector y, N_Vector ydot, void *f_data);");
	    

	    // main
	    out.println("int main(int argc, char *argv[]){");

	    out.println("int ktest;");
	    out.println("if(argc>1) ktest=atoi(argv[1]); else {printf(\"input ktest!\\n\"); exit(0);}");

	    out.println("double lb["+freeParIDs.length+"];");
	    out.println("double ub["+freeParIDs.length+"];");

	    if(scale>0){  
		/* 
			scale is used to define the upper-lower bounds of parameters
			e.g. for parameter k, scale = 10 means the range of k is [k/10, 10k]
			if scale = 0, then we use the bound specified in the par.csv file.
		*/ 
		for(int i=0; i<freeParIDs.length; i++)
		    out.println("lb["+i+"]="+pars[freeParIDs[i]].init/(scale*1.0)
				+"; ub["+i+"]="+pars[freeParIDs[i]].init*scale+";");
	    // for(int i=0; i<freeVarIDs.length; i++)
// 		out.println("lb["+(i+freeParIDs.length)+"]="+vars[freeVarIDs[i]].lowerBound
// 				+"; ub["+(i+freeParIDs.length)+"]="+vars[freeVarIDs[i]].upperBound+";");

	    } else {
		for(int i=0; i<freeParIDs.length; i++)
		    out.println("lb["+i+"]="+pars[freeParIDs[i]].lowerBound
				+"; ub["+i+"]="+pars[freeParIDs[i]].upperBound+";");
		// for(int i=0; i<freeVarIDs.length; i++)
// 		    out.println("lb["+(i+freeParIDs.length)+"]="+vars[freeVarIDs[i]].lowerBound
// 				+"; ub["+(i+freeParIDs.length)+"]="+vars[freeVarIDs[i]].upperBound+";");
	    }
	    out.println("/* set ODE parameters */");
	    out.println("NEQ = "+(varNum-1)+";");
	    out.println("RTOL = 1e-6;");
	    out.println("ATOL = 1e-6;");
	    out.println("T0 = 0.0;");
	    out.println("Tm = "+dataTPs[dataTPs.length-1]+";");
	    out.println("tn="+dataTPs.length+";");
	    out.println("xn="+dataVarIDs.length+";");
	    out.println("");
	    out.println("");

	    // init
	    for(int i=1;i<varNum;i++) 
		out.println(vars[i].name+"="+vars[i].init+";");
	    out.println("");
	    for(int i=1;i<parNum;i++) 
		out.println(pars[i].name+"="+pars[i].init+";");
	    out.println("");


// 	    // --------- variable 
// 	    for(int i=1;i<varNum;i++){
// 		out.print("double "+vars[i].name+"i[]={");
// 		double tmp_r=vars[i].upperBound-vars[i].lowerBound;
// 		double tmp_b=vars[i].lowerBound;
// 		for(int j=0;j<vars[i].binNum;j++){
// 		    out.print(tmp_b+",");
// 		    tmp_b+=vars[i].binSizes[j]*tmp_r;
// 		}
// 		out.println(vars[i].upperBound+"};");
// 	    }
// 	    out.println("");  
// 	    for(int i=1;i<varNum;i++) 
// 		out.println("int "+vars[i].name+"init="+vars[i].binInit+";");

	    //loadExp
	    out.println("");
	    out.println("// load data");
	    out.println("char *cp;");
	    out.println("FILE *in;");
	    out.println("char buffer[256];");

	    out.println("snprintf(buffer, sizeof(buffer), \""+dataFile+"\");");
	    out.println("in = fopen(buffer, \"r\");");
	    out.println("fgets(buffer, sizeof(buffer), in); //header");
	    //out.println("fgets(buffer, sizeof(buffer), in); //init");
	    out.println("for(int j=0;j<"+dataTPs.length+";j++){");
	    out.println("int ctr2=0;");
	    out.println("fgets(buffer, sizeof(buffer), in);");
	    out.println("cp = (char *)strtok(buffer, \"\\t\"); // time");
	    out.println("for(int k=0;k<xn;k++){");
	    out.println("cp = (char *)strtok(NULL, \"\\t\");");
	    out.println("double tmp = strtod(cp, NULL);");
	    out.println("data[j][k]=tmp;");
	    out.println("}");
	    out.println("}");
	    out.println("");

	    // output
	    out.println("");
	    out.println("// output setting");
	    out.println("FILE *out;");
	    out.println("snprintf(buffer, sizeof(buffer), \"Nk%d.txt\",ktest);");
	    out.println("out=fopen(buffer, \"w\");");





	    Random seed=new Random();
	    out.println("");
	    out.println("");
	    out.println("");
	    out.println("");
	    out.println("int sampleNo="+sampleSize+";");
	    out.println("dsfmt_t dsfmt;");
	    out.println("int seed="+seed.nextInt(12345)+";");
	    out.println("dsfmt_init_gen_rand(&dsfmt,seed);");
	    out.println("");
	    out.println("double sumsum=0,threshold="+threshold+",freq=0;");
	    out.println("int colNo="+colNo+";");
	    out.println("int ctr[colNo];");
	    out.println("for(int i=0;i<colNo;i++) ctr[i]=0;");
	    out.println("");
	    out.println("");
	    out.println("");


	    out.println("for(int b=0;b<colNo;b++){");
	    out.println("int validsample=0;");
	    out.println("for(int i=0;i<sampleNo;i++){");

	    for(int i=0; i<freeParIDs.length; i++){
		out.println(pars[freeParIDs[i]].name+"=lb["+i+
			    "]+dsfmt_genrand_close_open(&dsfmt)*(ub["+i+"]-lb["+i+"]);");
	    }

// 	    for(int i=0; i<freeParIDs.length; i++){
// 		out.println(pars[freeParIDs[i]].name+"="+pars[freeParIDs[i]].lowerBound+
// 			    "+dsfmt_genrand_close_open(&dsfmt)*"+
// 			    (pars[freeParIDs[i]].upperBound-pars[freeParIDs[i]].lowerBound)+";");
// 	    }
	    


	    /*
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0)

		    out.println(vars[i].name+"="+vars[i].name+"i["+vars[i].name+
				"init]+dsfmt_genrand_close_open(&dsfmt)*("+
				vars[i].name+"i["+vars[i].name+"init+1]-"+
				vars[i].name+"i["+vars[i].name+"init]);");

		    //out.println(vars[i].name+"="+vars[i].init+";");
		else
		    out.println(vars[i].name+"=f"+vars[i].name+"(0);");
	    }
	    out.println("");
	    */
	    
	    out.println("(*P[ktest])=lb[ktest]+(b+0.5)*(ub[ktest]-lb[ktest])/(colNo*1.0);");

	    //out.println(pars[ktest].name+"="+pars[ktest].lowerBound+"+(b+0.5)*"+(pars[ktest].upperBound-pars[ktest].lowerBound)+"/(colNo*1.0);");
	    out.println("");
	    out.println("");


	    // ODE solver
	    out.println("real ropt[OPT_SIZE], reltol, t, tout, ttmp;");
	    out.println("long int iopt[OPT_SIZE];");
	    out.println("N_Vector y;");
	    out.println("real abstol;");
	    out.println("void *cvode_mem;");
	    out.println("int iout, flag, fail=0;");
	    out.println("int i,k,j,l;");
	    out.println("int t1;");
	    out.println("double value = 0.0;");
	    out.println("double sum,err,max,wmax;");
	    out.println("double w[xn],m[xn];");
	    out.println("");

	    out.println("");
	    out.println("y=N_VNew(NEQ, NULL);");
	    for(int j=1;j<varNum;j++) 
		out.println("Ith(y,"+j+")="+vars[j].name+";");
	    out.println("");
	    out.println("reltol = RTOL;");
	    out.println("abstol = ATOL;");
	   
	    out.println("");
	    out.println("cvode_mem = CVodeMalloc(NEQ, difeq, T0, y, BDF, NEWTON, SS, &reltol, &abstol, NULL, NULL, FALSE, iopt, ropt, NULL);");
	    out.println("");
 
	    out.println("if(cvode_mem == NULL){");
	    out.println("printf(\"CVodeMalloc failed.\\n\");");
	    out.println("exit(1);");
	    out.println("}");
	    out.println("CVDense(cvode_mem, NULL, NULL);");
	    out.println("");
	    out.println("ttmp=100;");
	    out.println("for(iout=0;iout<tn;iout++){");
	    out.println("tout=T[iout];");
	    out.println("for(;ttmp<=tout;ttmp+=100){");
	    out.println("flag = CVode(cvode_mem, ttmp, y, &t, NORMAL);");
	    out.println("if (flag != SUCCESS){");
	    out.println("printf(\"CVode failed, flag=%d.\\n\", flag);");
	    out.println("fail=1;break;");
	    out.println("}");
	    out.println("}");
	    out.println("if(fail==1) break;");
	    for(int j=0;j<dataVarIDs.length;j++)
		out.println("result["+j+"][iout] = Ith(y,"+dataVarIDs[j]+");");
	    out.println("ttmp=tout+100;");
	    out.println("}");

	    out.println("N_VFree(y);");
	    out.println("CVodeFree(cvode_mem);");
	    
	    out.println("if(fail==0){");
	    out.println("validsample++;");
	    out.println("// evaluate");
	    out.println("t1 = tn;");
	    //out.println("for(i=0;i<NEPR;i++){");
	    //out.println("nepr=i;");
	    out.println("for(l=0,wmax=-1;l<xn;l++){");	    
	    out.println("for(k=0,sum=0.0,max=-1; k<t1; k++){");
	    //out.println("sum+=DataEpr[nepr][l][k]*DataEpr[nepr][l][k];");
	    out.println("sum+=data[l][k];"); // [lb]100124 mean weight
	    if(normalize==1){ // normalization
		out.println("if(result[l][k]>max) max=result[l][k];");
	    }
	    out.println("}");
	    out.println("w[l]=t1/sum;");
	    out.println("if(w[l]>wmax) wmax=w[l];");
	    if(normalize==1){ // normalization
	    out.println("m[l]=max;");
	    }
	    out.println("}");
	    out.println("");
	    out.println("for(l=0;l<xn;l++){");
	    out.println("for(k=0,sum=0.0; k<t1; k++){");
	    if(normalize==0){ // no normalization
		out.println("err=data[l][k]-result[l][k];");
	    } else if(normalize==1){ // normalization with max
		out.println("err=data[l][k]-result[l][k]/m[l];");
	    } 
	    out.println("sum +=err*err;");
	    out.println("}");
	    out.println("value+=sum*w[l]/wmax;");
	    out.println("}");
	    out.println("");

	    out.println("if(value<threshold){");
	    out.println("ctr[b]++;");
	    out.println("}");
	    out.println("}");
	    out.println("}");
	    out.println("freq=ctr[b]/(validsample*1.0);");
	    out.println("fprintf(out,\"%f\\t\",freq);");
	    out.println("}");
	    out.println("fprintf(out,\"\\n\");");
	    out.println("");
	    
	    out.println("fclose(out);");
	    out.println("return 0;");
	    
	    out.println("}");

	    // functions
	    out.println("/* differential equations */");
	    out.println("static void difeq(integer N, real t, N_Vector y, N_Vector ydot, void *f_data){");
	    for(int i=1;i<varNum;i++) 
		out.println("real "+vars[i].name+", d"+vars[i].name+";");
	    for(int i=1;i<varNum;i++) 
		out.println(vars[i].name+" = Ith(y,"+i+");");
	    for(int i=1;i<varNum;i++)
		out.println("d"+vars[i].eqn+";");
	    for(int i=1;i<varNum;i++) 
		out.println("Ith(ydot,"+i+")=d"+vars[i].name+";");
	    out.println("");
	    out.println("return;");
	    out.println("}");
	    out.println("");




	    
	    out.flush();
	    
	    
	} catch(IOException e){
	    e.printStackTrace();
	}
    }


 
    //------------- point-based ODE simulation (Runge-Kutta) --------------


    public void genTESTcode(int timeptNum,double dt,double simTime,String outputDir){
	try{
	    FileOutputStream outfile=new FileOutputStream(modelName+"Test.c");
	    PrintWriter out=new PrintWriter(outfile);
	    
	    // header
	    out.println("#include <stdio.h>");
	    out.println("#include <math.h>");
	    out.println("");
	    
	    // global variable
	    for(int i=1;i<varNum;i++) 
		out.println("double "+vars[i].name+","+vars[i].name+"p;");
	    for(int i=1;i<parNum;i++)
		out.println("double "+pars[i].name+";");
	    out.println("double result["+(timeptNum+1)+"]["+varNum+"];");
	    
	    out.println("");
	    
	    // functions
	    for(int i=1;i<varNum;i++){
		StringTokenizer st = new StringTokenizer(vars[i].eqn, "=");
		String name=st.nextToken(); //name
		String body=st.nextToken(); // body
		out.println("double f"+name+"(double "+name+"){");
		out.println("return "+body+";");
		out.println("}");	
	    }
	    out.println("");

	    // main
	    out.println("int main(int argc, char *argv[]){");
	    out.println("double dt="+dt+",halfdt=dt/2.0;");
	    out.println("int tps=(int)("+simTime+"/dt),t,i;");
	    out.println("int block=tps/"+timeptNum+",tb;");
	    out.println("double halfF1,halfF2,F3,F4;");

	    out.println("");
	    for(int i=1;i<parNum;i++){
		out.println(pars[i].name+"="+pars[i].init+";");
	    }
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0)
		    out.println(vars[i].name+"="+vars[i].init+";");
		else
		    out.println(vars[i].name+"=f"+vars[i].name+"(0);");
	    }
	    out.println("");

	    // ODE solver
	    out.println("for(int t=1;t<=tps+1;t++){");
	    out.println("if(t%block==1){ // collect results");
	    out.println("tb=t/block;");
	    for(int i=1;i<varNum;i++){ //collect results
		out.println("result[tb]["+i+"]="+vars[i].name+";");
	    }
	    out.println("}");

	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0){
		    out.println("// "+vars[i].name);
		    out.println("halfF1=halfdt*f"+vars[i].name+"("+vars[i].name+");");
		    out.println("halfF2=halfdt*f"+vars[i].name+"("+vars[i].name+"+halfF1);");
		    out.println("F3=dt*f"+vars[i].name+"("+vars[i].name+"+halfF2);");
		    out.println("F4=dt*f"+vars[i].name+"("+vars[i].name+"+F3);");
		    out.println(vars[i].name+"p="+vars[i].name+"+(2*halfF1+4*halfF2+2*F3+F4)/6.0;");
		}
	    }
	    
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0)
		    out.println(vars[i].name+"="+vars[i].name+"p;");
		else
		    out.println(vars[i].name+"=f"+vars[i].name+"(0);");
	    }
	    out.println("}");
	    
	    out.println("");
	    
	    // output
	    out.println("");
	    out.println("// output");
	    out.println("FILE *out;");
	    out.println("char buffer[256];");
	    out.println("snprintf(buffer, sizeof(buffer), \"dummy.txt\");");
	    out.println("int idx=0;");
	    out.println("");
	    
	    out.println("snprintf(buffer, sizeof(buffer), \""+outputDir+"/"+modelName+"Test.txt\");");
	    out.println("out=fopen(buffer, \"w\");");
	    out.println("fprintf(out,\"Time,\");"); // header
	    out.println("for(int i=1;i<"+varNum+";i++)");
	    out.println("fprintf(out,\"x%d,\",i);");
	    out.println("fprintf(out,\"\\n\");");
	    out.println("for(int tb=0;tb<"+(timeptNum+1)+";tb++){");
	    out.println("fprintf(out,\"%d,\",tb);");
	    out.println("for(int i=1;i<"+varNum+";i++)");
	    out.println("fprintf(out,\"%E,\",result[tb][i]);");
	    out.println("fprintf(out,\"\\n\");");
	    out.println("}");
	    out.println("fclose(out);");

	    out.println("return 0;");
	    
	    out.println("}");
	    
	    out.flush();
	    
	    
	} catch(IOException e){
	    e.printStackTrace();
	}
    }

 
    //------------- point-based ODE simulation (CVODE) --------------
    // todo

    public void genTESTCcode(int timeptNum,double dt,double simTime,String outputDir){
	try{
	    FileOutputStream outfile=new FileOutputStream(modelName+"Test.c");
	    PrintWriter out=new PrintWriter(outfile);
	    
	    // header
	    out.println("#include <stdio.h>");
	    out.println("#include <math.h>");
	    out.println("");
	    
	    // global variable
	    for(int i=1;i<varNum;i++) 
		out.println("double "+vars[i].name+","+vars[i].name+"p;");
	    for(int i=1;i<parNum;i++)
		out.println("double "+pars[i].name+";");
	    out.println("double result["+(timeptNum+1)+"]["+varNum+"];");
	    
	    out.println("");
	    
	    // functions
	    for(int i=1;i<varNum;i++){
		StringTokenizer st = new StringTokenizer(vars[i].eqn, "=");
		String name=st.nextToken(); //name
		String body=st.nextToken(); // body
		out.println("double f"+name+"(double "+name+"){");
		out.println("return "+body+";");
		out.println("}");	
	    }
	    out.println("");

	    // main
	    out.println("int main(int argc, char *argv[]){");
	    out.println("double dt="+dt+",halfdt=dt/2.0;");
	    out.println("int tps=(int)("+simTime+"/dt),t,i;");
	    out.println("int block=tps/"+timeptNum+",tb;");
	    out.println("double halfF1,halfF2,F3,F4;");

	    out.println("");
	    for(int i=1;i<parNum;i++){
		out.println(pars[i].name+"="+pars[i].init+";");
	    }
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0)
		    out.println(vars[i].name+"="+vars[i].init+";");
		else
		    out.println(vars[i].name+"=f"+vars[i].name+"(0);");
	    }
	    out.println("");

	    // ODE solver
	    out.println("for(int t=1;t<=tps+1;t++){");
	    out.println("if(t%block==1){ // collect results");
	    out.println("tb=t/block;");
	    for(int i=1;i<varNum;i++){ //collect results
		out.println("result[tb]["+i+"]="+vars[i].name+";");
	    }
	    out.println("}");

	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0){
		    out.println("// "+vars[i].name);
		    out.println("halfF1=halfdt*f"+vars[i].name+"("+vars[i].name+");");
		    out.println("halfF2=halfdt*f"+vars[i].name+"("+vars[i].name+"+halfF1);");
		    out.println("F3=dt*f"+vars[i].name+"("+vars[i].name+"+halfF2);");
		    out.println("F4=dt*f"+vars[i].name+"("+vars[i].name+"+F3);");
		    out.println(vars[i].name+"p="+vars[i].name+"+(2*halfF1+4*halfF2+2*F3+F4)/6.0;");
		}
	    }
	    
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0)
		    out.println(vars[i].name+"="+vars[i].name+"p;");
		else
		    out.println(vars[i].name+"=f"+vars[i].name+"(0);");
	    }
	    out.println("}");
	    
	    out.println("");
	    
	    // output
	    out.println("");
	    out.println("// output");
	    out.println("FILE *out;");
	    out.println("char buffer[256];");
	    out.println("snprintf(buffer, sizeof(buffer), \"dummy.txt\");");
	    out.println("int idx=0;");
	    out.println("");
	    
	    out.println("snprintf(buffer, sizeof(buffer), \""+outputDir+"/"+modelName+"Test.txt\");");
	    out.println("out=fopen(buffer, \"w\");");
	    out.println("fprintf(out,\"Time,\");"); // header
	    out.println("for(int i=1;i<"+varNum+";i++)");
	    out.println("fprintf(out,\"x%d,\",i);");
	    out.println("fprintf(out,\"\\n\");");
	    out.println("for(int tb=0;tb<"+(timeptNum+1)+";tb++){");
	    out.println("fprintf(out,\"%d,\",tb);");
	    out.println("for(int i=1;i<"+varNum+";i++)");
	    out.println("fprintf(out,\"%E,\",result[tb][i]);");
	    out.println("fprintf(out,\"\\n\");");
	    out.println("}");
	    out.println("fclose(out);");

	    out.println("return 0;");
	    
	    out.println("}");
	    
	    out.flush();
	    
	    
	} catch(IOException e){
	    e.printStackTrace();
	}
    }



    //------------- point-based Parameter estiamtion (CVODE solver) -------------- 
    // option 0: single node; option 1: MPI
    public void genSRESCcode(int[]freeParIDs, int[]freeVarIDs, int[]dataTPs,int[]dataVarIDs, String dataFile,int[]condVarIDs, int[]condParIDs, double[][]conds, int normalize, int option){
	try{
	    FileOutputStream outfile=new FileOutputStream(modelName+"SRESc.C");
	    PrintWriter out=new PrintWriter(outfile);
	    // header
	    out.println("#include <stdio.h>");
	    out.println("#include <stdlib.h>");
	    out.println("#include <math.h>");
	    out.println("");
	    out.println("/* header files for SRES */");
	    out.println("#include \"sharefunc.h\"");
	    out.println("#include \"ESSRSort.h\"");
	    out.println("#include \"ESES.h\"");
	    out.println("");
	    out.println("/* header files and macros for CVODE */");
	    out.println("#include \"llnltyps.h\"");
	    out.println("#include \"cvode.h\"");
	    out.println("#include \"cvdense.h\"");
	    out.println("#include \"nvector.h\"");
	    out.println("#include \"dense.h\"");
	    out.println("#define Ith(v,i) N_VIth(v,i-1)");
	    out.println("#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)");
	    out.println("");

	    // global variable
	    out.println("int NEQ;");
	    out.println("double RTOL, ATOL;");
	    out.println("double T0, T1, Tm;");
	    out.println("int SEP;");
	    out.println("int NEPR;");
	    out.println("static int nepr;");
	    out.println("double ***DataEpr, ***DataPrd;");
	    out.println("double *DataTPs;");
	    out.println("int *DataIDs;");
	    
	    for(int i=1;i<varNum;i++) 
		out.println("double "+vars[i].name+";");
	    for(int i=1;i<parNum;i++)
		out.println("double "+pars[i].name+";");
	    
	    for(int i=0;i<condVarIDs.length;i++)
		out.println("double *c"+vars[condVarIDs[i]].name+";");
	    
	    out.println("int tn,xn;");
	    out.println("");
	    out.println("void fitness(double *, double *, double *);");
	    out.println("ESfcnTrsfm *trsfm;");
	    out.println("double transform(double);");
	    out.println("void funReadEpr(void);");
	    out.println("static void difeq(integer N, real t, N_Vector y, N_Vector ydot, void *f_data);");
	    out.println("static void PrintFinalStats(long int iopt[]);");

	    // main
	    out.println("int main(int argc, char ** argv){");
	    out.println("FILE *out,*input2,*out2;");
	    out.println("char buffer[64];");
	    out.println("int i;");
	    out.println("ESParameter *param;");
	    out.println("ESPopulation *population;");
	    out.println("ESStatistics *stats;");
	    out.println("unsigned int seed;");
	    out.println("int es;");
	    out.println("int constraint, dim;");
	    out.println("double *ub, *lb;");
	    out.println("int miu, lambda, gen;");
	    out.println("double gamma, alpha, varphi;");
	    out.println("int retry;");
	    out.println("double pf;");
	    out.println("seed = shareDefSeed;");
	    out.println("gamma = esDefGamma;");
	    out.println("alpha = esDefAlpha;");
	    out.println("varphi = esDefVarphi;");
	    out.println("retry = esDefRetry;");
	    out.println("pf = essrDefPf;");
	    out.println("es = esDefESSlash;");
	    out.println("constraint = 1;"); // constraints??
	    out.println("dim = "+(freeParIDs.length+freeVarIDs.length)+";");
	    out.println("");
	    out.println("input2=fopen(\"input2.dat\", \"r\");");
	    out.println("fgets(buffer, sizeof(buffer), input2);");
	    out.println("miu=strtol(buffer,NULL,0);");
	    out.println("fgets(buffer, sizeof(buffer), input2);");
	    out.println("lambda=strtol(buffer,NULL,0);");
	    out.println("fgets(buffer, sizeof(buffer), input2);");
	    out.println("gen=strtol(buffer,NULL,0);");
	    out.println("");
	    out.println("ub = NULL;");
	    out.println("lb = NULL;");
	    out.println("ub = ShareMallocM1d(dim);");
	    out.println("lb = ShareMallocM1d(dim);");
	    for(int i=0; i<freeParIDs.length; i++)
		out.println("lb["+i+"]="+pars[freeParIDs[i]].lowerBound
			    +"; ub["+i+"]="+pars[freeParIDs[i]].upperBound+";");
	    for(int i=0; i<freeVarIDs.length; i++)
		out.println("lb["+(i+freeParIDs.length)+"]="+vars[freeVarIDs[i]].lowerBound
			    +"; ub["+(i+freeParIDs.length)+"]="+vars[freeVarIDs[i]].upperBound+";");
	    
	    out.println("trsfm = (ESfcnTrsfm *)ShareMallocM1c(dim*sizeof(ESfcnTrsfm));");
	    out.println("for(i=0;i<dim;i++) trsfm[i] = transform;");
	    out.println("");
	    out.println("/* set ODE parameters */");
	    out.println("NEQ = "+(varNum-1)+";");
	    out.println("RTOL = 1e-6;");
	    out.println("ATOL = 1e-6;");
	    out.println("T0 = 0.0;");
	    //out.println("T1 = "+dt+";");
	    out.println("Tm = "+dataTPs[dataTPs.length-1]+";");
	    //int NEPR=1;
	    //for(int i=0;i<conds.length;i++)
	    //NEPR*=conds[i].length;
	    //out.println("NEPR = "+NEPR+";");
	    out.println("NEPR = "+conds.length+";");
	    out.println("");

	    // init
	    for(int i=1;i<varNum;i++) 
		out.println(vars[i].name+"="+vars[i].init+";");
	    out.println("");
	    for(int i=1;i<parNum;i++) 
		out.println(pars[i].name+"="+pars[i].init+";");
	    out.println("");

	    out.println("");
	    out.println("tn="+dataTPs.length+";");
	    out.println("xn="+dataVarIDs.length+";");
	    out.println("");
	    out.println("DataEpr = ShareMallocM3d(NEPR,NEQ,tn);");
	    out.println("DataPrd = ShareMallocM3d(NEPR,NEQ,tn);");
	    out.println("DataTPs = ShareMallocM1d(tn);");

	    for(int i=0; i<dataTPs.length; i++)
		out.print("DataTPs["+i+"]="+dataTPs[i]+";");
	    out.println();

	    out.println();
	    out.println("");
	    out.println("/* read experimental data */");
	    out.println("funReadEpr();");
	    out.println("");

	    if(option==0){
		out.println("ESInitial(seed, &param, trsfm, fitness, es, constraint, dim, ub, lb, miu, lambda, gen, gamma, alpha, varphi, retry, &population, &stats);");
	    } else if(option==1){
		out.println("ESInitial(&argc,&argv,seed, &param, trsfm, fitness, es, constraint, dim, ub, lb, miu, lambda, gen, gamma, alpha, varphi, retry, &population, &stats);");
	    }

	    //out.println("printf(\"\\n======\\nseed=%u\\n======\\n\", param->seed);");
	    //out.println("");
	    //out.println("while(stats->curgen < param->gen)");
	    //out.println("ESStep(population, param, stats, pf);");
	    out.println("");
	    out.println("out=fopen(\"output1.dat\", \"w\");");
	    out.println("out2=fopen(\"output2.dat\", \"w\");");
	    out.println("");
	    out.println("while(stats->curgen < param->gen){");
	    out.println("ESStep(population, param, stats, pf);");
	    if(option==0){
		out.println("if(stats->curgen % 10 ==0)");
		out.println("fprintf(out2,\"gen=%d,dt=%d,bestgen=%d,bestfitness=%f,phi=%f,\\n\",stats->curgen,stats->dt,stats->bestgen,stats->bestindvdl->f,stats->bestindvdl->phi);");
	    }
	    out.println("}");
	    out.println("fprintf(out2,\"gen=%d,dt=%d,bestgen=%d,bestfitness=%f,phi=%f,\\n\",stats->curgen,stats->dt,stats->bestgen,stats->bestindvdl->f,stats->bestindvdl->phi);");
	    out.println("for(i=0;i<dim;i++)");
	    out.println("fprintf(out,\"%E,%E\\n\",stats->bestindvdl->op[i],stats->bestindvdl->sp[i]);");
	    out.println("fclose(out);");
	    out.println("fclose(out2);");
	    out.println("");
	    out.println("");


	    out.println("ESDeInitial(param, population, stats);");
	    out.println("");
	    out.println("ShareFreeM1d(ub);");
	    out.println("ub = NULL;");
	    out.println("ShareFreeM1d(lb);");
	    out.println("lb = NULL;");
	    out.println("ShareFreeM1c((char*)trsfm);");
	    out.println("trsfm = NULL;");

	    out.println("ShareFreeM3d(DataEpr, NEPR, NEQ);");
	    out.println("ShareFreeM3d(DataPrd, NEPR, NEQ);");
	    out.println("ShareFreeM1d(DataTPs);");
	    out.println("DataTPs = NULL;");
	    out.println("");
	    out.println("return 0;");
	    out.println("}");
	    out.println("");
	    out.println("");
	    // fitness
	    out.println("/* fitness fucntion */");
	    out.println("void fitness(double *x, double *f, double *g){");
	    out.println("real ropt[OPT_SIZE], reltol, t, tout, ttmp;");
	    out.println("long int iopt[OPT_SIZE];");
	    out.println("N_Vector y;");
	    out.println("real abstol;");
	    out.println("void *cvode_mem;");
	    out.println("int iout, flag;");
	    out.println("int i,k,j,l,i0,i1,i2,i3,i4,i5;");
	    out.println("int t1;");
	    out.println("double value = 0.0;");
	    out.println("double sum,err,max,wmax;");
	    out.println("double w[xn],m[xn];");
	    out.println("");
	    for(int i=0; i<freeParIDs.length; i++)
		out.println(pars[freeParIDs[i]].name+"=(trsfm["+i+"])(x["+i+"]);");
	    for(int i=0; i<freeVarIDs.length; i++)
		out.println(vars[freeVarIDs[i]].name+"=(trsfm["+(i+freeParIDs.length)+"])(x["+(i+freeParIDs.length)+"]);");
	    out.println("");

	    for(int i=0; i<conds.length;i++){
		out.println("nepr="+i+";");
		for(int j=0;j<condVarIDs.length;j++) 
		    out.println(vars[condVarIDs[j]].name+"="+conds[i][j]+";");
		for(int j=0;j<condParIDs.length;j++)
		    out.println(pars[condParIDs[j]].name+"="+conds[i][j+condVarIDs.length]+";");
		
		out.println("");
		out.println("y=N_VNew(NEQ, NULL);");
		for(int j=1;j<varNum;j++) 
		    out.println("Ith(y,"+j+")="+vars[j].name+";");
		out.println("");
		out.println("reltol = RTOL;");
		out.println("abstol = ATOL;");
		out.println("");
		out.println("cvode_mem = CVodeMalloc(NEQ, difeq, T0, y, BDF, NEWTON, SS, &reltol, &abstol, NULL, NULL, FALSE, iopt, ropt, NULL);");
		out.println("");
		out.println("if(cvode_mem == NULL){");
		//out.println("printf(\"CVodeMalloc failed.\\n\");");
		out.println("exit(1);");
		out.println("}");
		out.println("CVDense(cvode_mem, NULL, NULL);");
		out.println("");
		out.println("ttmp=100;");
		out.println("for(iout=0;iout<tn;iout++){");
		out.println("tout=DataTPs[iout];");
		out.println("for(;ttmp<=tout;ttmp+=100){");
		out.println("flag = CVode(cvode_mem, ttmp, y, &t, NORMAL);");
		out.println("if (flag != SUCCESS){");
		//out.println("printf(\"CVode failed, flag=%d.\\n\", flag);");
		out.println("(*f) = 80000;");
		out.println("g[0] = 0.0;");
		out.println("return;");
		out.println("}");
		out.println("}");
		for(int j=0;j<dataVarIDs.length;j++)
		    out.println("DataPrd[nepr]["+j+"][iout] = Ith(y,"+dataVarIDs[j]+");");
		out.println("ttmp=tout+100;");
		out.println("}");
		out.println("N_VFree(y);");
		out.println("CVodeFree(cvode_mem);");
		out.println("//PrintFinalStats(iopt);");
	    }
	    out.println("");
	    out.println("t1 = tn;");
	    out.println("for(i=0;i<NEPR;i++){");
	    out.println("nepr=i;");
	    out.println("for(l=0,wmax=-1;l<xn;l++){");	    
	    out.println("for(k=0,sum=0.0,max=-1; k<t1; k++){");
	    out.println("sum+=DataEpr[nepr][l][k]*DataEpr[nepr][l][k];");
	    if(normalize==1){ // normalization
		out.println("if(DataPrd[nepr][l][k]>max) max=DataPrd[nepr][l][k];");
	    }
	    out.println("}");
	    out.println("w[l]=t1/sum;");
	    out.println("if(w[l]>wmax) wmax=w[l];");
	    if(normalize==1){ // normalization
	    out.println("m[l]=max;");
	    }
	    out.println("}");
	    out.println("");
	    out.println("for(l=0;l<xn;l++){");
	    out.println("for(k=0,sum=0.0; k<t1; k++){");
	    if(normalize==0){ // no normalization
		out.println("err=DataEpr[nepr][l][k]-DataPrd[nepr][l][k];");
	    } else if(normalize==1){ // normalization with max
		out.println("err=DataEpr[nepr][l][k]-DataPrd[nepr][l][k]/m[l];");
	    } 
	    out.println("sum +=err*err;");
	    out.println("}");
	    out.println("value+=sum*w[l]/wmax;");
	    out.println("}");
	    out.println("}");
	    out.println("");
	    out.println("(*f) = value;");
	    out.println("g[0] = 0.0;");
	    out.println("");
	    out.println("return;");
	    out.println("}");
	    out.println("");
	    out.println("");
	    out.println("");


	    out.println("void funReadEpr(void){");
	    out.println("char file[] = \""+dataFile+"\";"); // <--
	    out.println("char buf[shareDefMaxLine];");
	    out.println("char **s1;");
	    out.println("FILE *fp;");
	    out.println("int i,j,k;");
	    out.println("int t,n;");
	    out.println("");
	    out.println("if((fp = fopen(file, \"r\")) == NULL){");
	    out.println("printf(\"fopen %s failed!\\n\", file);");
	    out.println("exit(1);");
	    out.println("}");
	    out.println("");
	    out.println("t = tn;");
	    out.println("for(i=0; i<NEPR; i++)");
	    out.println("for(j=0; j<t; j++){");
	    out.println("fgets(buf, shareDefMaxLine, fp);");
	    out.println("ShareChop(buf);");
	    out.println("s1 = ShareSplitStr(buf, \"\\t\", &n, shareDefNullNo);");
	    out.println("for(k=1;k<n;k++)");
	    out.println("DataEpr[i][k-1][j] = atof(s1[k]);");
	    out.println("ShareFreeM2c(s1, n);");
	    out.println("}");
	    out.println("");
	    out.println("fclose(fp);");
	    out.println("");
	    out.println("return;");
	    out.println("}");
	    out.println("");
	    out.println("/* differential equations */");
	    out.println("static void difeq(integer N, real t, N_Vector y, N_Vector ydot, void *f_data){");
	    for(int i=1;i<varNum;i++) 
		out.println("real "+vars[i].name+", d"+vars[i].name+";");
	    for(int i=1;i<varNum;i++) 
		out.println(vars[i].name+" = Ith(y,"+i+");");
	    for(int i=1;i<varNum;i++)
		out.println("d"+vars[i].eqn+";");
	    for(int i=1;i<varNum;i++) 
		out.println("Ith(ydot,"+i+")=d"+vars[i].name+";");
	    out.println("");
	    out.println("return;");
	    out.println("}");
	    out.println("");
	    out.println("static void PrintFinalStats(long int iopt[]){");
	    out.println("}");
	    out.println("");
	    out.println("double transform(double x){");
	    out.println("double y;");
	    out.println("y = x;");
	    out.println("return y;");
	    out.println("}");
	    out.flush();

	} catch(IOException e){
	    e.printStackTrace();
	}
    }
    

   //------------- point-based Parameter estiamtion (Runge-Kutta solver) --------------
    
    public void genSRESRcode(int timeptNum,double dt,double simTime,String outputDir){
	try{
	    FileOutputStream outfile=new FileOutputStream(modelName+"Test.c");
	    PrintWriter out=new PrintWriter(outfile);
	    
	    // header
	    out.println("#include <stdio.h>");
	    out.println("#include <math.h>");
	    out.println("");
	    
	    // global variable
	    for(int i=1;i<varNum;i++) 
		out.println("double "+vars[i].name+","+vars[i].name+"p;");
	    for(int i=1;i<parNum;i++)
		out.println("double "+pars[i].name+";");
	    out.println("double result["+(timeptNum+1)+"]["+varNum+"];");
	    
	    out.println("");
	    
	    // functions
	    for(int i=1;i<varNum;i++){
		StringTokenizer st = new StringTokenizer(vars[i].eqn, "=");
		String name=st.nextToken(); //name
		String body=st.nextToken(); // body
		out.println("double f"+name+"(double "+name+"){");
		out.println("return "+body+";");
		out.println("}");	
	    }
	    out.println("");

	    // main
	    out.println("int main(int argc, char *argv[]){");
	    out.println("double dt="+dt+",halfdt=dt/2.0;");
	    out.println("int tps=(int)("+simTime+"/dt),t,i;");
	    out.println("int block=tps/"+timeptNum+",tb;");
	    out.println("double halfF1,halfF2,F3,F4;");

	    out.println("");
	    for(int i=1;i<parNum;i++){
		out.println(pars[i].name+"="+pars[i].init+";");
	    }
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0)
		    out.println(vars[i].name+"="+vars[i].init+";");
		else
		    out.println(vars[i].name+"=f"+vars[i].name+"(0);");
	    }
	    out.println("");

	    // ODE solver
	    out.println("for(int t=1;t<=tps+1;t++){");
	    out.println("if(t%block==1){ // collect results");
	    out.println("tb=t/block;");
	    for(int i=1;i<varNum;i++){ //collect results
		out.println("result[tb]["+i+"]="+vars[i].name+";");
	    }
	    out.println("}");

	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0){
		    out.println("// "+vars[i].name);
		    out.println("halfF1=halfdt*f"+vars[i].name+"("+vars[i].name+");");
		    out.println("halfF2=halfdt*f"+vars[i].name+"("+vars[i].name+"+halfF1);");
		    out.println("F3=dt*f"+vars[i].name+"("+vars[i].name+"+halfF2);");
		    out.println("F4=dt*f"+vars[i].name+"("+vars[i].name+"+F3);");
		    out.println(vars[i].name+"p="+vars[i].name+"+(2*halfF1+4*halfF2+2*F3+F4)/6.0;");
		}
	    }
	    
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0)
		    out.println(vars[i].name+"="+vars[i].name+"p;");
		else
		    out.println(vars[i].name+"=f"+vars[i].name+"(0);");
	    }
	    out.println("}");
	    
	    out.println("");
	    
	    // output
	    out.println("");
	    out.println("// output");
	    out.println("FILE *out;");
	    out.println("char buffer[256];");
	    out.println("snprintf(buffer, sizeof(buffer), \"dummy.txt\");");
	    out.println("int idx=0;");
	    out.println("");
	    
	    out.println("snprintf(buffer, sizeof(buffer), \""+outputDir+"/"+modelName+"Test.txt\");");
	    out.println("out=fopen(buffer, \"w\");");
	    out.println("fprintf(out,\"Time,\");"); // header
	    out.println("for(int i=1;i<"+varNum+";i++)");
	    out.println("fprintf(out,\"x%d,\",i);");
	    out.println("fprintf(out,\"\\n\");");
	    out.println("for(int tb=0;tb<"+(timeptNum+1)+";tb++){");
	    out.println("fprintf(out,\"%d,\",tb);");
	    out.println("for(int i=1;i<"+varNum+";i++)");
	    out.println("fprintf(out,\"%f,\",result[tb][i]);");
	    out.println("fprintf(out,\"\\n\");");
	    out.println("}");
	    out.println("fclose(out);");

	    out.println("return 0;");
	    
	    out.println("}");
	    
	    out.flush();
	    
	    
	} catch(IOException e){
	    e.printStackTrace();
	}
    }





    //------------- probabilistic simulation (Factored-Frotier inference)-----------

    //------------- generating simulation code (c) -----------

    public void genFFcodeC(int timeptNum, int simTime, String outputDir, int option, int[]freeParIDs, int[]dataTPs,int[]dataVarIDs, double[]weights, String dataFile, boolean fulldata){
	try{
	    FileOutputStream outfile=new FileOutputStream(modelName+"_ff"+option+".c");
	    PrintWriter out=new PrintWriter(outfile);

	    // header
	    out.println("#include <stdio.h>");
	    out.println("#include <stdlib.h>");
	    out.println("#include <string.h>");
	    out.println("");
	    out.println("void load"+option+"(void);");
	    out.println("void simulate(int kvar[]);");
	    out.println("void output(void);");
	    out.println("void outputMarginal(void);");
	    out.println("double eval(int kvar[]);");
	    out.println("");
	    out.println("int tps="+timeptNum+",tb,MAX_LINE=1024;");
	    out.println("double dt="+(simTime/(timeptNum*1.0))+";");

	    // ----- gloabal variables
	    int binMax=0;
	    for(int i=1;i<varNum;i++)		
		if(vars[i].binNum>binMax)
		    binMax=vars[i].binNum;

	    // mean of each bin 
	    double [][]binM=new double[varNum][];
	    out.print("double M[]["+binMax+"]={{"); // updated 072109
	    for(int i=0;i<binMax;i++) out.print("-1,");
	    out.println("},");
	    for(int i=1;i<varNum;i++){
		int B=vars[i].binNum;
		binM[i]=new double[B];
		out.print("{");
		for(int j=0;j<B;j++){
		    binM[i][j]=((vars[i].binThresholds)[j+1]+(vars[i].binThresholds)[j])*0.5;
		    out.print(binM[i][j]+",");
		}
		for(int j=B;j<binMax;j++) out.print("-1,");
		out.println("},");
	    }
	    out.println("};");
	    out.println("");

	    // cpd
	    for(int i=1;i<varNum;i++){
		int n=(vars[i].varIDs).length+(vars[i].parIDs).length+1;
		out.print("double x"+i+"prob["+timeptNum+"]");
		int tmp_v=(vars[i].varIDs).length;
		int tmp_p=(vars[i].parIDs).length;
		for(int j=0;j<tmp_p;j++)
		    out.print("["+pars[(vars[i].parIDs)[j]].binNum+"]");
		for(int j=0;j<tmp_v;j++)
		    out.print("["+vars[(vars[i].varIDs)[j]].binNum+"]");
		out.print("["+vars[i].binNum+"]");
		out.println(";");
	    }

	    out.println("double result["+varNum+"]["+(1+timeptNum)+"];");
	    out.println("double x["+varNum+"]["+(1+timeptNum)+"]["+binMax+"];");

	    // parameter estimation
	    out.println("");
	    out.print("int P[]={");
	    for(int i=0;i<freeParIDs.length;i++)
		out.print(freeParIDs[i]+",");
	    out.println("};");

	    out.print("int T[]={");
	    for(int i=0;i<dataTPs.length;i++)
		out.print(dataTPs[i]+",");
	    out.println("};");
	    
	    out.print("int X[]={");
	    for(int i=0;i<dataVarIDs.length;i++)
		out.print(dataVarIDs[i]+",");
	    out.println("};");

	    out.print("double weights[]={");
	    for(int i=0;i<weights.length;i++)
		out.print(weights[i]+",");
	    out.println("};");

	    out.print("int k[]={-1,");
	    for(int i=1;i<parNum;i++)
		out.print(pars[i].binInit+",");	    
	    out.println("};");

	    out.print("int xbinNum[]={-1,");
	    for(int i=1;i<varNum;i++)
		out.print(vars[i].binNum+",");	    
	    out.println("};");

	    out.println("double expdata["+dataTPs.length+"]["+dataVarIDs.length+"];");

	    out.println("");
	    // -------------------


	    // ------- initialization
	    out.println("void load"+option+"(){");
	    //loadExp	    
	    out.println("");
	    out.println("//loadExp");	    

	    out.println("char buffer[MAX_LINE];");
	    out.println("char *cp;");
	    out.println("FILE *fin;");

	    out.println("");
	    out.println("double tmp;");
	    out.println("int i,j,k,ctr;");

	    out.println("fin = fopen(\""+dataFile+"\",\"r\");");
	    out.println("fgets(buffer, sizeof(buffer), fin); //header");
	    out.println("fgets(buffer, sizeof(buffer), fin); //init");
	    out.println("for(j=0;j<"+dataTPs.length+";j++){");
	    out.println("ctr=0;");
	    out.println("fgets(buffer, sizeof(buffer), fin);");
	    out.println("cp = (char *)strtok(buffer, \"\\t, \"); //time");
	    if(fulldata){
		out.println("for(k=1;k<"+varNum+";k++){");
		out.println("if(ctr=="+dataVarIDs.length+") break;");
		out.println("cp = (char *)strtok(NULL, \"\\t, \");");
		out.println("tmp=strtod(cp, NULL);");
		out.println("if(k==X[ctr]){");
		out.println("expdata[j][ctr]=tmp;");
		out.println("ctr++;");
		out.println("}");
		out.println("}");
	    } else {
		out.println("for(k=0;k<"+dataVarIDs.length+";k++){");
		out.println("cp = (char *)strtok(NULL, \"\\t, \");");
		out.println("tmp=strtod(cp, NULL);");
		out.println("expdata[j][ctr]=tmp;");
		out.println("}");
	    }
	    out.println("}");
	    out.println("fclose(fin);");


	    out.println("");
	    // load cpd
	    out.println("// load CPD");
	    out.println("int idx;double p;");
	    out.println("int vb,vi0,vi1,vi2,vi3,vi4,vi5,vi6,vi7,vi8,vi9,vi10;");
	    out.println("int pi0,pi1,pi2,pi3,pi4,pi5,pi6,pi7,pi8,pi9,pi10;");
	    out.println("for(tb=0;tb<"+timeptNum+";tb++){");
	    String tag="P";
	    if(modelName.equals("brown")) tag="XBNT"; // [patch] brown's model use XBNT lable
	    for(int i=1;i<varNum;i++){
		out.println("");
		out.println("snprintf(buffer, sizeof(buffer), \""+outputDir+"/tables/"+modelName+""+tag+"x"+i+"T%d.txt\", tb);"); 		
		out.println("fin = fopen(buffer, \"r\");");
		out.println("while (fgets(buffer, sizeof(buffer), fin)) {");
		out.println("cp = (char *)strtok(buffer, \"\\t, \");");
		out.println("idx=strtol(cp, NULL, 0);");
		out.println("cp = (char *)strtok(NULL, \"\\t, \");");
		out.println("p=strtod(cp, NULL);");
		int tmp_v=(vars[i].varIDs).length;
		int tmp_p=(vars[i].parIDs).length;
		out.println("vb=idx%"+vars[i].binNum+";");
		out.println("idx=idx/"+vars[i].binNum+";");	
		for(int j=tmp_v-1;j>=0;j--){
		    int tmp_b=vars[(vars[i].varIDs)[j]].binNum;
		    out.println("vi"+j+"=idx%"+tmp_b+";");
		    out.println("idx=idx/"+tmp_b+";");	    
		}
		for(int j=tmp_p-1;j>=0;j--){
		    int tmp_b=pars[(vars[i].parIDs)[j]].binNum;
		    out.println("pi"+j+"=idx%"+tmp_b+";");
		    out.println("idx=idx/"+tmp_b+";");	    
		}
		out.print("x"+i+"prob[tb]");
		for(int j=0;j<tmp_p;j++) out.print("[pi"+j+"]");
		for(int j=0;j<tmp_v;j++) out.print("[vi"+j+"]");
		out.println("[vb]=p;");
		out.println("}");
		out.println("");
		out.println("fclose(fin);");
	    }

	    
	    out.println("}");
	    out.println("}");
	    

	    // -------- simulation
	    out.println("void simulate(int kvar[]){");
	    out.println("int t,b,tb,i,i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15;");
	    out.println("double sum=0,margin=1,norm=0;");
	    out.println("");

	    out.println("for(i=0;i<"+freeParIDs.length+";i++) k[P[i]]=kvar[i];");

	    // init
	    out.println("");
	    out.println("//init");

	    /* bug fixed [lb] 122309 */
	    out.println("for(i=0;i<"+varNum+";i++)");
	    out.println("for(t=0;t<tps+1;t++)");
	    out.println("for(b=0;b<"+binMax+";b++)");
	    out.println("x[i][t][b]=0;");
	    out.println("");
	    out.println("for(i=0;i<"+varNum+";i++)");
	    out.println("for(t=0;t<tps+1;t++)");
	    out.println("result[i][t]=0;");
	    out.println("");
	    out.println("");

	    int N=1;
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0){
		    out.println("x["+i+"][0]["+vars[i].binInit+"]=1;");
		    N++;
		}
	    }
	    // expect
	    out.println("");
	    out.println("for(t=0;t<"+timeptNum+";t++){"); //<--
	    //out.println("tb=(int)(t/block);");
	    out.println("tb=t;");
	    out.println("// normalization");
	    out.println("for(i=1;i<"+N+";i++){");
	    out.println("norm=0; ");
	    out.println("for(b=0;b<xbinNum[i];b++) norm+=x[i][t][b]; if(norm==0) norm=1;");
	    out.println("for(b=0;b<xbinNum[i];b++) x[i][t][b]/=norm;");
	    out.println("}");
	    out.println("");
	    
	    // /* for reduced model !!!

	    for(int i=N;i<varNum;i++){
		out.println("// x"+i);
		out.println("for(b=0;b<"+vars[i].binNum+";b++){");
		//out.println("result["+i+"][t]+=x["+i+"][t][b]*b;");
		int tmp_p=(vars[i].parIDs).length;
		int tmp_v=(vars[i].varIDs).length;
		for(int j=0;j<tmp_v;j++) 
		    out.println("for(i"+j+"=0;i"+j+"<"+vars[(vars[i].varIDs)[j]].binNum+";i"+j+"++)");
      		out.print("x["+i+"][t][b]+=x"+i+"prob[tb]");
		for(int j=0;j<tmp_p;j++) out.print("[k["+(vars[i].parIDs)[j]+"]]");
		for(int j=0;j<tmp_v;j++) out.print("[i"+j+"]");
		out.print("[b]");
		for(int j=0;j<tmp_v;j++) 
		    out.print("*x["+(vars[i].varIDs)[j]+"][t][i"+j+"]");
		out.println(";");
		out.println("result["+i+"][t]+=x["+i+"][t][b]*M["+i+"][b];");
		out.println("}");
		
	    }

	    out.println("");	
	    out.println("// normalization2");
	    out.println("for(i="+N+";i<"+varNum+";i++){");
	    out.println("norm=0; ");
	    out.println("for(b=0;b<xbinNum[i];b++) norm+=x[i][t][b];if(norm==0) norm=1;");
	    out.println("for(b=0;b<xbinNum[i];b++) x[i][t][b]/=norm;");
	    out.println("}");
	    out.println("");	
	    // ----------- */

	    
	    for(int i=1;i<N;i++){
		out.println("// x"+i);
		out.println("for(b=0;b<"+vars[i].binNum+";b++){");
		out.println("result["+i+"][t]+=x["+i+"][t][b]*M["+i+"][b];");
		//out.println("result["+i+"][t]+=x["+i+"][t][b]*b;");
		int tmp_p=(vars[i].parIDs).length;
		int tmp_v=(vars[i].varIDs).length;
		for(int j=0;j<tmp_v;j++) 
		    out.println("for(i"+j+"=0;i"+j+"<"+vars[(vars[i].varIDs)[j]].binNum+";i"+j+"++)");
      		out.print("x["+i+"][t+1][b]+=x"+i+"prob[tb]");
		for(int j=0;j<tmp_p;j++) out.print("[k["+(vars[i].parIDs)[j]+"]]");
		for(int j=0;j<tmp_v;j++) out.print("[i"+j+"]");
		out.print("[b]");
		for(int j=0;j<tmp_v;j++) 
		    out.print("*x["+(vars[i].varIDs)[j]+"][t][i"+j+"]");
		out.println(";");
		out.println("}");
		
	    }
	    out.println("}");
	    out.println("}");
	    out.println("");
	    
	    //----------  output
	    out.println("// output");
	    out.println("void output(){");
	    out.println("FILE *out;");
	    out.println("char buffer[256];");
	    out.println("snprintf(buffer, sizeof(buffer), \"dummy.txt\");");
	    out.println("int idx=0,i,tb;");
	    out.println("");
	    
	    out.println("snprintf(buffer, sizeof(buffer), \""+modelName+"_T.csv\");");
	    out.println("out=fopen(buffer, \"w\");");
	    out.println("fprintf(out,\"Time,\");"); // header
	    out.println("for(i=1;i<"+varNum+";i++)");
	    out.println("fprintf(out,\"x%d,\",i);");
	    out.println("fprintf(out,\"\\n\");");
	    out.println("for(tb=0;tb<tps;tb++){");
	    out.println("fprintf(out,\"%d,\",tb);");
	    out.println("for(i=1;i<"+varNum+";i++)");
	    out.println("fprintf(out,\"%E,\",result[i][tb]);");
	    out.println("fprintf(out,\"\\n\");");
	    out.println("}");
	    out.println("fclose(out);");

	    out.println("}");

	    //-----------  output marginal probabilities
	    out.println("// output marginal probabilities");
	    out.println("void outputMarginal(){");
	    out.println("FILE *out;");
	    out.println("char buffer[256];");
	    out.println("snprintf(buffer, sizeof(buffer), \"dummy.txt\");");
	    out.println("int idx=0,t,i,b;");
	    out.println("");
	    
	    out.println("snprintf(buffer, sizeof(buffer), \""+modelName+"_M.csv\");");
	    out.println("out=fopen(buffer, \"w\");");
	    out.println("for(t=0;t<tps;t++)");
	    out.println("for(i=1;i<"+varNum+";i++)");
	    out.println("for(b=0;b<xbinNum[i];b++)");
	    out.println("fprintf(out,\"%d,%d,%d,%E\\n\",t,i,b,x[i][t][b]);");
	    out.println("fclose(out);");
	    out.println("}");



	    //-----------  eval the score
	    out.println("double eval(int kvar[]){");
	    out.println("double err=0,value=0;int i,j;");
	    out.println("simulate(kvar);");
	    out.println("for(i=0;i<"+dataTPs.length+";i++){");
	    out.println("for(j=0;j<"+dataVarIDs.length+";j++){");
	    out.println("err=result[X[j]][T[i]]-expdata[i][j];");
	    out.println("value+=weights[j]*err*err;");
	    out.println("}");
	    out.println("}");
	    out.println("return value;");
	    out.println("}");

	    /*
	    // main
	    out.println("int main(int argc, char *argv[]){");
	    out.println("load"+option+"();"); 
	    out.println("if(argc==1){");
	    out.print("int k[]={");
	    for(int i=0;i<freeParIDs.length;i++)
		out.print(pars[freeParIDs[i]].binInit+",");    
	    out.println("};");
	    out.println("printf(\"%E\\n\",eval(k));");
	    out.println("output();");
	    out.println("} else {");
	    out.println("int k["+freeParIDs.length+"];");
	    out.println("for(int i=0;i<"+freeParIDs.length+";i++)");
	    out.println("k[i]=atoi(argv[1][i]);");
	    out.println("");
	    out.println("printf(\"%E\\n\",eval(k));");
	    out.println("output();");
	    out.println("}");
	    out.println("");


	    //out.println("int k[]={3,3,3};"); // carefuel!!
	    //out.println("tff.simulate(k);");
	    out.println("");
	    out.println("}");
	    */

	    out.flush();


	} catch(IOException e){
	    e.printStackTrace();
	}
    }    



    //------------- generating simulation code (Java) -----------

    public void genFFcode(int timeptNum, int simTime, String outputDir, int option, int[]freeParIDs, int[]dataTPs,int[]dataVarIDs, double[]weights, String dataFile, boolean fulldata){
	try{
	    FileOutputStream outfile=new FileOutputStream(modelName+"FF"+option+".java");
	    PrintWriter out=new PrintWriter(outfile);

	    out.println("import java.io.*;");
	    out.println("import java.util.*;");
	    out.println("class "+modelName+"FF"+option+"{");   

	    out.println("int tps="+timeptNum+",tb;");
	    out.println("double dt="+(simTime/(timeptNum*1.0))+";");

	    // ----- gloabal variables

	    // mean of each bin
	    double [][]binM=new double[varNum][];
	    out.println("double [][]M={{-1},");
	    for(int i=1;i<varNum;i++){
		int B=vars[i].binNum;
		binM[i]=new double[B];
		out.print("{");
		for(int j=0;j<B;j++){
		    binM[i][j]=((vars[i].binThresholds)[j+1]+(vars[i].binThresholds)[j])*0.5;
		    out.print(binM[i][j]+",");
		}
		out.println("},");
	    }
	    out.println("};");
	    out.println("");

	    // cpd
	    for(int i=1;i<varNum;i++){
		int n=(vars[i].varIDs).length+(vars[i].parIDs).length+1;
		out.print("double x"+i+"prob");
		for(int j=0;j<n+1;j++) out.print("[]");
		out.print("=new double["+timeptNum+"]");
		//for(int j=0;j<n;j++) out.print("["+vars[i].binNum+"]"); //bug july13
		int tmp_v=(vars[i].varIDs).length;
		int tmp_p=(vars[i].parIDs).length;
		for(int j=0;j<tmp_p;j++)
		    out.print("["+pars[(vars[i].parIDs)[j]].binNum+"]");
		for(int j=0;j<tmp_v;j++)
		    out.print("["+vars[(vars[i].varIDs)[j]].binNum+"]");
		out.print("["+vars[i].binNum+"]");
		out.println(";");
	    }
	    out.println("double result[][];");
	    out.println("double x[][][];");

	    // parameter estimation
	    out.println("");
	    out.print("static int P[]={");
	    for(int i=0;i<freeParIDs.length;i++)
		out.print(freeParIDs[i]+",");
	    out.println("};");

	    out.print("int T[]={");
	    for(int i=0;i<dataTPs.length;i++)
		out.print(dataTPs[i]+",");
	    out.println("};");
	    
	    out.print("int X[]={");
	    for(int i=0;i<dataVarIDs.length;i++)
		out.print(dataVarIDs[i]+",");
	    out.println("};");

	    out.print("double weights[]={");
	    for(int i=0;i<weights.length;i++)
		out.print(weights[i]+",");
	    out.println("};");

	    out.print("int k[]={-1,");
	    for(int i=1;i<parNum;i++)
		out.print(pars[i].binInit+",");	    
	    out.println("};");

	    out.print("int xbinNum[]={-1,");
	    for(int i=1;i<varNum;i++)
		out.print(vars[i].binNum+",");	    
	    out.println("};");

	    out.println("double exp[][];");
	    out.println("");
	    // -------------------


	    // ------- constructor
	    out.println("public "+modelName+"FF"+option+"() throws IOException{");
	    //loadExp	    
	    out.println("");
	    out.println("//loadExp");
	    out.println("FileInputStream fis;");
	    out.println("BufferedReader stdin;");
	    out.println("StringTokenizer st;");
	    out.println("String line;");

	    out.println("double tmp;");

	    out.println("exp=new double["+dataTPs.length+"]["+dataVarIDs.length+"];");
	    out.println("fis = new FileInputStream(\""+dataFile+"\");");
	    out.println("stdin = new BufferedReader(new InputStreamReader(fis));");
	    out.println("line=stdin.readLine(); //header");
	    out.println("line=stdin.readLine(); //init");
	    out.println("for(int j=0;j<"+dataTPs.length+";j++){");
	    out.println("int ctr=0;");
	    out.println("line=stdin.readLine();");
	    out.println("st = new StringTokenizer(line, \"\\t, \");");
	    out.println("st.nextToken(); //time");
	    if(fulldata){
		out.println("for(int k=1;k<"+varNum+";k++){");
		out.println("if(ctr=="+dataVarIDs.length+") break;");
		out.println("tmp=Double.parseDouble(st.nextToken());");
		out.println("if(k==X[ctr]){");
		out.println("exp[j][ctr]=tmp;");
		out.println("ctr++;");
		out.println("}");
		out.println("}");
	    } else {
		out.println("for(int k=0;k<"+dataVarIDs.length+";k++){");
		out.println("tmp=Double.parseDouble(st.nextToken());");
		out.println("exp[j][ctr]=tmp;");
		out.println("}");
	    }
	    out.println("}");
	    out.println("");


	    out.println("");
	    // load cpd
	    out.println("// load CPD");
	    out.println("int idx;double p;");
	    out.println("for(tb=0;tb<"+timeptNum+";tb++){");
	    String tag="P";
	    if(modelName.equals("brown")) tag="XBNT"; // [patch] brown's model use XBNT lable
	    for(int i=1;i<varNum;i++){
		out.println("fis=new FileInputStream(\""+outputDir+"/tables/"+modelName+""+tag+"x"+i+"T\"+tb+\".txt\");");
		out.println("stdin = new BufferedReader(new InputStreamReader(fis));");
		out.println("while((line=stdin.readLine())!=null){");
		out.println("st = new StringTokenizer(line, \" \");");
		out.println("idx=Integer.parseInt(st.nextToken());");
		out.println("p=Double.parseDouble(st.nextToken());");
		int tmp_v=(vars[i].varIDs).length;
		int tmp_p=(vars[i].parIDs).length;
		out.println("int vb=idx%"+vars[i].binNum+";");
		out.println("idx=idx/"+vars[i].binNum+";");	
		for(int j=tmp_v-1;j>=0;j--){
		    int tmp_b=vars[(vars[i].varIDs)[j]].binNum;
		    out.println("int vi"+j+"=idx%"+tmp_b+";");
		    out.println("idx=idx/"+tmp_b+";");	    
		}
		for(int j=tmp_p-1;j>=0;j--){
		    int tmp_b=pars[(vars[i].parIDs)[j]].binNum;
		    out.println("int pi"+j+"=idx%"+tmp_b+";");
		    out.println("idx=idx/"+tmp_b+";");	    
		}
		out.print("x"+i+"prob[tb]");
		for(int j=0;j<tmp_p;j++) out.print("[pi"+j+"]");
		for(int j=0;j<tmp_v;j++) out.print("[vi"+j+"]");
		out.println("[vb]=p;");
		out.println("}");
		out.println("");
		out.println("stdin.close();fis.close();");
	    }
	    out.println("result=new double["+varNum+"][tps+1];");
	    
	    out.println("}");
	    out.println("}");
	    

	    // -------- simulation
	    out.println("public void simulate(int kvar[]){");
	    out.println("for(int i=0;i<"+freeParIDs.length+";i++) k[P[i]]=kvar[i];");

	    // init
	    int binMax=0;
	    for(int i=1;i<varNum;i++)		
		if(vars[i].binNum>binMax)
		    binMax=vars[i].binNum;
	    
	    out.println("");
	    out.println("//init");
	    out.println("x=new double["+varNum+"][tps+1]["+binMax+"];");	    
	    out.println("result=new double["+varNum+"][tps+1];");
	    out.println("");
	    	    

	    int N=1;
	    for(int i=1;i<varNum;i++){
		if(vars[i].init>=0){
		    out.println("x["+i+"][0]["+vars[i].binInit+"]=1;");
		    N++;
		}
	    }
	    // expect
	    out.println("");
	    out.println("int b,i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15;");
	    out.println("double sum=0,margin=1,norm=0;");
	    out.println("");
	    out.println("for(int t=0;t<"+timeptNum+";t++){"); //<--
	    //out.println("tb=(int)(t/block);");
	    out.println("tb=t;");
	    out.println("// normalization");
	    out.println("for(int i=1;i<"+N+";i++){");
	    out.println("norm=0; ");
	    out.println("for(b=0;b<xbinNum[i];b++) norm+=x[i][t][b]; if(norm==0) norm=1;");
	    out.println("for(b=0;b<xbinNum[i];b++) x[i][t][b]/=norm;");
	    out.println("}");
	    out.println("");
	    
	    // /* for reduced model !!!

	    for(int i=N;i<varNum;i++){
		out.println("// x"+i);
		out.println("for(b=0;b<"+vars[i].binNum+";b++){");
		//out.println("result["+i+"][t]+=x["+i+"][t][b]*b;");
		int tmp_p=(vars[i].parIDs).length;
		int tmp_v=(vars[i].varIDs).length;
		for(int j=0;j<tmp_v;j++) 
		    out.println("for(i"+j+"=0;i"+j+"<"+vars[(vars[i].varIDs)[j]].binNum+";i"+j+"++)");
      		out.print("x["+i+"][t][b]+=x"+i+"prob[tb]");
		for(int j=0;j<tmp_p;j++) out.print("[k["+(vars[i].parIDs)[j]+"]]");
		for(int j=0;j<tmp_v;j++) out.print("[i"+j+"]");
		out.print("[b]");
		for(int j=0;j<tmp_v;j++) 
		    out.print("*x["+(vars[i].varIDs)[j]+"][t][i"+j+"]");
		out.println(";");
		out.println("result["+i+"][t]+=x["+i+"][t][b]*M["+i+"][b];");
		out.println("}");
		
	    }

	    out.println("");	
	    out.println("// normalization2");
	    out.println("for(int i="+N+";i<"+varNum+";i++){");
	    out.println("norm=0; ");
	    out.println("for(b=0;b<xbinNum[i];b++) norm+=x[i][t][b];if(norm==0) norm=1;");
	    out.println("for(b=0;b<xbinNum[i];b++) x[i][t][b]/=norm;");
	    out.println("}");
	    out.println("");	
	    // ----------- */

	    
	    for(int i=1;i<N;i++){
		out.println("// x"+i);
		out.println("for(b=0;b<"+vars[i].binNum+";b++){");
		out.println("result["+i+"][t]+=x["+i+"][t][b]*M["+i+"][b];");
		//out.println("result["+i+"][t]+=x["+i+"][t][b]*b;");
		int tmp_p=(vars[i].parIDs).length;
		int tmp_v=(vars[i].varIDs).length;
		for(int j=0;j<tmp_v;j++) 
		    out.println("for(i"+j+"=0;i"+j+"<"+vars[(vars[i].varIDs)[j]].binNum+";i"+j+"++)");
      		out.print("x["+i+"][t+1][b]+=x"+i+"prob[tb]");
		for(int j=0;j<tmp_p;j++) out.print("[k["+(vars[i].parIDs)[j]+"]]");
		for(int j=0;j<tmp_v;j++) out.print("[i"+j+"]");
		out.print("[b]");
		for(int j=0;j<tmp_v;j++) 
		    out.print("*x["+(vars[i].varIDs)[j]+"][t][i"+j+"]");
		out.println(";");
		out.println("}");
		
	    }
	    out.println("}");
	    out.println("}");
	    out.println("");
	    
	    //----------  output
	    out.println("// output");
	    out.println("public void output(String postfix) throws IOException{");
	    out.println("FileOutputStream outfile = new FileOutputStream(\""+modelName+"_Te\"+postfix+\".csv\");");
	    out.println("PrintWriter out = new PrintWriter(outfile);");
	    out.println("FileOutputStream outfile2 = new FileOutputStream(\""+modelName+"_De\"+postfix+\".txt\");");
	    out.println("PrintWriter out2 = new PrintWriter(outfile2);");

	    out.println("out.print(\"Time,\");");
	    out.println("for(int i=1;i<"+varNum+";i++)");
	    out.println("out.print(\"x\"+i+\",\");");
	    out.println("out.println();");
	    out.println("for(int t=0;t<tps;t++){");
	    out.println("out.print(t+\",\");");
	    out.println("for(int i=1;i<"+varNum+";i++){");
	    out.println("out.print(result[i][t]+\",\");");
	    out.println("out2.print((int)result[i][t]+\" \");");
	    out.println("}");
	    out.println("out.println();");
	    out.println("out2.println();");
	    out.println("}");
	    out.println("out.flush();");
	    out.println("out2.flush();");
	    
	    out.println("}");

	    //-----------  output marginal probabilities
	    out.println("// output marginal probabilities");
out.println("public void outputMarginal(String postfix) throws IOException{");
	    out.println("FileOutputStream outfile = new FileOutputStream(\""+modelName+"_Me\"+postfix+\".csv\");");
	    out.println("PrintWriter out = new PrintWriter(outfile);");
	    out.println("for(int t=0;t<tps;t++)");
	    out.println("for(int i=1;i<"+varNum+";i++)");
	    out.println("for(int b=0;b<xbinNum[i];b++)");
	    out.println("out.println(t+\",\"+i+\",\"+b+\",\"+x[i][t][b]);");
	    out.println("out.flush();");
	    out.println("}");



	    //-----------  eval the score
	    out.println("public double eval(int kvar[]){");
	    out.println("simulate(kvar);");
	    out.println("double err=0,value=0;");
	    out.println("for(int i=0;i<"+dataTPs.length+";i++){");
	    out.println("for(int j=0;j<"+dataVarIDs.length+";j++){");
	    out.println("err=result[X[j]][T[i]]-exp[i][j];");
	    out.println("value+=weights[j]*err*err;");
	    out.println("}");
	    out.println("}");
	    out.println("return value;");
	    out.println("}");

	    // main
	    out.println("public static void main(String [] args) throws Exception{");
	    out.println(modelName+"FF"+option+" tff=new "+modelName+"FF"+option+"();");	    
	    out.println("if(args.length==0){");
	    out.print("int k[]={");
	    for(int i=0;i<freeParIDs.length;i++)
		out.print(pars[freeParIDs[i]].binInit+",");    
	    out.println("};");
	    out.println("System.out.println(tff.eval(k));");
	    out.println("tff.output(\"\");");
	    out.println("} else {");
	    out.println("int k[]=new int["+freeParIDs.length+"];");
	    out.println("for(int i=0;i<"+freeParIDs.length+";i++)");
	    out.println("k[i]=Integer.parseInt(args[0].substring(i,i+1));");
	    out.println("");
	    out.println("System.out.println(tff.eval(k));");
	    out.println("tff.output(args[0]);");
	    out.println("}");
	    out.println("");


	    //out.println("int k[]={3,3,3};"); // carefuel!!
	    //out.println("tff.simulate(k);");
	    out.println("");
	    out.println("}");
	    out.println("}");

	    out.flush();


	} catch(IOException e){
	    e.printStackTrace();
	}
    }


    // -------- utils ------------
    //private


    // -------- testing ----------
    public static void main(String[] args) throws IOException{
	SPAmodel m=new SPAmodel("complement");


	m.loadVar("./complement/var-slim-e-pc.csv");
 	m.loadPar("./complement/par-slim-e-pc.csv");
 	m.loadEqn("./complement/ode-slim-e-pc.txt");
	m.updateEqn("./complement/ode-pc.txt");
	System.out.println("varNum:"+m.varNum);
	System.out.println("parNum:"+m.parNum);
	for(int i=1;i<m.varNum;i++)
	    System.out.println(m.vars[i]);
	for(int i=1;i<m.parNum;i++)
	    System.out.println(m.pars[i]);

	//m.genTESTcode(100,0.01,12600,"."); 	// determin the range of e-nodes
	//m.genMPIcode(10,1000,100,0.01,12600,"/home/course/cs6280/dSFMT-src-2.0/dSFMT.c","/home/course/cs6280/complement/mpict");
	//m.genBATcode(1000000,100,0.01,200,"/home/course/cs6280/dSFMT-src-2.0/dSFMT.c","/home/course/cs6280/complement/batct");
	//m.genBATcode2(1000,100,0.01,12600,"/home/course/cs6280/dSFMT-src-2.0/dSFMT.c","/home/course/cs6280/complement/batct");
	//m.processN(1,50,100,"/home/course/cs6280/complement/","/home/course/cs6280/complement/ctn",0);
	m.genNORMcode(1000,100,0.01,12600,"/home/course/cs6280/dSFMT-src-2.0/dSFMT.c",".");
	//m.genNORMcode2(1000,100,0.01,500,"/home/course/cs6280/dSFMT-src-2.0/dSFMT.c",".");
	int P[]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
		 21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,
		 41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,
		 61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,
		 81,82,83,84,85,86,87,88,89,90,91,92};
	int T[]={1800/126,3600/126,5400/126,7200/126,9000/126,10800/126,12600/126};
	int X[]={1,2,3,4};
	double weights[]={1,1,1,1};
	
	//m.genMPSAcode(10,100,0.01,500,"/home/course/cs6280/dSFMT-src-2.0/dSFMT.c",".",P,T,X,weights,"/home/course/cs1101s/spa/m201Noisy-tight.txt",34,1.0);
	m.genFFcode(100,12600,"/home/svu/g0601796/spa/complement/ctn",0,P,T,X,weights,"/home/svu/g0601796/spa/complement/complement-data.txt",false);	
	//m.genFFcode(100,10,"/home/course/cs6280/m201/ctn",0,P,T,X,weights);


    }

}

class Variable{
    int ID;
    String name;
    double init; // init < 0 means intermediat variables (assignment type)
    double lowerBound;
    double upperBound;
    int binNum;
    int binInit;
    double []binSizes;
    double []binThresholds;
    String eqn; // ODE expresion
    int varIDs[]; // variable parents
    int parIDs[]; // parameter parants

    public Variable(String line, int ID){
	this.ID=ID;
	StringTokenizer st=new StringTokenizer(line, ",");
	name=st.nextToken();
	init=Double.parseDouble(st.nextToken());
	lowerBound=Double.parseDouble(st.nextToken());
	upperBound=Double.parseDouble(st.nextToken());
	binNum=Integer.parseInt(st.nextToken());
	binSizes=new double[binNum];
	if(st.hasMoreTokens()){
	    // add normalization code 071909
	    double sum=0.0;
	    for(int i=0;i<binNum;i++){
		binSizes[i]=Double.parseDouble(st.nextToken());
		sum+=binSizes[i];
	    }
	    if(sum!=0)
		for(int i=0;i<binNum;i++)
		    binSizes[i]=binSizes[i]/sum;
	    else{
		System.out.println("wrong bin size");
	    }
	} else {
	    for(int i=0;i<binNum;i++)
		binSizes[i]=1.0/binNum; // corrected @ Oct 06, 09
	}

	double tmp_init=(init-lowerBound)/(upperBound-lowerBound);
	double tmp_b=0;
	binInit=binNum-1;
	for(int j=0;j<binNum-1;j++){
	    tmp_b+=binSizes[j];
	    if(tmp_init<tmp_b){
		binInit=j;
		break;
	    }
	}

	binThresholds=new double[binNum+1];
	double tmp_r=upperBound-lowerBound;
	tmp_b=lowerBound;
	for(int j=0;j<binNum;j++){
	    binThresholds[j]=tmp_b;
	    tmp_b+=binSizes[j]*tmp_r;
	}
	binThresholds[binNum]=upperBound;
    }
    
    public void setEquation(String exp, int[]v, int[]p){
	eqn=exp;
	varIDs=v;
	parIDs=p;
    }

    public void updateEquation(String exp){
	eqn=exp;
    }

    public String toString(){
	String s=ID+" "+name+" "+init+" "+lowerBound+" "+upperBound+" "+binNum;
	for(int i=0;i<binNum;i++) s+=" "+binSizes[i];
	//s+=" "+eqn;
	for(int i=0;i<varIDs.length;i++) s+=" x"+varIDs[i];
	for(int i=0;i<parIDs.length;i++) s+=" k"+parIDs[i];
	int tmp=varIDs.length+parIDs.length;
	s+=" "+tmp;
	return s;
    }

}

class Parameter{
    int ID;
    String name;
    double init;
    double lowerBound;
    double upperBound;
    int binNum;
    int binInit;
    double []binSizes;
    public Parameter(String line, int ID){
	this.ID=ID;
	StringTokenizer st=new StringTokenizer(line, ",");
	name=st.nextToken();
	init=Double.parseDouble(st.nextToken());
	lowerBound=Double.parseDouble(st.nextToken());
	upperBound=Double.parseDouble(st.nextToken());
	binNum=Integer.parseInt(st.nextToken());
	binSizes=new double[binNum];
	if(st.hasMoreTokens()){
	    for(int i=0;i<binNum;i++)
		binSizes[i]=Double.parseDouble(st.nextToken());
	} else {
	    for(int i=0;i<binNum;i++)
		binSizes[i]=1.0/binNum; // corrected @ Oct 06, 09
	}

	double tmp_init=(init-lowerBound)/(upperBound-lowerBound);
	double tmp_b=0;
	binInit=binNum-1;
	for(int j=0;j<binNum-1;j++){
	    tmp_b+=binSizes[j];
	    if(tmp_init<tmp_b){
		binInit=j;
		break;
	    }
	}


    }
    public String toString(){
	String s=name;
	s+=" "+binInit;
	return s;
    }
}
