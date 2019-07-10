package strbal;

import java.util.ArrayList;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;

import strbal.Var.Type;

public class Algorithm {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		FileInputStream fstream = new FileInputStream(args[0]);
	    // Get the object of DataInputStream
	    DataInputStream in = new DataInputStream(fstream);
	        BufferedReader br = new BufferedReader(new InputStreamReader(in));
	    String strLine;
		ArrayList <String> fomula = new ArrayList <String>();
	    //Read File Line By Line
	    while ((strLine = br.readLine()) != null)   {
	      // Print the content on the console
	    	if (strLine.startsWith("#")) continue;
	    	if (strLine.startsWith("EQU") || strLine.startsWith("AEQU")) {
	    		if (!Equation.define(strLine,fomula)) {
	    			System.out.println("Parsing error: " + strLine);
	    			return;
	    		}
	    	}
	    	else if (!Var.define(strLine)) {
	    		System.out.println("Parsing error: " + strLine);
	    		return;
	    	}
	    }

	    //Close the input stream
	    in.close();
		int warps = Integer.parseInt(args[1]);
		int nrAccess = Integer.parseInt(args[4]);
		
		int thread_step = (16384 * 3 - 2048) / (args[3].equals("double") ? 8 : 4) / 32 - 1;
		int fullThreshold = Integer.parseInt(args[5]);
		Schedule sched = Schedule.getSchedule(thread_step, fullThreshold, warps, nrAccess);
		if (sched == null) {
			System.out.println("Failed to generate a schedule");
			return;
		}
		dump(sched, args[2], thread_step);
	}
	
	private static void dump(Schedule sched, String fileName, int thread_step) throws IOException { 		
		int i;
		
		PrintWriter out = new PrintWriter(new FileOutputStream(fileName));

		out.println("#define THREAD_STEP " + thread_step);
	
		out.print("__device__ real const_min[" + Var.count(Type.CONST) + "] = {");
		
		for (i = 0; i < Var.count(Type.CONST); i++) {
			out.print(Var.get(Type.CONST, i).getMin());
			if (i < Var.count(Type.CONST) - 1) out.print(", ");
		}
		out.println("};");
		
		out.print("__device__ real const_max[" + Var.count(Type.CONST) + "] = {");
		
		for (i = 0; i < Var.count(Type.CONST); i++) {
			out.print(Var.get(Type.CONST, i).getMax());
			if (i < Var.count(Type.CONST) - 1) out.print(", ");
		}
		out.println("};");

		out.print("__device__ int const_bins[" + Var.count(Type.CONST) + "] = {");
		
		for (i = 0; i < Var.count(Type.CONST); i++) {
			out.print(Var.get(Type.CONST, i).getBins());
			if (i < Var.count(Type.CONST) - 1) out.print(", ");
		}
		out.println("};");

		out.print("__device__ real var_min[" + (Var.count(Type.VARPRE) + Var.count(Type.FAT)) + "] = {");
		
		for (i = 0; i < Var.count(Type.VARPRE) + Var.count(Type.FAT); i++) {
			if (i < Var.count(Type.VARPRE)) out.print(Var.get(Type.VARPRE, i).getMin());
			else out.print(Var.get(Type.FAT, i).getMin());
			if (i < Var.count(Type.VARPRE) + Var.count(Type.FAT) - 1) out.print(", ");
		}
		out.println("};");
		
		out.print("__device__ real var_max[" + (Var.count(Type.VARPRE) + Var.count(Type.FAT)) + "] = {");
		
		for (i = 0; i < Var.count(Type.VARPRE) + Var.count(Type.FAT); i++) {
			if (i < Var.count(Type.VARPRE)) out.print(Var.get(Type.VARPRE, i).getMax());
			else out.print(Var.get(Type.FAT, i).getMax());
			if (i < Var.count(Type.VARPRE) + Var.count(Type.FAT) - 1) out.print(", ");
		}
		out.println("};");

		out.print("__device__ int var_bins[" + (Var.count(Type.VARPRE) + Var.count(Type.FAT)) + "] = {");
		
		for (i = 0; i < Var.count(Type.VARPRE) + Var.count(Type.FAT); i++) {
			if (i < Var.count(Type.VARPRE)) out.print(Var.get(Type.VARPRE, i).getBins());
			else out.print(Var.get(Type.FAT, i).getBins());
			if (i < Var.count(Type.VARPRE) + Var.count(Type.FAT) - 1) out.print(", ");
		}
		out.println("};");

		out.print("__device__ int var_initbin[" + Var.count(Type.VARPRE) + "] = {");
		
		for (i = 0; i < Var.count(Type.VARPRE); i++) {
			out.print(Var.get(Type.VARPRE, i).getInitBin());
			if (i < Var.count(Type.VARPRE) - 1) out.print(", ");
		}
		out.println("};");
		
		out.println("\n__device__ void compute(int t, volatile struct struct_trace_data *this_trace, volatile real *local_data) {");
		out.println("\treal xx_delta, halfF1, halfF2, F3, F4;");
		for (int warp = 0; warp < sched.computeWarps; warp++) {
			out.println("\tif (threadIdx.y == " + warp + ") {");
			for (Step step : sched.computeWarpSchedules.get(warp)) {
				out.println(step);
			}
			out.println("\t}");
		}
		for (int warp = 0; warp < sched.accessWarps; warp++) {
			out.println("\tif (threadIdx.y == " + (warp + sched.computeWarps) + ") {");
			for (Step step : sched.accessWarpSchedules.get(warp)) {
				out.println(step);
			}
			out.println("\t}");
		}
		
		out.println("}");

		out.close();
	}
}
