package strbal;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.StringTokenizer;
import java.util.Vector;

import strbal.Var.Type;

public class Equation {
	private String formula;
	private Var returnVar;
	private ArrayList<Var> dependencies;
	private ArrayList<Var> binObjs;

	public static ArrayList<Equation> equs = new ArrayList<Equation>(), lateEqus = new ArrayList<Equation>();

	public String toString() {
		return returnVar + " = " + formula;
	}
	
	public ArrayList<Var> getUsedDependencies() {
		return new ArrayList<Var>(dependencies);
	}

	public ArrayList<Var> getAllAccessed() {
		ArrayList<Var> allAccessed= new ArrayList<Var>(dependencies);
		if (!allAccessed.contains(returnVar)) allAccessed.add(returnVar);
		return allAccessed;
	}

	public Var getResultVar() {
		return returnVar;
	}
	
	private Equation(String form, Var ret, ArrayList<Var> dep, ArrayList<Var> bins) {
		formula = form;
		returnVar = ret;
		dependencies = dep;
		binObjs = bins;
	}
	
	public static boolean define(String strLine,ArrayList <String> form) {
		StringTokenizer st = new StringTokenizer(strLine, ",");
		String firstToken = st.nextToken();
		int equNr = -1;
		Var ret = null;
		if (firstToken.startsWith("AEQUATION")) {
			equNr = Integer.parseInt(firstToken.substring(10));
			ret = Var.get(Type.FAT, equNr);
		}
		if (firstToken.startsWith("EQUATION")) {
			equNr = Integer.parseInt(firstToken.substring(9));
			ret = Var.get(Type.VARPOST, equNr);
		}
		if (equNr < 0) return false; 
		if (!st.hasMoreTokens()) return false;
		String formula = st.nextToken().replaceAll("XX", "xx_delta");
		String formula_t = formula.replaceAll("C\\(", "(real)(");
		
		if (firstToken.startsWith("EQUATION")){
				formula = "(xx_delta = X(" + equNr + "), halfF1 = ((real)HLFDT) * (" + formula + " ), " +
				  "xx_delta = X(" + equNr + ") + halfF1, halfF2 = ((real) HLFDT) * (" + formula + " ), " + 
				  "xx_delta = X(" + equNr + ") + halfF2, F3 = ((real) DT) * (" + formula + " ), " + 
				  "xx_delta = X(" + equNr + ") + F3, F4 = ((real)DT) * (" + formula + " ), " + 
				  "X(" + equNr + ") + (((real)2) * (halfF1 + F3) + (real)4 * halfF2 + F4) * ((real)1.0/6.0))"; 
		}	
		formula = formula.replaceAll("PLUS", "+");
		formula = formula.replaceAll("DIV", "/");
		formula = formula.replaceAll("MUL", "*");
		formula = formula.replaceAll("MINUS", "-");
		formula = formula.replaceAll("ZERO", "0");
		formula = formula.replaceAll("OPEN", "(");
		formula = formula.replaceAll("CLOSE", ")");
		formula = formula.replaceAll("AX", "A");
		
		String binObjs = st.nextToken();
		binObjs = binObjs.substring(0, binObjs.length() - 1);
		binObjs = binObjs.replaceAll("AX", "A");
		
		ArrayList<Var> dep = new ArrayList<Var>();

		for (String token : formula.split("[^\\d()AXK]")) {
			if (token.startsWith("A")) {
				int nr = Integer.parseInt(token.substring(2, token.length() - 1));
				dep.add(Var.get(Type.FAT, nr));
			} else if (token.startsWith("K")) {
				int nr = Integer.parseInt(token.substring(2, token.length() - 1));
				dep.add(Var.get(Type.CONST, nr));
			} else if (token.startsWith("X")) {
				int nr = Integer.parseInt(token.substring(2, token.length() - 1));
				dep.add(Var.get(Type.VARPRE, nr));
			}			
		
		}

		ArrayList<Var> bins = new ArrayList<Var>();
		
		for (String token : binObjs.split("[^\\d()AXK]")) {
			if (token.startsWith("A")) {
				int nr = Integer.parseInt(token.substring(2, token.length() - 1));
				bins.add(Var.get(Type.FAT, nr));
			} else if (token.startsWith("K")) {
				int nr = Integer.parseInt(token.substring(2, token.length() - 1));
				bins.add(Var.get(Type.CONST, nr));
			} else if (token.startsWith("X")) {
				int nr = Integer.parseInt(token.substring(2, token.length() - 1));
				bins.add(Var.get(Type.VARPRE, nr));
			}			
		}
		if (firstToken.startsWith("AEQUATION")) {
			bins.add(Var.get(Type.FAT, equNr));
		} else {
			bins.add(Var.get(Type.VARPOST, equNr));
		}
		HashSet<Var> h = new HashSet<Var>(dep);
		dep.clear();
		dep.addAll(h);
		
		equs.add(new Equation(formula, ret, dep, bins));
		ArrayList<Var> arRet = new ArrayList<Var>();
		arRet.add(ret);
		if (ret.type == Type.FAT) {
			lateEqus.add(new Equation("if (t == 0) bin_pack(this_trace,  " + ret.nr + ", " + ret + ")", null, arRet, null));
		} else {
			lateEqus.add(new Equation("if (t == STEPS - 1) bin_pack(this_trace,  " + ret.nr + ", " + ret + ")", null, arRet, null));
		}
/*		if (firstToken.startsWith("EQUATION")) {
			ArrayList<Var> arCopy = new ArrayList<Var>();
			arCopy.add(Var.get(Type.VARPOST, equNr));
			arCopy.add(Var.get(Type.VARPRE, equNr));
			lateEqus.add(new Equation("X(" + equNr + ") = Y(" + equNr + ");", null, arCopy, null));
		}*/
		return true;
	}

	public String bindForm(Vector<Var> varBinding) {
		String modifiedFormula = ((returnVar != null) ? returnVar + " = "  : "") + formula;
		int i;
		for (i = 0; i < varBinding.size(); i++) {
			Var var = varBinding.get(i);
			if (var == null) continue;
			modifiedFormula = modifiedFormula.replace(var.toString(), "local_data[" + i + "]"); //!!
		}
		modifiedFormula = modifiedFormula.replace("C", "(real)");
		return ((returnVar != null && returnVar.type == Type.FAT) ? "if (t == 0) " : "") + modifiedFormula; //!!
	}

	public static ArrayList<Equation> getEqus() {
		return new ArrayList<Equation>(equs);
	}

	public static ArrayList<Equation> getLateEqus() {
		return new ArrayList<Equation>(lateEqus);
	}
	
}
