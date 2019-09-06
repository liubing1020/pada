package strbal;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.StringTokenizer;

public class Var {

	public enum Type {CONST, FAT, VARPRE, VARPOST};
	
	protected static HashMap<Integer, Var> 	consts = new HashMap<Integer, Var>(), 
											fats = new HashMap<Integer, Var>(),
											preVars = new HashMap<Integer, Var>(),
											postVars = new HashMap<Integer, Var>();
	
	public final Type type;
	public final int nr;
	private double min, max;
	private int bins, initBin;
	
	public String toString() {
		switch (type) {
			case CONST: return "K(" + nr + ")";
			case FAT: return "A(" + nr + ")";
			case VARPRE: return "X(" + nr + ")";
			case VARPOST: return "Y(" + nr + ")";
			default: return null;
		}
	}
	
	private Var(Type t, int n, double mi, double ma, int b, int ib) {
		type = t;
		nr = n;
		min = mi;
		max = ma;
		bins = b;
		initBin = ib;		
	}

	public static int count(Type type) {
		switch (type) {
			case CONST: return consts.size();
			case FAT: return fats.size();
			case VARPRE: return preVars.size();
			case VARPOST: return postVars.size();
		}
		return 0;
	}
	
	public Var getInitial() {
		if (type != Type.VARPOST) return this;
		return preVars.get(nr);
	}
	
	public Var getFinal() {
		if (type != Type.VARPRE) return this;
		return postVars.get(nr);
	}
	
	public static Var get(Type type, int nr) {
		switch (type) {
			case CONST: return consts.get(nr);
			case FAT: return fats.get(nr);
			case VARPRE: return preVars.get(nr);
			case VARPOST: return postVars.get(nr);
		}
		return null;
	}
	
	public static boolean define(String strLine) {
		StringTokenizer st = new StringTokenizer(strLine, " (,)");
		if (!st.hasMoreTokens()) return true;
		String type = st.nextToken();
		if (!st.hasMoreTokens()) return false;
		int nr = Integer.parseInt(st.nextToken());
		if (!st.hasMoreTokens()) return false;
		double min = Double.parseDouble(st.nextToken());
		if (!st.hasMoreTokens()) return false;
		double max = Double.parseDouble(st.nextToken());
		if (!st.hasMoreTokens()) return false;
		int bins = Integer.parseInt(st.nextToken());
		if (type.startsWith("CONST")) {
			if (consts.containsKey(nr)) return false;
			Var cst = new Var(Type.CONST, nr, min, max, bins, -1);
			consts.put(nr, cst);
			return true;
		} else if (type.startsWith("AVAR")) {
			if (fats.containsKey(nr)) return false;
			Var var = new Var(Type.FAT, nr, min, max, bins, -1);
			fats.put(nr, var);
			return true;
		} else if (type.startsWith("VAR")) {
			if (!st.hasMoreTokens()) return false;
			int initBin = Integer.parseInt(st.nextToken());
			if (preVars.containsKey(nr) || postVars.containsKey(nr)) return false;
			Var var = new Var(Type.VARPRE, nr, min, max, bins, initBin);
			preVars.put(nr, var);
			var = new Var(Type.VARPOST, nr, min, max, bins, initBin);
			postVars.put(nr, var);
			return true;
		}
		return false;
	}

	public static int problemSize() {
		return consts.size() + postVars.size();
	}

	
/*	public static void placeAll(old_Memory memory) {
		ArrayList<Var> list = new ArrayList<Var>(consts.values());
		list.addAll(preVars.values());
		int warp = 0, reg = 0, total = 0;
		for (Var var : list) {
			memory.putReg(warp, reg, var);
			reg++;
			total++;
			if (total >= (warp + 1) * list.size() / memory.nrWarps) {
				reg = 0;
				warp++;
			}
		}
	}*/

	public int hashCode() {
		return nr;
	}

	public double getMin() {
		return min;
	}
	
	public double getMax() {
		return max;
	}

	public int getBins() {
		// TODO Auto-generated method stub
		return bins;
	}

	public int getInitBin() {
		return initBin;
	}

	public String globalLoc() {
		switch (type) {
			case VARPRE: return "this_trace -> traceVars[" + nr + "]"; 
			case VARPOST: return "this_trace -> tracePost[" + nr + "]"; 
			case CONST: return "this_trace -> traceConsts[" + nr + "]";
			case FAT: return "this_trace -> traceFats[" + nr + "]";
			default: return "?????";
		}
	}
}
