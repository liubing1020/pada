package strbal;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Vector;

import strbal.Var.Type;

public class Schedule {
	private HashMap<Var, Integer> varLiveness;
	public Vector<ArrayList<Step>> computeWarpSchedules, accessWarpSchedules;
	private Vector<Var> memory;
	boolean[] used, loadable;
	int computeWarps, accessWarps;
	
	private Schedule(int computeWarps, int accessWarps, int mem_size, HashMap<Var, Integer> live) {
		varLiveness = live;
		computeWarpSchedules = new Vector<ArrayList<Step>>(computeWarps);
		computeWarpSchedules.setSize(computeWarps);
		accessWarpSchedules = new Vector<ArrayList<Step>>(accessWarps);
		accessWarpSchedules.setSize(accessWarps);
		for (int warp = 0; warp < computeWarps; warp++) computeWarpSchedules.set(warp, new ArrayList<Step>());
		for (int warp = 0; warp < accessWarps; warp++) accessWarpSchedules.set(warp, new ArrayList<Step>());
		memory = new Vector<Var>();
		memory.setSize(mem_size);
		used = new boolean[mem_size];
		loadable = new boolean[mem_size];
		for (int i = 0; i < memory.size(); i++) {
			memory.set(i, null);
			used[i] = false;
			loadable[i] = true;
		}
		this.computeWarps = computeWarps;
		this.accessWarps = accessWarps;
	}
	
	/* counts var utilizations, only for fat variables */
	
	private static HashMap<Var, Integer> computeLiveness(ArrayList<Equation> equs, ArrayList<Equation> lateEqus) {
		HashMap<Var, Integer> limitedLife = new HashMap<Var, Integer>();
		for (Equation equ : equs) {
			// count var utilizations once per equation
			for (Var var : equ.getUsedDependencies()) {
				if (limitedLife.containsKey(var)) {
					limitedLife.put(var, limitedLife.get(var) + 1);
				} else {
					if (var.type == Type.FAT) limitedLife.put(var, 1);
				}
			}
		}
		for (Equation equ : lateEqus) {
			// count var utilizations once per equation
			for (Var var : equ.getUsedDependencies()) {
				if (limitedLife.containsKey(var)) {
					limitedLife.put(var, limitedLife.get(var) + 1);
				} else {
					if (var.type == Type.FAT) limitedLife.put(var, 1);
				}
			}
		}
		return limitedLife;
	}

	private int computeMoves(Equation equ) {
		int moves = 0;
		for (Var v : equ.getUsedDependencies()) {
			if (!memory.contains(v)) moves++;
		}
		return moves;
	}
	
	private void addSyncPoint() {
		for (ArrayList<Step> warpSchedule : computeWarpSchedules) {
			warpSchedule.add(Step.newSyncPoint());
		}
		for (ArrayList<Step> warpSchedule : accessWarpSchedules) {
			warpSchedule.add(Step.newSyncPoint());
		}
		for (int i = 0; i < memory.size(); i++) {
			if (!used[i]) loadable[i] = true;
			used[i] = false;
		}
	}
	
	private int occupancy() {
		int ocup = 0;
		for (int i = 0; i < memory.size(); i++) {
			if (used[i]) ocup++;
		}
		return ocup;
	}
	
	public static Schedule getSchedule(int mem_size, int threshold, int warps, int mem_warps) {
		System.out.println("Enter scheduling step");
		ArrayList<Equation> equs = Equation.getEqus();

		Schedule sched = new Schedule(warps - mem_warps, mem_warps, mem_size, computeLiveness(equs, Equation.getLateEqus()));
		
		for (ArrayList<Step> warpSchedule : sched.computeWarpSchedules) {
			warpSchedule.add(Step.newSyncPoint());
		}
		
		while (equs.size() > 0) {
			Equation best = null;
			int moves = Integer.MAX_VALUE;
			for (Equation equ : equs) {
				int this_moves = sched.computeMoves(equ);
				if (this_moves < moves) best = equ;
				if (this_moves == 0) break;
			}
			sched.scheduleEqu(best);
			equs.remove(best);
			if (sched.occupancy() > threshold * mem_size / 100) sched.addSyncPoint();
		}

		sched.addSyncPoint();

		equs = Equation.getLateEqus();
		
		while (equs.size() > 0) {
			Equation best = null;
			int moves = Integer.MAX_VALUE;
			for (Equation equ : equs) {
				int this_moves = sched.computeMoves(equ);
				if (this_moves < moves) best = equ;
				if (this_moves == 0) break;
			}
			sched.scheduleEqu(best);
			equs.remove(best);
			if (sched.occupancy() > threshold * mem_size / 100) sched.addSyncPoint();
		}
		
		sched.addSyncPoint();

		for (ArrayList<Step> warpSchedule : sched.accessWarpSchedules) {
			warpSchedule.add(Step.newSyncPoint());
		}

		sched.commitAll();
		
		sched.addSyncPoint();
		
		return sched;
	}
	
	// ****************************************************************************
	// Actual scheduling
	// ****************************************************************************
	
	int whichCompute = 0, whichAccess = 0;
	
	int scanPointer = 0;
	
	private void commitAll() {
		for (int i = 0; i < memory.size(); i++) {
			if (memory.get(i) == null) continue;
			if (memory.get(i).type == Type.VARPRE || memory.get(i).type == Type.CONST) continue;
			accessWarpSchedules.get(whichAccess).add(Step.newShmToMem(memory.get(i), i));
			whichAccess = (whichAccess + 1) % accessWarps;
		}
	}
	
	void loadVar(Var v, boolean out, ArrayList<Var> deps) {
		int origin = scanPointer;
		for (scanPointer = origin; scanPointer != (origin + memory.size() - 1) % memory.size(); scanPointer = (scanPointer + 1) % memory.size()) {
			if (!loadable[scanPointer]) continue;
			if (memory.get(scanPointer) == null) {
				memory.set(scanPointer, v);
				used[scanPointer] = true;
				loadable[scanPointer] = false;
				if (!out) {
					accessWarpSchedules.get(whichAccess).add(Step.newMemToShm(v, scanPointer));
					whichAccess = (whichAccess + 1) % accessWarps;
				}
				return;
			}
		}
		for (scanPointer = origin; scanPointer != (origin + memory.size() - 1) % memory.size(); scanPointer = (scanPointer + 1) % memory.size()) {
			if (!loadable[scanPointer]) continue;
			if ((memory.get(scanPointer).type == Type.VARPRE || memory.get(scanPointer).type == Type.CONST) && !deps.contains(memory.get(scanPointer))) {
				memory.set(scanPointer, v);
				used[scanPointer] = true;
				loadable[scanPointer] = false;
				if (!out) {
					accessWarpSchedules.get(whichAccess).add(Step.newMemToShm(v, scanPointer));
					whichAccess = (whichAccess + 1) % accessWarps;
				}
				return;
			}
		}
		for (scanPointer = origin; scanPointer != (origin + memory.size() - 1) % memory.size(); scanPointer = (scanPointer + 1) % memory.size()) {
			if (!loadable[scanPointer]) continue;
			if (deps.contains(memory.get(scanPointer))) continue;
			accessWarpSchedules.get(whichAccess).add(Step.newShmToMem(memory.get(scanPointer), scanPointer));
			memory.set(scanPointer, v);
			used[scanPointer] = true;
			loadable[scanPointer] = false;
			if (!out) {
				accessWarpSchedules.get(whichAccess).add(Step.newMemToShm(v, scanPointer));
				whichAccess = (whichAccess + 1) % accessWarps;
			}
			return;
		}
		System.out.println("Errrrorrrrrr: " + v);
	}
	
	void useVar(Var v) {
		int where = memory.indexOf(v);
		if (where == -1) System.out.println("Error, not allocated: " + v);
		used[where] = true;
		loadable[where] = false;
		if (varLiveness.containsKey(v)) {
			int count = varLiveness.get(v);
			if (count > 1) varLiveness.put(v, count - 1);
			else {
				varLiveness.remove(v);
				memory.set(where, null);
			}
		}
		
	}
	
	private void scheduleEqu(Equation equ) {
		for (Var v : equ.getUsedDependencies()) {
			if (!memory.contains(v)) {
				loadVar(v, false, equ.getUsedDependencies());
			}
		}
		if (equ.getResultVar() != null) {
			loadVar(equ.getResultVar(), true, equ.getUsedDependencies());
		}
		computeWarpSchedules.get(whichCompute).add(Step.newComputation(equ, memory));
		whichCompute = (whichCompute + 1) % computeWarps;
		for (Var v : equ.getUsedDependencies()) {
			useVar(v);
		}
	}
	

}
