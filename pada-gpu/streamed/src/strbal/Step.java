package strbal;

import java.util.HashMap;
import java.util.Vector;

public class Step {
	enum StepType {SYNCPOINT, COMPUTE, MEMTOSHM, SHMTOMEM};
	private StepType type;
	private Equation equ = null;
	private int shmLoc;
	private Var v;
	Vector<Var> varBinding;
	
	public String toString() {
		switch (type) {
			case SYNCPOINT:
				return "\t\t__syncthreads();__threadfence();";
			case COMPUTE:
				return "\t\t" + equ.bindForm(varBinding) + ";";
			case SHMTOMEM:
				return "\t\t" + v.globalLoc() + " = local_data[" + shmLoc + "];";
			case MEMTOSHM:
				return "\t\tlocal_data[" + shmLoc + "] = " + v.globalLoc() + ";"; 
			default: 
				return "";
		}
	}
	
	public boolean isSyncPoint() {
		return type == StepType.SYNCPOINT;
	}
	
	private Step(StepType type) {
		// creates a sync. point in the schedule
		assert (type == StepType.SYNCPOINT);
		this.type = type;
	}
	
	private Step(StepType type, Equation e, Vector<Var> mem_map) {
		// adds an equation to the schedule
		assert (type == StepType.COMPUTE);
		equ = e;
		varBinding = new Vector<Var>(mem_map);
		this.type = type;
	}

	private Step(StepType type, Var v, int loc) {
		assert (type == StepType.MEMTOSHM|| type == StepType.SHMTOMEM);
		shmLoc = loc;
		this.v = v;
		this.type = type;
	}
	
	public static Step newSyncPoint() {
		return new Step(StepType.SYNCPOINT);
	}
	
	public static Step newMemToShm(Var v, int scanPointer) {
		return new Step(StepType.MEMTOSHM, v, scanPointer);
	}

	public static Step newShmToMem(Var v, int scanPointer) {
		return new Step(StepType.SHMTOMEM, v, scanPointer);
	}

	public static Step newComputation(Equation equ, Vector<Var> memory) {
		return new Step(StepType.COMPUTE, equ, memory);
	}

}
