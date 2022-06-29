package solutions;

import java.util.ArrayList;

@SuppressWarnings("serial")
public class Solution<E> extends ArrayList<E> {
	
	public Double cost = Double.POSITIVE_INFINITY;
	
	public Solution() {
		super();
	}
	
	public Solution(Solution<E> sol) {
		super(sol);
		cost = sol.cost;
	}

	@Override
	public String toString() {
		return "Solution: cost=[" + cost + "], size=[" + this.size() + "], elements=" + super.toString();
	}

}

