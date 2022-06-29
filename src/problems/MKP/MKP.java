package problems.MKP;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.io.StreamTokenizer;
import java.util.Arrays;
import java.util.TreeMap;

import problems.Evaluator;
import solutions.Solution;

public class MKP implements Evaluator<Integer> {

	/**
	 * Dimension of the domain.
	 */
	private final Integer size;

	/**
	 * Dimension of the knapsack.
	 */
	private Integer numKnapsacks;

	/**
	 * Optimal solution.
	 */
	private Double optimalValue;

	/**
	 * The array of numbers representing the domain.
	 */
	public final Double[] variables;

	/**
	 * The matrix r of coefficients for the MKP (Weights)
	 */
	private Double[][] r;

	/**
	 * The array p of coefficients for the MKP (Values)
	 */
	private Double[] p;


	/**
	 * The array b of constraints for the MKP (Knapsack limit)
	 */
	private Double[] b;

	/**
	 * The array of surrogate constraint coefficients
	 */
	private Double[] w;

	/**
	 * The treemap of pseudo utility ratios -> j 
	 */
	private TreeMap<Double, Integer> u;
	
	public Integer getSize() {
		return size;
	}

	public Integer getNumKnapsacks() {
		return numKnapsacks;
	}

	public Double getOptimalValue() {
		return optimalValue;
	}

	public Double[][] getR() {
		return r;
	}

	public Double[] getP() {
		return p;
	}

	public Double[] getB() {
		return b;
	}

	public Double[] getW() {
		return w;
	}

	public TreeMap<Double, Integer> getU() {
		return u;
	}

	public MKP(String filename) throws IOException {
		size = readInput(filename);
		variables = allocateVariables();
		u = calculatePseudoUtilityRatio();
	}

	public void setVariables(Solution<Integer> sol) {

		resetVariables();
		if (!sol.isEmpty()) {
			for (Integer elem : sol) {
				variables[elem] = 1.0;
			}
		}

	}

	@Override
	public Integer getDomainSize() {
		return size;
	}

	@Override
	public Double evaluate(Solution<Integer> sol) {

		setVariables(sol);
		return sol.cost = evaluateMKP();

	}

	public Double evaluateMKP() {

		Double sum = (double) 0;

		for (int j = 0; j < size; j++) {
			sum += variables[j] * p[j];
		}
		
		return sum;
	}

	@Override
	public Double evaluateInsertionCost(Integer elem, Solution<Integer> sol) {

		setVariables(sol);
		return evaluateInsertionMKP(elem);

	}

	public Double evaluateInsertionMKP(int i) {

		if (variables[i] == 1)
			return 0.0;

		return evaluateContributionMKP(i);
	}

	@Override
	public Double evaluateRemovalCost(Integer elem, Solution<Integer> sol) {

		setVariables(sol);
		return evaluateRemovalMKP(elem);

	}
	public Double evaluateRemovalMKP(int i) {

		if (variables[i] == 0)
			return 0.0;

		return -evaluateContributionMKP(i);

	}

	@Override
	public Double evaluateExchangeCost(Integer elemIn, Integer elemOut, Solution<Integer> sol) {

		setVariables(sol);
		return evaluateExchangeMKP(elemIn, elemOut);

	}

	public Double evaluateExchangeMKP(int in, int out) {

		Double sum = 0.0;

		if (in == out)
			return 0.0;
		if (variables[in] == 1)
			return evaluateRemovalMKP(out);
		if (variables[out] == 0)
			return evaluateInsertionMKP(in);

		sum += evaluateContributionMKP(in);
		sum -= evaluateContributionMKP(out);

		return sum;
	}

	private Double evaluateContributionMKP(int j) {
		return p[j];
	}

	protected Integer readInput(String filename) throws IOException {

		Reader fileInst = new BufferedReader(new FileReader(filename));
		StreamTokenizer stok = new StreamTokenizer(fileInst);

		stok.nextToken();
		Integer _size = (int) stok.nval;
		stok.nextToken();
		numKnapsacks = (int) stok.nval;
		stok.nextToken();
		optimalValue = stok.nval;

		r = new Double[numKnapsacks][_size];
		p = new Double[_size];
		b = new Double[numKnapsacks];
		w = new Double[numKnapsacks];

		for (int j = 0; j < _size; j++) {
			stok.nextToken();
			p[j] = stok.nval;
		}

		for (int i = 0; i < numKnapsacks; i++) {
			for (int j = 0; j < _size; j++) {
				stok.nextToken();
				r[i][j] = stok.nval;
			}
		}
		
		for (int i = 0; i < numKnapsacks; i++) {
			stok.nextToken();
			b[i] = stok.nval;
		}

		for (int i = 0; i < numKnapsacks; i++) {
			stok.nextToken();
			w[i] = stok.nval;
		}

		return _size;
	}

	protected Double[] allocateVariables() {
		Double[] _variables = new Double[size];
		return _variables;
	}

	public void resetVariables() {
		Arrays.fill(variables, 0.0);
	}

	// Return pseudo utility ratios in ascending order
	private TreeMap<Double, Integer> calculatePseudoUtilityRatio() {
		
		TreeMap<Double, Integer> u = new TreeMap<>();

		for (int j=0; j < size; j++) {
			Double sum = 0.0;
			for (int i=0; i < numKnapsacks; i++) {
				sum += r[i][j];
			}
			u.put( p[j]/sum, j);
		}
		return u;
	}
	public static void main(String[] args) throws IOException {

		MKP mkp = new MKP("instances/mkp0.txt");
		Double maxVal = Double.NEGATIVE_INFINITY;
		
		// evaluates randomly generated values for the domain, saving the best
		// one.
		for (int i = 0; i < 10000000; i++) {
			for (int j = 0; j < mkp.size; j++) {
				if (Math.random() < 0.5)
					mkp.variables[j] = 0.0;
				else
					mkp.variables[j] = 1.0;
			}
			//System.out.println("x = " + Arrays.toString(qbf.variables));
			Double eval = mkp.evaluateMKP();
			//System.out.println("f(x) = " + eval);
			if (maxVal < eval)
				maxVal = eval;
		}
		System.out.println("maxVal = " + maxVal);

		// evaluates the zero array.
		for (int j = 0; j < mkp.size; j++) {
			mkp.variables[j] = 0.0;
		}
		System.out.println("x = " + Arrays.toString(mkp.variables));
		System.out.println("f(x) = " + mkp.evaluateMKP());

		// evaluates the all-ones array.
		for (int j = 0; j < mkp.size; j++) {
			mkp.variables[j] = 1.0;
		}
		System.out.println("x = " + Arrays.toString(mkp.variables));
		System.out.println("f(x) = " + mkp.evaluateMKP());

	}

}
