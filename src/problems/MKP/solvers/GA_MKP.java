package problems.MKP.solvers;

import java.io.IOException;
import java.util.Map;

import metaheuristics.ga.AbstractGA;
import problems.MKP.MKP;
import solutions.Solution;

public class GA_MKP extends AbstractGA<Integer, Integer> {

	public GA_MKP(Integer generations, Integer popSize, Double mutationRate, String filename) throws IOException {
		super(new MKP(filename), generations, popSize, mutationRate);
	}

	@Override
	public Solution<Integer> createEmptySol() {
		Solution<Integer> sol = new Solution<Integer>();
		sol.cost = 0.0;
		return sol;
	}

	@Override
	protected Solution<Integer> decode(Chromosome chromosome) {

		Solution<Integer> solution = createEmptySol();
		for (int locus = 0; locus < chromosome.size(); locus++) {
			if (chromosome.get(locus) == 1) {
				solution.add(new Integer(locus));
			}
		}

		ObjFunction.evaluate(solution);
		return solution;
	}

	@Override
	protected Chromosome generateRandomChromosome() {

		Chromosome chromosome = new Chromosome();
		for (int i = 0; i < chromosomeSize; i++) {
			chromosome.add(rng.nextInt(2));
		}
		repairChromosome(chromosome);
		return chromosome;
	}

	@Override
	protected Double fitness(Chromosome chromosome) {

		return decode(chromosome).cost;

	}

	@Override
	protected void mutateGene(Chromosome chromosome, Integer locus) {

		chromosome.set(locus, 1 - chromosome.get(locus));
		repairChromosome(chromosome);

	}

	@Override
	protected void repairChromosome(Chromosome chromosome) {
		MKP mkp = (MKP) this.ObjFunction;
		Double[] accumulatedResources = new Double[mkp.getNumKnapsacks()];
		for (int i=0; i < mkp.getNumKnapsacks(); i++) {
			Double sum = 0.0;
			for (int j=0; j < mkp.getDomainSize(); j++) {
				sum += mkp.getR()[i][j]*chromosome.get(j);
			}
			accumulatedResources[i] = sum;
		}

		if (isValidChromosome(accumulatedResources)) {
			return;
		}

		// DROP Phase
		for (Map.Entry<Double, Integer> entry : mkp.getU().entrySet()) {
			int j = entry.getValue();
			if ((chromosome.get(j) == 1) && (!isValidChromosome(accumulatedResources))) {
				chromosome.set(j, 0);
				updateAccumulatedResourcesDrop(accumulatedResources, j);
			}
		}

		// ADD Phase
		for (Map.Entry<Double, Integer> entry : mkp.getU().descendingMap().entrySet()) {
			int j = entry.getValue();
			if ((chromosome.get(j) == 0) && (isValidChromosome(accumulatedResources, j))) {
				chromosome.set(j, 1);
				updateAccumulatedResourcesAdd(accumulatedResources, j);
			}
		}
	}

	private void updateAccumulatedResourcesDrop(Double[] accumulatedResources, int j) {
		MKP mkp = (MKP) this.ObjFunction;
		for (int i=0; i < mkp.getNumKnapsacks(); i++) {
			accumulatedResources[i] -= mkp.getR()[i][j];
		}
	}

	private void updateAccumulatedResourcesAdd(Double[] accumulatedResources, int j) {
		MKP mkp = (MKP) this.ObjFunction;
		for (int i=0; i < mkp.getNumKnapsacks(); i++) {
			accumulatedResources[i] += mkp.getR()[i][j];
		}
	}

	private boolean isValidChromosome(Double[] accumulatedResources) {
		MKP mkp = (MKP) this.ObjFunction;
		for (int i=0; i < mkp.getNumKnapsacks(); i++) {
			if (accumulatedResources[i] > mkp.getB()[i]) {
				return false;
			}
		}
		return true;
	}

	private boolean isValidChromosome(Double[] accumulatedResources, int j) {
		MKP mkp = (MKP) this.ObjFunction;
		for (int i=0; i < mkp.getNumKnapsacks(); i++) {
			if (accumulatedResources[i] + mkp.getR()[i][j] > mkp.getB()[i]) {
				return false;
			}
		}
		return true;
	}

	public static void main(String[] args) throws IOException {

		long startTime = System.currentTimeMillis();
		GA_MKP ga = new GA_MKP(1000, 100, 1.0 / 100.0, "instances/mkp42.txt");
		Solution<Integer> bestSol = ga.solve();
		System.out.println("maxVal = " + bestSol);
		System.out.println(("Optimal solution: " + ((MKP) ga.ObjFunction).getOptimalValue()));
		long endTime = System.currentTimeMillis();
		long totalTime = endTime - startTime;
		System.out.println("Time = " + (double) totalTime / (double) 1000 + " seg");

	}

}
