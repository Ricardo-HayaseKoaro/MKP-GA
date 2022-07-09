package problems.MKP.solvers;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import apple.laf.JRSUIUtils.InternalFrame;
import metaheuristics.ga.AbstractGA;
import problems.MKP.MKP;
import solutions.Solution;

public class GA_MKP extends AbstractGA<Integer, Integer> {

	public enum MutationStrategy {
		STANDARD, ADAPTIVE
	}

	public enum CrossoverOperator {
		MULTIPOINT, UNIFORM
	}

	public enum RepairFunction {
		DROPADD, RANDOM
	}


	private final MutationStrategy mutationStrategy;
	private final CrossoverOperator crossoverOperator;
	private final RepairFunction repairFunction;

	public GA_MKP(Integer generations, Integer popSize, Double mutationRate, String filename, Boolean verbose, MutationStrategy mutationStrategy, CrossoverOperator crossoverOperator, RepairFunction repairFunction, Boolean localSearchEnabled) throws IOException {
		super(new MKP(filename), generations, popSize, mutationRate, verbose, localSearchEnabled);
		this.mutationStrategy = mutationStrategy;
		this.crossoverOperator = crossoverOperator;
		this.repairFunction = repairFunction;
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
	protected Population mutate(Population offsprings) {
		double diversity = 0.0;
		if (mutationStrategy.equals(MutationStrategy.ADAPTIVE)) {
			Double worseSolCost = decode(getWorseChromosome(offsprings)).cost;
			worseSolCost += Math.abs(worseSolCost) + 1;
			Double bestSolCost = decode(getBestChromosome(offsprings)).cost + Math.abs(worseSolCost) + 1;
			diversity = 1.0 - (worseSolCost / bestSolCost);
		}

		for (Chromosome c : offsprings) {
			for (int locus = 0; locus < chromosomeSize; locus++) {
				if (applyMutationStrategy(diversity)) {
					mutateGene(c, locus);
				}
			}
		}

		return offsprings;
	}

	private Boolean applyMutationStrategy(Double diversity) {
		switch (mutationStrategy) {
			case STANDARD:
				return rng.nextDouble() < mutationRate;
			case ADAPTIVE:
				// When difference between worst and best cost is higher than 0.40, doubled the mutationRate
				double rateBasedOnDiversity = diversity < 0.40 ? mutationRate * 2 : mutationRate;
				return rng.nextDouble() < rateBasedOnDiversity;
			default:
				return false;
		}
	}

	protected void repairChromosome(Chromosome chromosome) {
		switch (repairFunction) {
			case RANDOM:
				repairChromosomeRandom(chromosome);
			default:
				repairChromosomeDropAdd(chromosome);
		}
	}

	protected void repairChromosomeDropAdd(Chromosome chromosome) {
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

	protected void repairChromosomeRandom(Chromosome chromosome) {
		MKP mkp = (MKP) this.ObjFunction;
		Double[] accumulatedResources = new Double[mkp.getNumKnapsacks()];
		for (int i=0; i < mkp.getNumKnapsacks(); i++) {
			Double sum = 0.0;
			for (int j=0; j < mkp.getDomainSize(); j++) {
				sum += mkp.getR()[i][j]*chromosome.get(j);
			}
			accumulatedResources[i] = sum;
		}

		int randomIdx = rng.nextInt(chromosome.size());		
		while (!isValidChromosome(accumulatedResources)) {
			while (chromosome.get(randomIdx) != 1) {
				randomIdx = rng.nextInt(chromosome.size());
			}
			chromosome.set(randomIdx, 0);
			updateAccumulatedResourcesDrop(accumulatedResources, randomIdx);
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

	private boolean isValidChromosome(Double[] accumulatedResources, int candIn, int candOut) {
		MKP mkp = (MKP) this.ObjFunction;
		for (int i=0; i < mkp.getNumKnapsacks(); i++) {
			if (accumulatedResources[i] + mkp.getR()[i][candIn] - mkp.getR()[i][candOut] > mkp.getB()[i]) {
				return false;
			}
		}
		return true;
	}


	@Override
	protected Population crossover(Population parents) {

		Population offsprings = new Population();

		for (int i = 0; i < popSize; i = i + 2) {

			Chromosome parent1 = parents.get(i);
			Chromosome parent2 = parents.get(i + 1);

			int crosspoint1 = rng.nextInt(chromosomeSize + 1);
			int crosspoint2 = crosspoint1 + rng.nextInt((chromosomeSize + 1) - crosspoint1);

			Chromosome offspring1 = new Chromosome();
			Chromosome offspring2 = new Chromosome();

			for (int j = 0; j < chromosomeSize; j++) {
				if(crossoverOperator.equals(CrossoverOperator.MULTIPOINT)) {
					if (j >= crosspoint1 && j < crosspoint2) {
						offspring1.add(parent2.get(j));
						offspring2.add(parent1.get(j));
					} else {
						offspring1.add(parent1.get(j));
						offspring2.add(parent2.get(j));
					}
				} else {
					if(rng.nextDouble() > 0.5) {
						offspring1.add(parent1.get(j));
						offspring2.add(parent2.get(j));
					} else {
						offspring1.add(parent2.get(j));
						offspring2.add(parent1.get(j));
					}
				}
			}

			repairChromosome(offspring1);
			repairChromosome(offspring2);

			offsprings.add(offspring1);
			offsprings.add(offspring2);

		}

		return offsprings;

	}

	public void localSearch(Chromosome chromosome) {
		MKP mkp = (MKP) this.ObjFunction;
		Double[] accumulatedResources = new Double[mkp.getNumKnapsacks()];
		for (int i=0; i < mkp.getNumKnapsacks(); i++) {
			Double sum = 0.0;
			for (int j=0; j < mkp.getDomainSize(); j++) {
				sum += mkp.getR()[i][j]*chromosome.get(j);
			}
			accumulatedResources[i] = sum;
		}

		Double maxDeltaCost;
		Integer bestCandIn = null, bestCandOut = null;

		do {
			maxDeltaCost = 0.0;
			Solution<Integer> sol = decode(chromosome);

			// Evaluate insertions
			for (int j=0; j < mkp.getDomainSize(); j++ ) {
				Integer candIn = chromosome.get(j);
				if (candIn == 0) {
					double deltaCost = ObjFunction.evaluateInsertionCost(j, sol);
					if (deltaCost > maxDeltaCost && isValidChromosome(accumulatedResources, j)) {
						maxDeltaCost = deltaCost;
						bestCandIn = j;
						bestCandOut = null;
					}
				}
			}
		
			// Evaluate exchanges
			for (int j=0; j < mkp.getDomainSize(); j++ ) {
				Integer candIn = chromosome.get(j);
				if (candIn == 0) {
					for (int j2=0; j2 < mkp.getDomainSize(); j2++ ) {
						Integer candOut = chromosome.get(j2);
						if (candOut == 1) {
							double deltaCost = ObjFunction.evaluateExchangeCost(j, j2, sol);
							if (j != j2 && deltaCost > maxDeltaCost && isValidChromosome(accumulatedResources, j, j2)) {
								maxDeltaCost = deltaCost;
								bestCandIn = j;
								bestCandOut = j2;
							}
						}
						
					}
				}
			}
			
			// Implement the best move, if it reduces the solution cost.
			if (maxDeltaCost > 0.0) {
				if (bestCandOut != null) {
					chromosome.set(bestCandOut, 0);
				}
				if (bestCandIn != null) {
					chromosome.set(bestCandIn, 1);
				}
				sol = decode(chromosome);
			}

		} while (maxDeltaCost != 0.0);

	}

	public static void main(String[] args) throws IOException {

		int generations = 1000;
		int popSize = 100;
		Double mutationRate = 0.001;
		Boolean csv = false;
		CrossoverOperator crossoverOperator = CrossoverOperator.UNIFORM;
		MutationStrategy mutationStrategy = MutationStrategy.ADAPTIVE;
		RepairFunction repairFunction = RepairFunction.DROPADD;
		Boolean localSearchEnabled = true;
		String solution_name = "solution4.csv";

		//solution1 = repair method Dro

		List<File> files= Arrays.asList(new File("/Users/ricardohk/Dev/GA-Framework/MKP-GA/instances/").listFiles());
		Collections.sort(files, (o1, o2) -> (((File) o1).getName().compareTo(((File) o2).getName())));

		List<String> lines = new ArrayList<>();
		lines.add("Filename, Time(ms), Solution, Optimal solution");

		for (File file : files) {
			if (file.getName().contains("DS")){
				continue;
			}

			long startTime = System.currentTimeMillis();
			GA_MKP ga = new GA_MKP(Integer.valueOf(generations), Integer.valueOf(popSize), mutationRate, file.getAbsolutePath(), !csv, mutationStrategy, crossoverOperator, repairFunction, localSearchEnabled);
			Solution<Integer> bestSol = ga.solve();
			if (csv) {
				System.out.println("filename = " + file.getName());
				long endTime = System.currentTimeMillis();
				long totalTime = endTime - startTime;
				lines.add(file.getName() + "," + totalTime + "," + bestSol.cost + "," + ((MKP) ga.ObjFunction).getOptimalValue().toString());
			} else {
				System.out.println("maxVal = " + bestSol);
				System.out.println(("Optimal solution: " + ((MKP) ga.ObjFunction).getOptimalValue()));
				long endTime = System.currentTimeMillis();
				long totalTime = endTime - startTime;
				System.out.println("Time = " + (double) totalTime / (double) 1000 + " seg");
			}
		}
		if (csv) {
			Path solution = Paths.get("/Users/ricardohk/Dev/GA-Framework/solutions/" + solution_name);
			Files.write(solution, lines, StandardCharsets.UTF_8);
		}
	}
}
