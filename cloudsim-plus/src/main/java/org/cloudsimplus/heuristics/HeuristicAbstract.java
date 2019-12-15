/*
 * CloudSim Plus: A modern, highly-extensible and easier-to-use Framework for
 * Modeling and Simulation of Cloud Computing Infrastructures and Services.
 * http://cloudsimplus.org
 *
 *     Copyright (C) 2015-2018 Universidade da Beira Interior (UBI, Portugal) and
 *     the Instituto Federal de Educação Ciência e Tecnologia do Tocantins (IFTO, Brazil).
 *
 *     This file is part of CloudSim Plus.
 *
 *     CloudSim Plus is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     CloudSim Plus is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with CloudSim Plus. If not, see <http://www.gnu.org/licenses/>.
 */
package org.cloudsimplus.heuristics;

import org.cloudbus.cloudsim.cloudlets.Cloudlet;
import org.cloudbus.cloudsim.distributions.ContinuousDistribution;
import org.cloudbus.cloudsim.vms.Vm;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.*;

/**
 * A base class for {@link Heuristic} implementations.
 *
 * @author Manoel Campos da Silva Filho
 * @param <S> The {@link HeuristicSolution class of solutions} the heuristic will deal with.
 *            It start with an initial
 *           solution (usually random, depending on each sub-class implementation)
 *           and executes the solution search in order
 *           to find a satisfying solution (defined by a stop criteria)
 * @since CloudSim Plus 1.0
 */
public abstract class HeuristicAbstract<S extends HeuristicSolution<?>>  implements Heuristic<S> {
	/**
	 * Reference to the generic class that will be used to instantiate objects.
	 */
	public static Double M;
    private final Class<S> solutionClass;

	private final ContinuousDistribution random;
	/**
	 * @see #getNeighborhoodSearchesByIteration()
	 */
    private int neighborhoodSearchesByIteration;
	/**
	 * @see #getBestSolutionSoFar()
	 */
    private S bestSolutionSoFar;
	/**
	 * @see #getNeighborSolution()
	 */
    private S neighborSolution;

	/**
	 * @see #getSolveTime()
	 */
	private double solveTime;

	/**
	 * Creates a heuristic.
	 *
	 * @param random a random number generator
	 * @param solutionClass reference to the generic class that will be used to instantiate heuristic solutions
	 */
	/* default */ HeuristicAbstract(final ContinuousDistribution random, final Class<S> solutionClass){
		this.solutionClass = solutionClass;
		this.random = random;
		this.neighborhoodSearchesByIteration = 1;
		setBestSolutionSoFar(newSolutionInstance());
		setNeighborSolution(bestSolutionSoFar);
	}

	@Override
	public double getSolveTime() {
		return solveTime;
	}

	/**
	 * Sets the time taken to solve the heuristic.
	 * @param solveTime the time to set (in seconds)
	 */
	protected void setSolveTime(final double solveTime) {
		this.solveTime = solveTime;
	}

	/**
	 *
	 * @return a random number generator
	 */
	protected ContinuousDistribution getRandom(){
		return random;
	}

	private S newSolutionInstance() {
	    try {
	        final Constructor<S> constructor = solutionClass.getConstructor(Heuristic.class);
	        return constructor.newInstance(this);
	    } catch (IllegalArgumentException | InvocationTargetException | NoSuchMethodException | SecurityException | InstantiationException | IllegalAccessException ex) {
	        throw new RuntimeException(ex);
	    }
	}

	/**
	 * Updates the state of the system in order to keep looking
	 * for a suboptimal solution.
	 */
	protected abstract void updateSystemState();

	@Override
	public int getRandomValue(final int maxValue){
		final double uniform = getRandom().sample();

        /*always get an index between [0 and size[,
        regardless if the random number generator returns
        values between [0 and 1[ or >= 1*/
		return (int)(uniform >= 1 ? uniform % maxValue : uniform * maxValue);
	}

	@Override
	public S solve() {
		final long startTime = System.currentTimeMillis();
//		setBestSolutionSoFar(getInitialSolution());
        List<Vm> vmList =  getVmList();
        List<Cloudlet> cloudletList = getCloudletList();
        //TODO Calculate M

        // defines initial population
        List<List<Vm>> population = getInitPop(vmList);

        // fitness calculation
        List<Double> fitVals = getFitnessValues(population);

        // finding the leader route
        List<Vm> leader = null;
        Double maxFitness = Double.MIN_VALUE;
        for(int i = 0; i < population.size(); i++) {
            if(fitVals.get(i) > maxFitness) {
                leader = population.get(i);
                maxFitness = fitVals.get(i);
            }
        }

        // iteration start here
        int iteration = 0;
        int maxIteration = 100;
		while (iteration < maxIteration) {
		    // new population that would replace the current population at the end of the iteration
            List<List<Vm>> newPopulation = new ArrayList<>();

            // for each search agent
            for (List<Vm> agent : population){
                // TODO: evaluate values
                // select k
                // update route
            }

            // update system state : (i) Replace population with newly generated population (ii) Recalculate fitVals, (iii) Update leader route, (iv) iteration++
            population = newPopulation;
            fitVals = getFitnessValues(population, cloudletList);
            maxFitness = Double.MIN_VALUE;
            for (int i = 0; i < population.size(); i++) {
                if (fitVals.get(i) > maxFitness) {
                    leader = population.get(i);
                    maxFitness = fitVals.get(i);
                }
            }
            iteration++;
        }

        // get S : cloudletToVmMapping
		S clVmMap = convertToSolution(leader, cloudletList);

        setSolveTime((System.currentTimeMillis() - startTime)/1000.0);

        return clVmMap;
	}

	public S convertToSolution(List<Vm> leader, List<Cloudlet> cloudletList){
        final CloudletToVmMappingSolution solution = new CloudletToVmMappingSolution(this);
        int i=0, j=0;
        long load = 0;
        while(i<cloudletList.size() && j<leader.size()){
            Cloudlet cl = cloudletList.get(i++);
            load+=cl.getLength();
            solution.bindCloudletToVm(cl, leader.get(j));

            if(load>M){
                j++;
                load=0;
            }
        }

        //TODO not sure about this typecase
        return (S)solution;
    }
	public List<List<Vm>> getInitPop(List<Vm> vmList){
	    // TODO return the population by permutation: 30
        return new ArrayList<>();
    }

    public Double getFitnessValue(List<Vm> perm, List<Cloudlet> cloudlets){
        int i=0, j=0;
        long[] loads = new long[perm.size()];
        while(i < cloudlets.size() && j<perm.size()){
            loads[j]+=cloudlets.get(i++).getLength();
            if(loads[j]>M){
                j++;
            }
        }
        long maxLoad = -1;
        for(long load: loads)
            maxLoad = Math.max(maxLoad,load);

        double avgUtilization = 0;
        for(long load: loads){
            avgUtilization+=(load/maxLoad);
        }
        avgUtilization/=perm.size();

        //assuming all processor queues are valid (just to get this to work, threshold can be added later)

        return avgUtilization/(double)maxLoad;

    }

    public List<Double> getFitnessValues(List<List<Vm>> population, List<Cloudlet> cloudlets){
	    List<Double> rv = new ArrayList<>();
	    for(List<Vm> vms: population){
	        rv.add(getFitnessValue(vms, cloudlets));
        }
	    return rv;
    }

    private void searchSolutionInNeighborhood() {
        for (int i = 0; i < getNeighborhoodSearchesByIteration(); i++) {
            setNeighborSolution(createNeighbor(getBestSolutionSoFar()));
            if (getAcceptanceProbability() > getRandomValue(1)) {
                setBestSolutionSoFar(getNeighborSolution());
            }
        }
    }

    private void updateSystemStates(S solution) {

    }
//    public
    @Override
	public S getBestSolutionSoFar() {
	    return bestSolutionSoFar;
	}

	@Override
	public S getNeighborSolution() {
	    return neighborSolution;
	}

	/**
	 * Sets a solution as the current one.
	 * @param solution the solution to set as the current one.
	 */
	protected final void setBestSolutionSoFar(final S solution) {
        this.bestSolutionSoFar = solution;
    }

	/**
	 * Sets a solution as the neighbor one.
	 * @param neighborSolution the solution to set as the neighbor one.
	 */
    protected final void setNeighborSolution(final S neighborSolution) {
        this.neighborSolution = neighborSolution;
    }

    /**
     * Gets the number of neighborhood searches by each iteration of the heuristic.
     * @return
     */
	public int getNeighborhoodSearchesByIteration() {
        return neighborhoodSearchesByIteration;
    }

    /**
     * Sets the number of neighborhood searches by each iteration of the heuristic.
     * @param neighborhoodSearches the number of neighborhood searches to set
     */
	public void setNeighborhoodSearchesByIteration(final int neighborhoodSearches) {
        this.neighborhoodSearchesByIteration = neighborhoodSearches;
    }

    public double getVectorLength(List<Double> vec){
	    double val = 0;
	    for(double d: vec)
	        val+=(d*d);
	    return Math.sqrt(val);
    }
    public List<Double> getA(List<Double> a, List<Double> r){
        List<Double> rv = new ArrayList<>();
        for(int i=0;i<a.size();i++){
            rv.add(2*a.get(i)*r.get(i) - a.get(i));
        }
        return rv;
    }
    public List<Double> geta(int curIt, int totIt, List<Vm> vms){
        double val = 2 - 2*((double)(curIt))/totIt;
        List<Double> rv = new ArrayList<>();
        for(Vm vm: vms)
            rv.add(val);

        return rv;
    }
    public List<Double> getC(List<Double> r){
	    List<Double> rv = new ArrayList<>();
	    for(double val: r){
	        rv.add(2*val);
        }
        return rv;
    }
    public List<Double> getr(List<Vm> vms){
	    List<Double> rVec = new ArrayList<>();
	    for(Vm vm: vms){
	        rVec.add(getRandomInRange(0,1));
        }
        return rVec;
    }
    public double getRandomInRange(int low, int high){
        return low + Math.random()*(high-low);
    }
}
