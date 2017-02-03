package sim;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.*;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;

import org.apache.commons.cli.*;

import sun.reflect.Reflection;
import utility.MyFunc;
import utility.MyMat;

public class SimulatorOld {

	// parameters for simulations
	double mutationRate = 0.0001;
	double growthRate = 0.0001;
	double deathRate =  0.0000001;
	double deathRateForNonStem =  0.01;
	double symmetricReplicationProbablity = 0.1;
	int gatewayDriver = 0;
	int initialPopulationSize = 10;
	int maxTime = 10000000;
	int maxPopulationSize = 1000000;
	int genomeSize = 30;
	int  driverSize = 10;
	double fitnessIncrease  = 5;
	
	//simulated tumor 
	List <Cell> population;
	
	//for simulation
	Random R = new Random();
	
	//parameters for analyzing results
	String header = "time\tpopulationSize\tentropy\tmutationCount\tmutatedDriverGeneCount\taverageFitness\tfounderMutationPropotion\tmajorSubpopulation";
	int mutationProfileMatrixSize = 1000;
	double subpopulationProportionCutoff = 0.01;
	double populationCoverageForfounderMutation = 0.99;
	
	//statistics for analyzing results 
	int time;
	int populationSize;
	double mutatedGeneCount;
	double mutatedDriverGeneCount;
	double founderMutationPropotion;
	double averageFitness;
	Map<String, Double> subpopulationProportion;
	Map<String, Double> majorSubpopulationProportion;
	Double populationEntropy;
	List <Double> eigenValue;
	
	
	//for simulation
	
	
	public class Genome {
		Set <Integer> mutatedGenes;
		
		Genome(){
			 mutatedGenes = new TreeSet<Integer>();
		}
		Genome(Genome genome){
			 mutatedGenes = new TreeSet<Integer>(genome.mutatedGenes);
		}
		
		boolean isMutated(int i){
			return mutatedGenes.contains(i);
		}
		
		void mutate(){
			for(int i=0; i < genomeSize; i++){
				if(!isMutated(i) & (R.nextDouble() < mutationRate)){
					mutatedGenes.add(i);
				}
			}
		}
	}
	
	
	public class Cell {
		Genome genome;
		double fitness;
		boolean isStem;
				
		Cell(){
			genome = new Genome();
			fitness = 1;
			isStem = true;
		}
			
		Cell(Cell C){
			fitness  = C.fitness;
			isStem = C.isStem;
			genome = new Genome(C.genome);
		}
		
		void mutateGenome(){
			genome.mutate();
			setPenotype();
		}
		void setPenotype(){
			fitness = 1;
			boolean gateway = true;
			if(gatewayDriver > 0){
		      for(int i=0; i< gatewayDriver; i++){	
				if(genome.isMutated(i)){
					 for(int j=0; j< i; j++){
						 if(!genome.isMutated(j)){
							 gateway = false;
						 }
					 }
				}
		      }
		      for(int i=gatewayDriver; i< driverSize; i++){
		    	  if(genome.isMutated(i)){
		    		  for(int j=0; j< gatewayDriver; j++){
							 if(!genome.isMutated(j)){
								 gateway = false;
							 }
		    		  }
		    	  }
		      }
			}
			if(!gateway){
				fitness = -1;
			}else{
				for(int i=0; i< driverSize; i++){
					if(genome.isMutated(i)){
						fitness *= fitnessIncrease;
					}
				}
			}		
		}
		
		String getGenomeString(){
			return genome.mutatedGenes.toString();
		}
		
	}
	
	
	public SimulatorOld(){}
	
	public SimulatorOld(SimulatorOld S){
		mutationRate = S.mutationRate;
		growthRate = S.growthRate;
		deathRate = S.deathRate;
		initialPopulationSize = S.initialPopulationSize;
		maxTime = S.maxTime;
		maxPopulationSize = S.maxPopulationSize;
		genomeSize = S.genomeSize;
		driverSize = S.driverSize;
		fitnessIncrease  = S.fitnessIncrease;
		subpopulationProportionCutoff = S.subpopulationProportionCutoff;
		deathRateForNonStem  = S.deathRateForNonStem;
		symmetricReplicationProbablity = S.symmetricReplicationProbablity;
	}
	
	
	void initialize(){
		population = new ArrayList<Cell>();
		for(int i = 0; i < initialPopulationSize; i++){
			population.add(new Cell());
		}
	}
	
	public void simulate(){
		initialize();
		for(time = 1; time <= maxTime; time++){
			growPopulation();
			//System.err.println(population.size());
			if(population.size() > maxPopulationSize){
				break;
			}
		}
	}
	
	
	public void simulateWhilePrintingStatistics(int n){
		initialize();
		System.out.println(header);
		for(time = 1; time <= maxTime; time++){
			growPopulation();
			if(Math.abs(((double)time)/n - Math.round(time/n)) < 0.0000000000001){
				getStatistics();
				System.out.println(this);
			}
			if(population.size() > maxPopulationSize){
				break;
			}
		}
		getStatistics();
		System.out.println(this);
	}
	
	void growPopulation(){
		int n = population.size();
		for(int i = 0; i < n & i < population.size(); i++){
			Cell C = population.get(i);
			if(R.nextDouble() < C.fitness*growthRate){
				C.mutateGenome();
				Cell newC = new Cell(C);
				if(C.isStem & R.nextDouble() > symmetricReplicationProbablity){
					newC.isStem = false;
				}
				population.add(newC);
			}
			if((C.isStem & R.nextDouble() < deathRate) | (!C.isStem & R.nextDouble() < deathRateForNonStem) | C.fitness < 0){
				population.remove(C);
				i--;
			}
		}
	}
	
	// for analysis of  simulation results
	
	void getStatistics(){
		populationSize = population.size();
		getSubpopulationProportion();
		getMajorSubpopulationProportion();
		getPopulationEntropy();
		getMutatedGeneCount();
		getMutatedGeneCount();
	}
	
	public String toString(){
		return time + "\t" + populationSize + "\t" + populationEntropy + "\t" + mutatedGeneCount  + "\t" + mutatedDriverGeneCount  + "\t" + averageFitness + "\t" + founderMutationPropotion + "\t" + majorSubpopulationProportion;
	}
	
	void  getSubpopulationProportion(){
		subpopulationProportion = new HashMap<String, Double>();
		for(Cell C:population){
			String k = C.getGenomeString();
			if(subpopulationProportion.containsKey(k)){
				subpopulationProportion.put(k, subpopulationProportion.get(k)+1.0);
			}else{
				subpopulationProportion.put(k, 1.0);
			}
		}
		for(String k: subpopulationProportion.keySet()){
			subpopulationProportion.put(k, subpopulationProportion.get(k)/populationSize);
		}
	}
	
	void getMajorSubpopulationProportion(){
		majorSubpopulationProportion = new LinkedHashMap<String, Double>();
		for(String k: MyFunc.sortKeysByDescendingOrderOfValues(subpopulationProportion)){
			if(subpopulationProportion.get(k) > subpopulationProportionCutoff){
				majorSubpopulationProportion.put(k, subpopulationProportion.get(k));
			}else{
				break;
			}
		}
	}
	
	/*void getPopulationEntropy(){
		populationEntropy = 0.0;
		for(String k: subpopulationProportion.keySet()){
			populationEntropy -= subpopulationProportion.get(k) * Math.log(subpopulationProportion.get(k));
		}
	}*/
	
	void getPopulationEntropy() {
		try {
			populationEntropy = getEigenValueEntropy(getMutaionSimilarityMatrix());
		} catch (NotConvergedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	
	MyMat getMutaionProfileMatrix(){
		MyMat M = new MyMat(genomeSize,mutationProfileMatrixSize);
		for(int i = 0; i < mutationProfileMatrixSize; i++){
			String s = sampleGenome();
			s = s.substring(1,s.length()-1);
			if(s.equals("")){
				continue;
			}
			List <String> tmp = MyFunc.split(", ", s);
			for(String t:tmp){
				int r = Integer.parseInt(t);
				M.set(r, i, 1.0);
			}
		}
		return M;
	}
	
	
	MyMat getMutaionSimilarityMatrix(){
		return getSimilarityMat(getMutaionProfileMatrix());
	}
	
	public MyMat getSimilarityMat(MyMat m){
		MyMat M = new MyMat(m.colSize(), m.colSize());
		for(int i = 0; i <  m.colSize(); i++){
			for(int j = 0; j <  i; j++){
				double tmp = 0;
				for(int k = 0; k <  m.rowSize(); k++){
					//tmp += Math.pow(m.get(k, i)-m.get(k, j), 2); 
					tmp += Math.abs(m.get(k, i)-m.get(k, j)); 
				}
				tmp = 1 - tmp/m.rowSize();
				M.set(i, j, tmp);
				M.set(j, i, tmp);
			}	
		}
		return M;
	}
	
	
	public  double getEigenValueEntropy(MyMat m) throws NotConvergedException{
		 DenseMatrix  M = m.getMTJDenseMatrix();
		 SVD svd = SVD.factorize(M);
		 double[] sv =  svd.getS();
		 double ic = shannonEntropy(sv);
		 return  ic;
	}
	
	public  double shannonEntropy(double[] v){
		   double e = 0;
		   double s = 0;
		   List <Double> v2 = new ArrayList<Double>();
		   for(int i = 0, n = v.length; i < n; i++){
			   s  +=  v[i];
			   v2.add(v[i]);
		   }
		   eigenValue =  v2;
		   for(Double d: v2){
			   if(d > 0){
			    double p = d/s;
			    e += p*Math.log(p);
			   }
		   }
		   return - e/Math.log(v.length);   
	 }
	
	
	
	public void getMutatedGeneCount(){
		mutatedDriverGeneCount = 0;
		mutatedGeneCount = 0;
		averageFitness = 0;
		Map<String, Double> populationCoverage = new HashMap <String ,Double>();
		for(String s: subpopulationProportion.keySet()){			
			double f = 1;
			String s2 = s.substring(1,s.length()-1);
			List<String> mutatedGenes = MyFunc.split(", ", s2);
			if(s2.equals("")){
				continue;
			}
			mutatedGeneCount += mutatedGenes.size()*subpopulationProportion.get(s);
			int driver =0;
			for(String g: mutatedGenes){
				int i = Integer.valueOf(g);
				if(i < driverSize){
					driver++;
					f *= fitnessIncrease;
				}
				if(populationCoverage.containsKey(g)){
					populationCoverage.put(g, populationCoverage.get(g) + subpopulationProportion.get(s));
				}else{
					populationCoverage.put(g, subpopulationProportion.get(s));
				}
			}
			mutatedDriverGeneCount += driver*subpopulationProportion.get(s);
			averageFitness  += f*subpopulationProportion.get(s);
		}
		double founder = 0;
		for(String g:  populationCoverage.keySet()){
			if(populationCoverage.get(g) >= populationCoverageForfounderMutation){
				founder++;
			}
		}
		founderMutationPropotion = founder/mutatedGeneCount;
		
	}
	
	String sampleGenome(){
		while(true){
			Double r = R.nextDouble();
			for(String s: subpopulationProportion.keySet()){
				Double w = subpopulationProportion.get(s);
				r -= w;
				if(r < 0){
					return s;
				}
			}
		}
	}
	
	
	public static void main(String [] args) throws Exception{	
		Options options = new Options();
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		SimulatorOld S = new SimulatorOld();
		options.addOption("m", "mut", true,  "mutaionRate (" + S.mutationRate  + ")" );
		options.addOption("g", "grow", true, "growthRate (" + S.growthRate  + ")");
		options.addOption("D", "death", true, "deathRate (" + S.deathRate  + ")");
		options.addOption("N", "deathns", true, "deathRate for non stem cell (" + S.deathRateForNonStem  + ")");
		options.addOption("P", "maxpop", true,  "maxPopulationSize (" + S.maxPopulationSize  + ")");
		options.addOption("G", "gen", true,  "genomeSize (" + S.genomeSize  + ")");
		options.addOption("d", "drv", true,  "driverSize (" + S.driverSize  + ")");
		options.addOption("f", "fit", true,  "fitnessIncrease (" + S.fitnessIncrease  + ")");
		options.addOption("p", "inipop", true,  "initialPopulationSize (" + S.initialPopulationSize  + ")");
		options.addOption("T", "maxtime", true,  "maxTime (" + S.maxTime  + ")");
		options.addOption("M", "mutp", true,  "get mutation profile");
		options.addOption("s", "stat", true,  "simulate while printing statistics");
		options.addOption("S", "sym", true,  "symmetric replication probablity (" +  S.symmetricReplicationProbablity + ")");
		options.addOption("e", "eigen", true,  "print eigen values");
		options.addOption("w", "gateway", true,  "the number of gateway drivers");
		
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] ", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 0){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] ", options);
			return;
		}
		

		if(commandLine.hasOption("m")){
			S.mutationRate = Double.valueOf(commandLine.getOptionValue("m"));
		}
		if(commandLine.hasOption("g")){
			S.growthRate = Double.valueOf(commandLine.getOptionValue("g"));
		}
		if(commandLine.hasOption("D")){
			S.deathRate = Double.valueOf(commandLine.getOptionValue("D"));
		}
		if(commandLine.hasOption("N")){
			S.deathRateForNonStem = Double.valueOf(commandLine.getOptionValue("N"));
		}
		if(commandLine.hasOption("p")){
			S.initialPopulationSize = Integer.valueOf(commandLine.getOptionValue("p"));
		}
		if(commandLine.hasOption("P")){
			S.maxPopulationSize = Integer.valueOf(commandLine.getOptionValue("P"));
		}
		if(commandLine.hasOption("G")){
			S.genomeSize = Integer.valueOf(commandLine.getOptionValue("G"));
		}
		if(commandLine.hasOption("d")){
			S.driverSize = Integer.valueOf(commandLine.getOptionValue("d"));
		}
		if(commandLine.hasOption("f")){
			S.fitnessIncrease = Double.valueOf(commandLine.getOptionValue("f"));
		}
		if(commandLine.hasOption("T")){
			S.maxTime = Integer.valueOf(commandLine.getOptionValue("T"));
		}
		if(commandLine.hasOption("S")){
			S.symmetricReplicationProbablity  = Double.valueOf(commandLine.getOptionValue("S"));
		}
		if(commandLine.hasOption("w")){
			S.gatewayDriver = Integer.valueOf(commandLine.getOptionValue("w"));
		}
		if(commandLine.hasOption("s")){
			S.simulateWhilePrintingStatistics(Integer.valueOf(commandLine.getOptionValue("s")));
		}else{
			S.simulate();
			S.getStatistics();
			System.out.println(S.header);
			System.out.println(S);
		}
		if(commandLine.hasOption("e")){
			PrintWriter os = new PrintWriter(new FileWriter(commandLine.getOptionValue("e")));
			os.println(MyFunc.join("\t", S.eigenValue));
			os.flush();
			os.close();
		}
		if(commandLine.hasOption("M")){
			PrintWriter os = new PrintWriter(new FileWriter(commandLine.getOptionValue("M")));
			os.println(S.getMutaionProfileMatrix());
			os.flush();
			os.close();
		}
	}
}
