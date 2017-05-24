package sim;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;


import sun.reflect.Reflection;
import utility.MyFunc;

public class SimulatorWithStemCellAndGateWayDriver extends SimulatorWithStemCell{
	public  SimulatorWithStemCellAndGateWayDriver (){}
	
	List <Cell> population;
	
	public  SimulatorWithStemCellAndGateWayDriver  (SimulatorWithStemCellAndGateWayDriver  S){
		mutaionRate = S.mutaionRate;
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
			boolean gateway = false; 
			if(genome.isMutated(0)){
				gateway = true;
			}
			for(int i=0; i< driverSize; i++){
				if(genome.isMutated(i)){
					fitness *= fitnessIncrease;
				}
			}
			if(!gateway  & fitness > 1){
				fitness = -1;
			}
		}
		
		String getGenomeString(){
			return genome.mutatedGenes.toString();
		}
		
	}
	
	void initialize(){
		population = new ArrayList<Cell>();
		for(int i = 0; i < initialPopulationSize; i++){
			population.add(new Cell());
		}
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
	
	
	void getStatistics(){
		populationSize = population.size();
		getSubpopulationProportion();
		getMajorSubpopulationProportion();
		getPopulationEntropy();
		getMutatedGeneCount();
		getMutatedDriverGeneCountAndAverageFitness();
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
	
		
	public static void main(String [] args) throws Exception{	
		Options options = new Options();
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		SimulatorWithStemCellAndGateWayDriver S = new SimulatorWithStemCellAndGateWayDriver();
		options.addOption("m", "mut", true,  "mutaionRate (" + S.mutaionRate  + ")" );
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
			S.mutaionRate = Double.valueOf(commandLine.getOptionValue("m"));
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
