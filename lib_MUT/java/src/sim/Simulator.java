package sim;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.*;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

import sun.reflect.Reflection;

public class Simulator {
	// parameters for simulations
	double mutationRate = 0.0001;
	double growthRate = 0.0001;
	double deathRate =  0.0000001;
	double deathRateForNonStem =  0.01;
	double symmetricReplicationProbablity = 0.1;
	int initialPopulationSize = 10;
	int maxTime = 10000000;
	int maxPopulationSize = 100000;
	int genomeSize = 30;
	int  driverSize = 10;
	double fitnessIncrease  = 5;
	
	//variable 
	List<Long> genomes;
	Map <Long, Double> genome2fitnessMap;
	int populationSize;
	int time;
	Map <Long, Set<Integer>> genome2mutatedGeneMap;	
		
	//for simulation
	Random R = new Random();

	static String decimal2binary(long d){
		return Long.toBinaryString(d);
	}
	

	static int binary2decimal(String numberString){
        int convertedNumber = 0; 
        for (int i = 0; i < numberString.length(); i++){
        	if (Integer.parseInt(Character.toString(numberString.charAt(i))) == 1){
        		convertedNumber = convertedNumber + (Integer.parseInt(Character.toString(numberString.charAt(i))) * (int)Math.pow(2, (numberString.length() - i - 1)));
        	}
        }
        return convertedNumber;	
	}
	
	static long binary2decimalLong(String numberString){
        long convertedNumber = 0L; 
        for (int i = 0; i < numberString.length(); i++){
        	if (Integer.parseInt(Character.toString(numberString.charAt(i))) == 1){
        		convertedNumber = convertedNumber + (Long.parseLong(Character.toString(numberString.charAt(i))) * (long)Math.pow(2, (numberString.length() - i - 1)));
        	}
        }
        return convertedNumber;	
	}
	
	List <Long> getGenomesAsList(){
		return genomes;
	}
	
	long mutatedGenes2genome(Set <Integer> mutatedGenes){
		String tmp = "";	
		for(int i = 0; i < genomeSize; i++){
			if(mutatedGenes.contains(i)){
				tmp = "1" + tmp;
			}else{
				tmp = "0" + tmp;
			}
		}
		long tmp2 =  binary2decimalLong(tmp)+1L;
		genome2mutatedGeneMap.put(tmp2, mutatedGenes);
		return tmp2;
	}
	
	Set <Integer>  genome2mutatedGenes(long genome){
		long g = Math.abs(genome);
		if(genome2mutatedGeneMap.containsKey(g)){
			return genome2mutatedGeneMap.get(g);
		}else{
			String tmp = decimal2binary(g-1L);
			Set <Integer> mutatedGenes = new TreeSet<Integer>();
			int l = tmp.length();
			for(int i = 0; i < l; i++){
				if(tmp.substring(l-i-1,l-i).equals("1")){
					mutatedGenes.add(i);
				}
			}
			genome2mutatedGeneMap.put(g, mutatedGenes);
			return mutatedGenes;	
		}
	}
	
	boolean isMutated(long genome, int g){
		return genome2mutatedGenes(genome).contains(g);
	}
	
	boolean isMutated(int i, int g){
		return isMutated(get(i), g);
	}
	
	double genome2fitness(long genome){
		if(genome2fitnessMap.containsKey(genome)){
			return genome2fitnessMap.get(genome);
		}else{
			double fitness = 1;
			for(int i=0; i< driverSize; i++){
				if(isMutated(genome,i)){
					fitness *= fitnessIncrease;
				}
			}
			genome2fitnessMap.put(genome, fitness);
			return fitness;
		}
	}
	
	void set(int i, long genome){
		genomes.set(i, genome);
	}
	
	long get(int i){
		return genomes.get(i);
	}
	
	public void clear(int i){
		genomes.remove(i);
	}
	

	double getFitness(int i){
		return genome2fitness(get(i));
	}	
	
	boolean isStem(int i){
		return (get(i) > 0);
	}
	
	
	void  mutate(int j){
		long genome = get(j);	 
		long sign  = (genome<0.1)?-1L:1L;
		Set <Integer> newMutatedGenes = null;	
		for(int i=0; i < genomeSize; i++){
			if(!isMutated(genome,i) & (R.nextDouble() < mutationRate)){
				if(newMutatedGenes==null){
					newMutatedGenes = new TreeSet <Integer>(genome2mutatedGenes(genome));
				}
				newMutatedGenes.add(i);
			}
		}
		if(newMutatedGenes==null){
			return;
		}
		set(j,  sign*mutatedGenes2genome(newMutatedGenes));
	}
	void growPopulation(){
		int n = genomes.size();
		for(int i = 0; i < n & i < genomes.size(); i++){
			if(R.nextDouble() < getFitness(i)*growthRate){
				mutate(i);
				replicate(i);
				populationSize++;
			}
			if((isStem(i) & R.nextDouble() < deathRate) | (!isStem(i) & R.nextDouble() < deathRateForNonStem)){
				clear(i);
				i--;
				populationSize--;
			}
			if(populationSize>= maxPopulationSize){
				 break;
			}
		}
	}
	
	void replicate(int i){
		if(symmetricReplicationProbablity>=1){
			genomes.add(get(i));
		}else{
			if(isStem(i)){
				if(R.nextDouble() < symmetricReplicationProbablity){
					genomes.add(get(i));
				}else{
					genomes.add(-get(i));
				}
			}else{
				genomes.add(get(i));
			}
		}
	}
	
	void initializeGenomes(){
		populationSize = initialPopulationSize;
		time = 0;
		genomes  = new ArrayList<Long>();
		for(int i = 0;i < initialPopulationSize; i++){
			genomes.add(1L);
		}
		genome2fitnessMap =  new HashMap <Long, Double>();
		genome2mutatedGeneMap = new HashMap <Long, Set<Integer>>();
	}
	
	
	public void simulate(){
		 for(time = 1; time <= maxTime; time++){
			 growPopulation();
			 if(time/10 == (double)time/10){
				 System.err.println("time=" + time +"\t" + "populationSize=" + populationSize) ;
			}
			if(populationSize>= maxPopulationSize){
				 break;
			}
		}
		 System.err.println(populationSize + "\t" + genomes.size());
	}
		
	
	public String toString(){
		StringBuffer S = new StringBuffer("");
		for(int i = 0; i< genomes.size(); i++){
			S.append(get(i) + "\n");	
		}
		return S.toString();
	}
	
	public String getParmetersAsString(){
		StringBuffer S = new StringBuffer("");
		S.append("mutationRate" + "\t" +  mutationRate + "\n");
		S.append("growthRate" + "\t" + growthRate + "\n");
		S.append("deathRate" + "\t" + deathRate + "\n");
		S.append("deathRateForNonStem" + "\t" + deathRateForNonStem + "\n");
		S.append("symmetricReplicationProbablity" + "\t" + symmetricReplicationProbablity + "\n");
		S.append("initialPopulationSize" + "\t" + initialPopulationSize + "\n");
		S.append("maxTime" + "\t" + maxTime + "\n");
		S.append("maxPopulationSize" + "\t" + maxPopulationSize + "\n");
		S.append("genomeSize" + "\t" + genomeSize + "\n");
		S.append("driverSize" + "\t" + driverSize + "\n");
		S.append("fitnessIncrease" + "\t" + fitnessIncrease + "\n");
		S.append("populationSize" + "\t" + populationSize + "\n");
		S.append("time" + "\t" + time + "\n");
		 return S.toString();
	}
	
	public static void main(String [] args) throws Exception{	
		Options options = new Options();
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		Simulator S = new Simulator();
		options.addOption("m", "mut", true,  "mutaionRate (" + S.mutationRate  + ")" );
		options.addOption("D", "death", true, "deathRate (" + S.deathRate  + ")");
		options.addOption("N", "deathns", true, "deathRate for non stem cell (" + S.deathRateForNonStem  + ")");
		options.addOption("P", "maxpop", true,  "maxPopulationSize (" + S.maxPopulationSize  + ")");
		options.addOption("G", "gen", true,  "genomeSize (" + S.genomeSize  + ")");
		options.addOption("g", "grow", true, "growthRate (" + S.growthRate  + ")");
		options.addOption("d", "drv", true,  "driverSize (" + S.driverSize  + ")");
		options.addOption("f", "fit", true,  "fitnessIncrease (" + S.fitnessIncrease  + ")");
		options.addOption("p", "inipop", true,  "initialPopulationSize (" + S.initialPopulationSize   + ")");
		options.addOption("T", "maxtime", true,  "maxTime (" + S.maxTime  + ")");
		options.addOption("S", "sym", true,  "symmetric replication probablity (" +  S.symmetricReplicationProbablity + ")");
		
		options.addOption("o", "pprm", true,  "print parameters");
		
		options.addOption("s", "stat", true,  "print statistics");
		
		
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
			S.initialPopulationSize  = Integer.valueOf(commandLine.getOptionValue("p"));
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
		S.initializeGenomes();
		S.simulate();
		System.out.println(S);
		
		if(commandLine.hasOption("o")){
			PrintWriter os = new PrintWriter(new FileWriter(commandLine.getOptionValue("o")));
			os.print(S.getParmetersAsString());
			os.flush();
			os.close();
		}
	}
}