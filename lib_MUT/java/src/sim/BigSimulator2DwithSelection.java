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

public class BigSimulator2DwithSelection extends BigSimulator2D{
	//parameter
	int selectiveDriverSize = 4;	
	
	void initializeGenomes(){
		populationSize = 0;
		time = 0;
		limXY = (int) (5*Math.pow(maxPopulationSize, 0.5));
		G = genomeSize/maxBinary + 1; 
		genomes  = new long[2*limXY+1][2*limXY+1][G];
		state= new boolean [2*limXY+1][2*limXY+1][3];
		genome2fitnessMap =  new HashMap <Long, Double>();
		genome2mutatedGeneMap = new HashMap <Long, Set<Integer>>();
		fillCore();
	}
	
	double genome2fitness(long[] genome){
		if(genome2fitnessMap.containsKey(genome[0])){
			return genome2fitnessMap.get(genome[0]);
		}else{
			double fitness = 1; 
			for(int i=0; i< driverSize - selectiveDriverSize; i++){
				if(isMutated(genome,i)){
					fitness *= fitnessIncrease;
				}
			}
			genome2fitnessMap.put(genome[0], fitness);
			return fitness;
		}
	}
	
	double getFitness(int x, int y){
		double fitness = genome2fitness(get(x,y));
		for(int i=0; i< selectiveDriverSize; i++){
			int g = i+ driverSize - selectiveDriverSize;
			/*if(!isMutated(x, y ,g) & hasSelection(x, y, g)){
				fitness = -1;
			}*/
			if(isMutated(x, y, g) & hasSelection(x, y, g)){
				fitness *= fitnessIncrease;
			}
		}
		
		return fitness;
	}
	
	
	boolean hasSelection(int x, int y,  int g){
		int i = g + selectiveDriverSize - driverSize;
		double midX = (minX + maxX)/2;
		double midY = (minY + maxY)/2;
		if(i==0){
			if(x > midX & y > midY){
				return  true;
			}else{
				return false;
			}
		}else if(i==1){
			if(x > midX & y < midY){
				return  true;
			}else{
				return false;
			}
		}else if(i==2){
			if(x < midX & y < midY){
				return  true;
			}else{
				return false;
			}
		}else if(i==3){
			if(x < midX & y > midY){
				return  true;
			}else{
				return false;
			}
		}
		return (Boolean) null;
	}
	
	
	public String getParmetersAsString(){
		StringBuffer S = new StringBuffer(super.getParmetersAsString());
		S.append("selectiveDriverSize" + "\t" + selectiveDriverSize + "\n");
		return S.toString();
	}
	
	
	
	public static void main(String [] args) throws Exception{	
		Options options = new Options();
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		BigSimulator2DwithSelection S = new  BigSimulator2DwithSelection();
		String outfile = "out";
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
		
		options.addOption("o", "outfile", true,  "out filename");
		options.addOption("s", "snap", true,  "take snap shot");
		
		
		
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
		if(commandLine.hasOption("o")){
			outfile = commandLine.getOptionValue("o");
		}		
		S.initializeGenomes();
		if(!commandLine.hasOption("s")){
			S.simulate();
		}else{
			S.simulateWhilePrinting(outfile, Integer.valueOf(commandLine.getOptionValue("s")));
		}
		S.printResults(outfile);
	}
	
	
}
