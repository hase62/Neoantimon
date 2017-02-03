package sim;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

import sun.reflect.Reflection;

public class BigSimulator2Ddrug extends BigSimulator2D{
	double deathRateDrug =  0.9;
	boolean drugTreatmant = false;
	
	boolean isResistant(int x, int y){		
		return isMutated(get(x,y), driverSize);
	}

	boolean growCell1(int x, int y){
		if(isEmpty(x,y)){
			return true;
		}
		if(Math.abs(x) ==  limXY | Math.abs(y) ==  limXY){
			return false;
		}
		if(!drugTreatmant | drugTreatmant & isResistant(x,y)){
			if((isStem(x,y) & R.nextDouble() < deathRate) | 
				(!isStem(x,y) & R.nextDouble() < deathRateForNonStem)){
				clear(x,y);
				populationSize--;
				return true;
			}
		}else{
			if(R.nextDouble() < deathRateDrug){
				clear(x,y);
				populationSize--;
				return true;
			}
			
		}
		if(!drugTreatmant | drugTreatmant & isResistant(x,y)){
		if(R.nextDouble() < getFitness(x,y)*growthRate){
			mutate(x,y);
			setRep(x,y,true);
		}else{
			setRep(x,y,false);
		}
		}else{
			setRep(x,y,false);
		}
		if(populationSize >= maxPopulationSize){
			return false;
		}else{
			return true;
		}
	}
	public void simulateWhilePrinting(String f, int I) throws IOException{
		int i = 0;
		printResults(f + "." + i);
		int poplationSizeCutoff = populationSize*I;
		 for(time = 1; time <= maxTime; time++){
			 growPopulation();
			 //if(time/10 == (double)time/10){
				 System.err.println("time=" + time +"\t" + "populationSize=" + populationSize) ;
			//}
			 if(Math.abs(maxX) ==  limXY | Math.abs(minX) ==  limXY | Math.abs(maxY) ==  limXY | Math.abs(minY) ==  limXY){
				 System.err.println("warn: reached space limit!");
				 break;
			}
			
			if(poplationSizeCutoff < 0){
				i++;
				printResults(f + "." + i);
				poplationSizeCutoff = populationSize*I;
				System.err.println("print result!");
			}
			if(!drugTreatmant & populationSize>= maxPopulationSize){
				System.err.println("drug treatment  starts!");
				drugTreatmant = true;
				i++;
				printResults(f + "." + i);
				poplationSizeCutoff = -1;
				System.err.println("print result!");
				continue;
			}else if(populationSize >= poplationSizeCutoff){
				i++;
				printResults(f + "." + i);
				poplationSizeCutoff *= I;
				System.err.println("print result!");
			}else if( populationSize>= maxPopulationSize){
				i++;
				printResults(f + "." + i);
				poplationSizeCutoff *= I;
				System.err.println("print result!");
			}
			
			 if(populationSize == 0 | populationSize>= maxPopulationSize){
				 break;
			}
		}
		 System.err.println(minX + " " + maxX + " " + minY + " " + maxY + " " + limXY + " " + count());
			
	}
	public void simulate(){
		 for(time = 1; time <= maxTime; time++){
			 growPopulation();
			 if(time/10 == (double)time/10){
				 System.err.println("time=" + time +"\t" + "populationSize=" + populationSize) ;
			}
			 if(Math.abs(maxX) ==  limXY | Math.abs(minX) ==  limXY | Math.abs(maxY) ==  limXY | Math.abs(minY) ==  limXY){
				 System.err.println("warn: reached space limit!");
				 break;
			}
			 if(!drugTreatmant & populationSize>= maxPopulationSize){
				 System.err.println("drug treatment  starts!");
					drugTreatmant = true;
				}
			 if(populationSize == 0 | populationSize>= maxPopulationSize){
				 break;
			}
		}
		System.err.println(minX + " " + maxX + " " + minY + " " + maxY + " " + limXY + " " + count());
	}
	public static void main(String [] args) throws Exception{	
		Options options = new Options();
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		BigSimulator2Ddrug  S = new BigSimulator2Ddrug();
		String outfile = "out";
		options.addOption("m", "mut", true,  "mutaionRate (" + S.mutationRate  + ")" );
		options.addOption("D", "death", true, "deathRate (" + S.deathRate  + ")");
		options.addOption("N", "deathns", true, "deathRate for non stem cell (" + S.deathRateForNonStem  + ")");
		options.addOption("P", "maxpop", true,  "maxPopulationSize (" + S.maxPopulationSize  + ")");
		options.addOption("G", "gen", true,  "genomeSize (" + S.genomeSize  + ")");
		options.addOption("g", "grow", true, "growthRate (" + S.growthRate  + ")");
		options.addOption("d", "drv", true,  "driverSize (" + S.driverSize  + ")");
		options.addOption("f", "fit", true,  "fitnessIncrease (" + S.fitnessIncrease  + ")");
		options.addOption("F", "fitlog", true,  "fitnessIncrease in log10 scale (" + Math.log10(S.fitnessIncrease)  + ")");
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
		if(commandLine.hasOption("F")){
			S.fitnessIncrease = Math.pow(10, Double.valueOf(commandLine.getOptionValue("F")));
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
