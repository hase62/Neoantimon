package sim;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

import sun.reflect.Reflection;
import utility.MyFunc;

public class MultiSimulator {
	
	List<Double> entropy;
	List<Double> time;
	List<Double> mutatedGenePropotion;
	SimulatorWithStemCell S;
	int n = 5;
	
	public MultiSimulator(SimulatorWithStemCell S){
		entropy = new ArrayList<Double>();
		time = new ArrayList<Double>();
		 mutatedGenePropotion = new ArrayList<Double>();
		this.S = S;
	}
	
	public void simulate(){
		for(int i =0; i < n; i++){
			SimulatorWithStemCell s = new SimulatorWithStemCell(S);
			System.err.println("the " + (i+1) + " th trial......");
			s.simulate();
			s.getStatistics();
			System.err.println(s.header);
			System.err.println(s);
			if(s.time < S.maxTime){
				entropy.add(s.populationEntropy);
				time.add((double)s.time);
				mutatedGenePropotion.add(s.mutatedGeneCount);
			}else{
				System.err.println("warn:time exceeded");
			}
			System.err.println();
		}
	}
	
	public String toString(){
		if(!entropy.isEmpty()){
			String tmp = "timeMean\ttimeSD\tentropyMean\tentropySD\tmutatedGeneProportionMean\tmutatedGeneProportionSD\n";
			tmp += MyFunc.mean(time) + "\t" + MyFunc.sd(time)  + "\t" + MyFunc.mean(entropy) + "\t" + MyFunc.sd(entropy)  + "\t" + MyFunc.mean(mutatedGenePropotion) + "\t" + MyFunc.sd(mutatedGenePropotion)+"\n";
			return tmp;
		}else{
			System.err.println("warn:time exceeded for all trials");
			return "";
		}
	}
	
	
	public static void main(String [] args) throws Exception{	
		Options options = new Options();
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		 SimulatorWithStemCell S = new SimulatorWithStemCell();
		options.addOption("m", "mut", true,  "mutaionRate (" + S.mutaionRate  + ")" );
		options.addOption("g", "grow", true, "growthRate (" + S.growthRate  + ")");
		options.addOption("D", "death", true, "deathRate (" + S.deathRate  + ")");
		options.addOption("N", "deathns", true, "deathRate for non stem cell (" + S.deathRateForNonStem  + ")");
		options.addOption("P", "maxpop", true,  "maxPopulationSize (" + S.maxPopulationSize  + ")");
		options.addOption("G", "gen", true,  "genomeSize (" + S.genomeSize  + ")");
		options.addOption("d", "drv", true,  "driverProportion (" + S.driverSize  + ")");
		options.addOption("f", "fit", true,  "fitnessIncrease (" + S.fitnessIncrease  + ")");
		options.addOption("p", "inipop", true,  "initialPopulationSize (" + S.initialPopulationSize  + ")");
		options.addOption("T", "maxtime", true,  "maxTime (" + S.maxTime  + ")");
		options.addOption("S", "sym", true,  "symmetric replication probablity (" +  S.symmetricReplicationProbablity + ")");
		options.addOption("n", "ntrial", true,  "numberOfTrial (" + 10 + ")");
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
		
		MultiSimulator MS = new MultiSimulator(S);
		if(commandLine.hasOption("n")){
			MS.n = Integer.valueOf(commandLine.getOptionValue("n"));
		}
		MS.simulate();
		System.out.print(MS);
		
	}

}
