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

public class Simulator2DwithSelection extends Simulator2D{

	//parameter
	int selectiveDriverSize = 4;
	double adaptation = 0;
	
	
	//variable
	int[][] selection;
	Map<Integer, Set<Integer>> selection2selectiveDriverGenes;
	
	
	int selectiveDriverGenes2selection(Set <Integer> selectiveDriverGenes){
		String tmp = "";	
		for(int i = 0; i < selectiveDriverSize; i++){
			if(selectiveDriverGenes.contains(i)){
				tmp += 1;
			}else{
				tmp += 0;
			}
		}
		int tmp2 =  binary2decimal(tmp);
		selection2selectiveDriverGenes.put(tmp2, selectiveDriverGenes);
		return tmp2;
	}
	
	Set <Integer>  selection2selectiveDriverGenes(int selection){
		if(selection2selectiveDriverGenes.containsKey(selection)){
			return selection2selectiveDriverGenes.get(selection);
		}else{
			String tmp = decimal2binary(selection);
			Set <Integer> selectiveDriverGenes = new TreeSet<Integer>();
			for(int i = 0; i < tmp.length(); i++){
				if(tmp.substring(i,i+1).equals("1")){
					selectiveDriverGenes.add(i);
				}
			}
			selection2selectiveDriverGenes.put(selection, selectiveDriverGenes);
			return selectiveDriverGenes;	
		}
	}
	
	void setSelection(int x, int y, int g){	
		int s = selection[coordinate2index(x)][coordinate2index(y)];
		Set <Integer >tmp = new TreeSet <Integer>(selection2selectiveDriverGenes(s));
		tmp.add(g);
		int s2  = selectiveDriverGenes2selection(tmp);
		selection[coordinate2index(x)][coordinate2index(y)] = s2;
	}
	
	int getSelection(int x, int y){
		return selection[coordinate2index(x)][coordinate2index(y)];
	}
	
	boolean hasSelection(int x, int y, int g){
		return selection2selectiveDriverGenes(selection[coordinate2index(x)][coordinate2index(y)]).contains(g);
	}
	
	void setSelectiveRegion(){
		selectiveDriverSize = 4;
		selection = new int[2*limXY+1][2*limXY+1];
		selection2selectiveDriverGenes = new HashMap<Integer, Set<Integer>>();
		
		for(int x = -limXY; x <= limXY; x++){
			for(int y = -limXY; y <= limXY; y++){
				if(x > 0 & y > 0){
					setSelection(x,y, 0);
				}else if(x > 0 & y < 0){
					setSelection(x,y, 1);
				}else if(x < 0 & y < 0){
					setSelection(x, y, 2);
				}else if(x < 0 & y > 0){
					setSelection(x,y, 3);
				}	
			}
			
		}
		
	}
	
	void calculateAdaptaion(){
		double N = 0;
		double n = 0;
		for(int x = -limXY; x <= limXY; x++){
			for(int y = -limXY; y <= limXY; y++){
				if(isEmpty(x,y)){
					continue;
				}
				for(int i = 0; i <  selectiveDriverSize; i++){
					int g = i + driverSize - selectiveDriverSize;
					N++;
					if(isMutated(x, y, g) & hasSelection(x,y,g)){
						n++;
					}else if(!isMutated(x, y, g) & !hasSelection(x,y,g)){
						n++;
					}
				}
			}
		}
		adaptation =  n/N;
	}
		
	
	double genome2fitness(long genome){
		if(genome2fitnessMap.containsKey(genome)){
			return genome2fitnessMap.get(genome);
		}else{
			double fitness = 1; 
			for(int i=0; i< driverSize - selectiveDriverSize; i++){
				if(isMutated(genome,i)){
					fitness *= fitnessIncrease;
				}
			}
			genome2fitnessMap.put(genome, fitness);
			return fitness;
		}
	}
	
	double getFitness(int x, int y){
		double fitness = genome2fitness(get(x,y));
		long genome = get(x,y);
		for(int i=0; i< selectiveDriverSize; i++){
			//if(!isMutated(genome,i+driverSize - selectiveDriverSize) & hasSelection(x, y, i)){
			//	fitness *= 0;
			//}
			
			if(isMutated(genome,i+driverSize - selectiveDriverSize) & hasSelection(x, y, i)){
				fitness *= fitnessIncrease;
			}
		}
		
		return fitness;
	}
	
	public String getParmetersAsString(){
		StringBuffer S = new StringBuffer(super.getParmetersAsString());
		S.append("selectiveDriverSize" + "\t" + selectiveDriverSize + "\n");
		calculateAdaptaion();
		S.append("adaptation" + "\t" + adaptation + "\n");
		return S.toString();
	}
	
	public static void main(String [] args) throws Exception{	
		Options options = new Options();
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		Simulator2DwithSelection S = new Simulator2DwithSelection();
		options.addOption("m", "mut", true,  "mutaionRate (" + S.mutationRate  + ")" );
		options.addOption("D", "death", true, "deathRate (" + S.deathRate  + ")");
		options.addOption("N", "deathns", true, "deathRate for non stem cell (" + S.deathRateForNonStem  + ")");
		options.addOption("P", "maxpop", true,  "maxPopulationSize (" + S.maxPopulationSize  + ")");
		options.addOption("G", "gen", true,  "genomeSize (" + S.genomeSize  + ")");
		options.addOption("d", "drv", true,  "driverSize (" + S.driverSize  + ")");
		options.addOption("f", "fit", true,  "fitnessIncrease (" + S.fitnessIncrease  + ")");
		options.addOption("p", "inipop", true,  "initialPopulationCoreSize (" + S.initialPopulationCoreSize   + ")");
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
		if(commandLine.hasOption("D")){
			S.deathRate = Double.valueOf(commandLine.getOptionValue("D"));
		}
		if(commandLine.hasOption("N")){
			S.deathRateForNonStem = Double.valueOf(commandLine.getOptionValue("N"));
		}
		if(commandLine.hasOption("p")){
			S.initialPopulationCoreSize  = Integer.valueOf(commandLine.getOptionValue("p"));
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
		S.equilibratePopulation();
		S.setSelectiveRegion();
		S.simulate();
		System.out.println(S);
		
		if(commandLine.hasOption("o")){
			PrintWriter os = new PrintWriter(new FileWriter(commandLine.getOptionValue("o")));
			os.print(S.getParmetersAsString());
			os.flush();
			os.close();
		}
		if(commandLine.hasOption("s")){
			PrintWriter os = new PrintWriter(new FileWriter(commandLine.getOptionValue("s")));
			StatisticsCalculator2D SC = new StatisticsCalculator2D(S);
			SC.calculateStatistics();
			os.print(SC);
			os.flush();
			os.close();
		}
	}
	
	
}
