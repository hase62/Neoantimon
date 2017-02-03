package eem;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;


import sun.reflect.Reflection;
import utility.MyFunc;
import utility.MyMat;

public class SimulationDataGeneratorUsingSineWaveModule extends SimulationDataGenerator {
	
	protected void simulateExpression(){
		for(int i = 1; i <= geneSize; i++){
			geneName.add("gene" + i);
		}
		sampleName = new ArrayList<String>();
		for(int i = 1; i <= sampleSize; i++){
			sampleName.add("sample" + i);
		}
		Exp = new MyMat(geneName, sampleName);
		for(int i = 0; i < geneSize; i++){
			for(int j = 0; j < sampleSize; j++){
				Exp.set(i, j, RG.nextGaussian());
			}
		}
		double phase_shift = Math.PI/8;
		for(int k = 0; k < moduleNumber; k++){
			List <Integer> index = new ArrayList<Integer>();
			for(int i = 0; i < sampleSize; i++){
				index.add(i);
			}
			index = MyFunc.sample(index, index.size());
			for(int i = 0; i < moduleSize; i++){
				for(int j = 0; j < sampleSize; j++){
					Exp.set(i + k*moduleSize, index.get(j), (signalStrength)*Math.sin(phase_shift*(i+j)) + (1-signalStrength)*RG.nextGaussian());
				}
				
			}
		}
	}
	
	private int moduleNumberInOneGeneSet = 1;

	public void setModuleNumberInOneGeneSet(int k){
		moduleNumberInOneGeneSet = k;
	}
	
	public void simulateGeneset(){
		Random rnd = new Random();
		int moduleGeneSubsetSize = (int) Math.round(geneSetSize*moduleGeneSubsetRate);
		for(int i = 0; i < positiveGeneSetNumber; i++){
			List<String> tmp = new ArrayList<String>();
			Set <Integer> seen = new HashSet<Integer>();
			for(int j = 0; j < moduleNumberInOneGeneSet;j++){
				int k = rnd.nextInt(moduleNumber);
				if(seen.contains(k)){
					j--;
					continue;
				}else{
					tmp.addAll(geneName.subList(k*moduleSize, k*moduleSize + moduleGeneSubsetSize/moduleNumberInOneGeneSet));
					seen.add(k);
				}
			}
			tmp.addAll(MyFunc.sample(MyFunc.diff(geneName, tmp), geneSetSize - tmp.size()));
			geneSet.put("positive" + (i+1),tmp);
		}
		for(int i = 0; i < negativeGeneSetNumber; i++){
			geneSet.put("negative"+(i+1), MyFunc.sample(geneName, geneSetSize));
		}
		
	}
	
	public static void main(String [] args) throws Exception{Options options = new Options();
	options.addOption("g", "gene", true, "module gene subset ratio");
	options.addOption("S", "signal", true, "signal strength");
	options.addOption("R", "rowsize", true, "row (gene) size of exp");
	options.addOption("C", "colsize", true, "col (sample) size of exp");
	options.addOption("G", "gssize", true, "gene set size");
	options.addOption("m", "emsize", true, "exp mod size");
	options.addOption("k", "gsno", true, "gene set number");
	options.addOption("n", "emno", true, "exp mod number");
	options.addOption("K", "compnum", true, "composite module number");
	
	HelpFormatter formatter = new HelpFormatter();
	CommandLineParser parser = new BasicParser();
	CommandLine commandLine;
	try{
		commandLine = parser.parse(options, args);
	}catch (Exception e) {
		formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options]  expfile gmtfile", options);
		return ;
	}
	List <String> argList = commandLine.getArgList();
	if(argList.size() != 2){
		formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] expfile gmtfile", options);
		return;
	}
	SimulationDataGeneratorUsingSineWaveModule S;
	S  = new SimulationDataGeneratorUsingSineWaveModule();
	if(commandLine.hasOption("g")){
		S.setModuleGeneSubsetRate(Double.valueOf(commandLine.getOptionValue("g")));
	}
	if(commandLine.hasOption("S")){
		S.setSignalStrength(Double.valueOf(commandLine.getOptionValue("S")));
	}
	if(commandLine.hasOption("K")){		
		S.setModuleNumberInOneGeneSet(Integer.valueOf(commandLine.getOptionValue("K")));
	}
	if(commandLine.hasOption("R")){
		S.geneSize = Integer.valueOf(commandLine.getOptionValue("R"));
	}
	if(commandLine.hasOption("C")){
		S.sampleSize = Integer.valueOf(commandLine.getOptionValue("C"));
	}
	if(commandLine.hasOption("G")){
		S.geneSetSize = Integer.valueOf(commandLine.getOptionValue("G"));
	}
	if(commandLine.hasOption("m")){
		S.moduleSize = Integer.valueOf(commandLine.getOptionValue("m"));
	}
	if(commandLine.hasOption("k")){
		S.positiveGeneSetNumber = Integer.valueOf(commandLine.getOptionValue("k"));
		S.negativeGeneSetNumber = Integer.valueOf(commandLine.getOptionValue("k"));
	}
	if(commandLine.hasOption("n")){
		S.moduleNumber = Integer.valueOf(commandLine.getOptionValue("n"));
	}
	
	
	S.simulate();
	S.printGeneSets(argList.get(1));
	S.getExpression().print(argList.get(0));	
}

}
