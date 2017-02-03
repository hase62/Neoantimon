package eem;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;


import sun.reflect.Reflection;
import utility.MyMat;

public class SimulationDataGeneratorUsingSingleSampleModule   extends SimulationDataGenerator {
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
		Random rnd = new Random();
		for(int k = 0; k < moduleNumber; k++){
			int j = rnd.nextInt(sampleSize);
			for(int i = 1; i < moduleSize; i++){
				if(absolute == false){
					Exp.set(i + k*moduleSize, j, Exp.get(i + k*moduleSize,j) + signalStrength);
				}else{
					Exp.set(i + k*moduleSize, j, Exp.get(i + k*moduleSize,j) + (RG.nextGaussian()>0?1:-1)*signalStrength);
				}
			}
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
	SimulationDataGeneratorUsingSingleSampleModule S = new SimulationDataGeneratorUsingSingleSampleModule();
	if(commandLine.hasOption("g")){
		S.setModuleGeneSubsetRate(Double.valueOf(commandLine.getOptionValue("g")));
	}
	if(commandLine.hasOption("S")){
		S.setSignalStrength(Double.valueOf(commandLine.getOptionValue("S")));
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
