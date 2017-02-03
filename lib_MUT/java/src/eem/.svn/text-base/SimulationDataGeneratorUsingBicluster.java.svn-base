package eem;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;


import sun.reflect.Reflection;
import utility.MyFunc;
import utility.MyMat;

public class SimulationDataGeneratorUsingBicluster extends SimulationDataGenerator{
	int moduleNumber = 50;
	double moduleSampleSubsetRate = 0.1;
	
	protected List <List<String>> moduleGeneList;
	
	public void setModuleSampleSubsetRate(double d){
		moduleSampleSubsetRate = d;
	}
	
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
		
		moduleGeneList = new ArrayList<List<String>>();
		for(int k = 0; k < moduleNumber; k++){
			List <String> moduleGene = new ArrayList<String>(MyFunc.sample(geneName, moduleSize));
			List <String> moduleSample = new ArrayList<String>(MyFunc.sample(sampleName, (int)Math.round(sampleSize*moduleSampleSubsetRate)));
 			moduleGeneList.add(moduleGene);
			for(String s: moduleGene){
				for(String t: moduleSample){
					Exp.set(s, t, Exp.get(s,t) + signalStrength);
				}
			}
		}
	}
	protected void simulateGeneset(){
		int moduleGeneSubsetSize = (int) Math.round(geneSetSize*moduleGeneSubsetRate);
		int nonModuleGeneSubsetSize = geneSetSize - moduleGeneSubsetSize;
		Random rnd = new Random();
		for(int i = 0; i < positiveGeneSetNumber; i++){
			List<String> tmp = new ArrayList<String>();
			int j = rnd.nextInt(moduleNumber);
			tmp.addAll(new ArrayList<String>(moduleGeneList.get(j).subList(0,moduleGeneSubsetSize)));
			tmp.addAll(MyFunc.sample(MyFunc.diff(geneName, tmp), nonModuleGeneSubsetSize));
			geneSet.put("positive" + (i+1),tmp);
		}
		for(int i = 0; i < negativeGeneSetNumber; i++){
			geneSet.put("negative"+(i+1), MyFunc.sample(geneName, geneSetSize));
		}
	}
	
	public static void main(String [] args) throws Exception{
	Options options = new Options();
	options.addOption("g", "gene", true, "module gene subset ratio");
	options.addOption("S", "signal", true, "signal strength");
	options.addOption("s", "sample", true, "module sample subset ratio");
	options.addOption("c", "composite", false, "use composite modules");
	options.addOption("K", "compnum", true, "composite module number");
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
	SimulationDataGeneratorUsingBicluster S;
	if(commandLine.hasOption("c") || commandLine.hasOption("K")){
		S = new SimulationDataGeneratorUsingBiclusterCompositeModule();
	}else{
		S= new SimulationDataGeneratorUsingBicluster();
	}
	if(commandLine.hasOption("K")){		
		((SimulationDataGeneratorUsingBiclusterCompositeModule)S).setModuleNumberInOneGeneSet(Integer.valueOf(commandLine.getOptionValue("K")));
	}
	if(commandLine.hasOption("g")){
		S.setModuleGeneSubsetRate(Double.valueOf(commandLine.getOptionValue("g")));
	}
	if(commandLine.hasOption("S")){
		S.setSignalStrength(Double.valueOf(commandLine.getOptionValue("S")));
	}
	if(commandLine.hasOption("s")){
		S.setModuleSampleSubsetRate(Double.valueOf(commandLine.getOptionValue("s")));
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
