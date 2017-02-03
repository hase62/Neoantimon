package eem;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.math.random.*;


import sun.reflect.Reflection;
import utility.*;

public class SimulationDataGenerator {
	int geneSize = 6000;
	int sampleSize = 100;
	int geneSetSize = 300;
	int negativeGeneSetNumber = 10;
	int positiveGeneSetNumber = 10;
	int moduleNumber = 20;
	int moduleSize = 300;
	double moduleGeneSubsetRate = 0.1;
	double signalStrength = 0.1;
	
	boolean absolute = false;	

	RandomGenerator RG;
	
	MyMat Exp;
	Map <String, List <String>> geneSet; 
	List <String> geneName;
	List <String> sampleName;
	
	
	public SimulationDataGenerator(){
		RG = new JDKRandomGenerator();
		geneSet = new HashMap<String, List<String>>();
		geneName = new ArrayList<String>();
	}
	
	public void setGeneSetSize(int i){
		geneSetSize = i;
	}
	public void setGeneSetNumber(int positiveGeneSetNumber, int negativeGeneSetNumber){
		this.positiveGeneSetNumber = positiveGeneSetNumber;
		this.negativeGeneSetNumber = negativeGeneSetNumber;
	}
	public void setModuleGeneSubsetRate(double d){
		moduleGeneSubsetRate = d;
	}
	public void setSignalStrength(double d){
		signalStrength= d;
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
		for(int k = 0; k < moduleNumber; k++){
			List <Double> v = Exp.getRow(k*moduleSize);
			for(int i = 1; i < moduleSize; i++){
				if(absolute == false){
					for(int j = 0; j < sampleSize; j++){
						Exp.set(i + k*moduleSize, j, (signalStrength)*v.get(j) + (1-signalStrength)*RG.nextGaussian());
					}
				}else{
					if(RG.nextGaussian()>0){
						for(int j = 0; j < sampleSize; j++){
							Exp.set(i + k*moduleSize, j, (signalStrength)*v.get(j) + (1-signalStrength)*RG.nextGaussian());
						}
					}else{
						for(int j = 0; j < sampleSize; j++){
							Exp.set(i + k*moduleSize, j, -(signalStrength)*v.get(j) + (1-signalStrength)*RG.nextGaussian());
						}
					}	
				}
			}
		}
	}
	protected void simulateGeneset(){
		int moduleGeneSubsetSize = (int) Math.round(geneSetSize*moduleGeneSubsetRate);
		Random rnd = new Random();
		int nonModuleGeneSubsetSize = geneSetSize - moduleGeneSubsetSize;
		for(int i = 0; i < positiveGeneSetNumber; i++){
			int j = rnd.nextInt(moduleNumber);
			List<String> tmp = new ArrayList<String>();
			tmp.addAll(geneName.subList(j*moduleSize, j*moduleSize + moduleGeneSubsetSize));
			tmp.addAll(MyFunc.sample(MyFunc.diff(geneName, tmp), nonModuleGeneSubsetSize));
			geneSet.put("positive" + (i+1),tmp);
		}
		for(int i = 0; i < negativeGeneSetNumber; i++){
			geneSet.put("negative"+(i+1), MyFunc.sample(geneName, geneSetSize));
		}
		
	}
	
	public void simulate(){
		simulateExpression();
		simulateGeneset();
	}
	
	
	public MyMat getExpression(){
		return Exp;
	}
	
	public Map <String, List<String>> getGeneSet(){
		return geneSet;
	}
	
	
	public void printGeneSets(String outfile) throws IOException {
		PrintWriter os = new PrintWriter(new FileWriter(outfile));
		for(String s : geneSet.keySet()){
			os.println(s + "\t\t" + MyFunc.join("\t",geneSet.get(s)));
		}
		os.close();
	}
	                                   
	public static void main(String [] args) throws Exception{Options options = new Options();
		options.addOption("g", "gene", true, "module gene subset ratio");
		options.addOption("S", "signal", true, "signal strength");
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
		SimulationDataGenerator S;
		if(commandLine.hasOption("c") || commandLine.hasOption("K")){
			S = new SimulationDataGeneratorUsingCompositeModule();	
		}else{
			S  = new SimulationDataGenerator();
		}
		if(commandLine.hasOption("K")){		
			((SimulationDataGeneratorUsingCompositeModule)S).setModuleNumberInOneGeneSet(Integer.valueOf(commandLine.getOptionValue("K")));
		}
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
