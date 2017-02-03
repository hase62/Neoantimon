package eem;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.*;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

import eem.BiclusterBasedEEMsearch.BiclusterType;
import sun.reflect.Reflection;
import utility.*;

public class ParameterOptimizingBiclusterBasedEEMsearch  extends AbstractParameterOptimizingEEMsearch{
	private List <Double> Tgs;
	private List <Double> Tcs;
	private List <BiclusterType> biclusterTypes;
	public ParameterOptimizingBiclusterBasedEEMsearch(MyMat E){
		super();
		EEMsearchTemplate  = new BiclusterBasedEEMsearch(E);
		Tgs = new ArrayList<Double>();
		Tgs.add(0.05);
		Tgs.add(0.1);
		Tgs.add(0.15);
		Tcs = new ArrayList<Double>();
		Tcs.add(0.1);
		Tcs.add(0.2);
		Tcs.add(0.3);
		biclusterTypes = new ArrayList<BiclusterType>();
		biclusterTypes.add(BiclusterType.UP);
	}
	
	public void setTgs(List <Double> Tgs){
		this.Tgs = new ArrayList<Double>(Tgs);
	}
	public void setTcs(List <Double> Tcs){
		this.Tcs = new ArrayList<Double>(Tcs);
	}
	
	
	public void setMaxItrForISA(int i){
		((BiclusterBasedEEMsearch)EEMsearchTemplate).setMaxItrForISA(i);
	}
	public void setExtrapolationCutoff(double d){
		((BiclusterBasedEEMsearch)EEMsearchTemplate).setExtrapolationCutoff( d);
	}


	public void setBiclusterType(BiclusterType type){
		biclusterTypes.clear();
		biclusterTypes.add(type);
	}
	public void useTwoBiclusterTypes(){
		biclusterTypes.clear();
		biclusterTypes.add(BiclusterType.UP);
		biclusterTypes.add(BiclusterType.DOWN);
	}
	public void useAllBiclusterTypes(){
		biclusterTypes.clear();
		biclusterTypes.add(BiclusterType.UP);
		biclusterTypes.add(BiclusterType.DOWN);
		biclusterTypes.add(BiclusterType.ABSOLUTE);
		biclusterTypes.add(BiclusterType.gUPcDOWN);
	}
	
	protected void initializeEEMsearches(){
		((BiclusterBasedEEMsearch)EEMsearchTemplate).setBiclusterType(biclusterTypes.get(0));
		for(int i = 0; i < Tgs.size(); i++){
			for(int j = 0; j < Tcs.size(); j++){
				BiclusterBasedEEMsearch eem  = new BiclusterBasedEEMsearch((BiclusterBasedEEMsearch) EEMsearchTemplate, false);
				eem.setTg(Tgs.get(i));
				eem.setTc(Tcs.get(j));
				EEMsearches.add(eem);
			}
		}
		for(int k = 1; k < biclusterTypes.size(); k++){
			BiclusterBasedEEMsearch EEMsearchTemplateDeep = new  BiclusterBasedEEMsearch((BiclusterBasedEEMsearch) EEMsearchTemplate,true);
			EEMsearchTemplateDeep.setBiclusterType(biclusterTypes.get(k));
			for(int i = 0; i < Tgs.size(); i++){
				for(int j = 0; j < Tcs.size(); j++){
					BiclusterBasedEEMsearch eem  = new BiclusterBasedEEMsearch(EEMsearchTemplateDeep, false);
					eem.setTg(Tgs.get(i));
					eem.setTc(Tcs.get(j));
					EEMsearches.add(eem);
				}
			}			
		}
	}
	public String toString(){
		List <String> tmp = new ArrayList<String>();
		tmp.add("biclusterType= " + biclusterTypes);
		tmp.add("Tgs= " + Tgs);
		tmp.add("Tcs= " + Tcs);
		tmp.add("maxGeneSetSize= " + ((BiclusterBasedEEMsearch)EEMsearchTemplate).maxSeedGeneSize);
		tmp.add("minGeneSetSize= " + ((BiclusterBasedEEMsearch)EEMsearchTemplate).minSeedGeneSize);	
		tmp.add("itrForPvalue2Calculation= " + ((BiclusterBasedEEMsearch)EEMsearchTemplate).itrForPvalue2Calculation);
		tmp.add("maxItrForISA= " + ((BiclusterBasedEEMsearch)EEMsearchTemplate).maxItrForISA);
		tmp.add("Pvalue1Cutoff= " + ((BiclusterBasedEEMsearch)EEMsearchTemplate).Pvalue1Cutoff);
		return MyFunc.join(" ", tmp); 
	}
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("g", "tg", true,  "Tg parameters (delimited by ':')");
		options.addOption("c", "tc", true,  "Tc parameters (delimited by ':')");
		options.addOption("p", "pcut", true,  "first P value cutoff");
		options.addOption("i", "itr", true,  "itration for second P value calculation");
		options.addOption("a", "absolute", false, "absolute bicluster type");
		options.addOption("G", "gupcdown", false, "GeneUpButConditionDown bicluster type");
		options.addOption("t", "two", false, "two bicluster types");
		options.addOption("A", "all", false, "all bicluster types");
		options.addOption("d", "down", false, "down bicluster type");
		options.addOption("o", "outfile", true, "output file");
		options.addOption("m", "mingeneset", true, "min geneset size");
		options.addOption("M", "maxgeneset", true, "max geneset size");
		options.addOption("e", "expmodset", true, "expression module set file");
		options.addOption("E", "extrapolation", true, "set a pvalue cutoff used for extrapolation");
		options.addOption("l", "log", true, "log file");
		options.addOption("r", "recnull", false, "recycle null distribution");
		options.addOption("S", "maxisa", true,  "maximum itration for ISA");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] expfile gmtfile", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 2){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] expfile gmtfile", options);
			return;
		}
		ParameterOptimizingBiclusterBasedEEMsearch BE = new ParameterOptimizingBiclusterBasedEEMsearch(new MyMat(argList.get(0)));	
		BE.setGeneSets(MyFunc.readGeneSetFromGmtFile(argList.get(1)));
		if(commandLine.hasOption("g")){
			BE.setTgs(MyFunc.toDouble(MyFunc.split(":", (commandLine.getOptionValue("g")))));
		}	
		if(commandLine.hasOption("c")){
			BE.setTcs(MyFunc.toDouble(MyFunc.split(":", (commandLine.getOptionValue("c")))));
		}
		if(commandLine.hasOption("p")){
			BE.setPvalue1Cutoff(Double.valueOf(commandLine.getOptionValue("p")));
		}
		if(commandLine.hasOption("i")){
			BE.setItrForPvalue2Calculation(Integer.valueOf(commandLine.getOptionValue("i")));
		}
		if(commandLine.hasOption("a")){
			BE.setBiclusterType(BiclusterType.ABSOLUTE);
		}
		if(commandLine.hasOption("d")){
			BE.setBiclusterType(BiclusterType.DOWN);
		}
		if(commandLine.hasOption("G")){
			BE.setBiclusterType(BiclusterType.gUPcDOWN);
		}
		if(commandLine.hasOption("t")){
			BE.useTwoBiclusterTypes();
		}
		if(commandLine.hasOption("A")){
			BE.useAllBiclusterTypes();
		}
		if(commandLine.hasOption("m")){
			BE.setMinGeneSetSize(Integer.valueOf(commandLine.getOptionValue("m")));
		}
		if(commandLine.hasOption("M")){
			BE.setMaxGeneSetSize(Integer.valueOf(commandLine.getOptionValue("M")));
		}
		if(commandLine.hasOption("r")){
			BE.recycleNullDistribution();
		}
		if(commandLine.hasOption("S")){
			BE.setMaxItrForISA(Integer.valueOf(commandLine.getOptionValue("S")));
		}
		if(commandLine.hasOption("E")){
			 BE.setExtrapolationCutoff(Double.valueOf(commandLine.getOptionValue("E")));
		}	
		BE.perform();
		ExpressionModuleSet ems = BE.getExpressionModuleSet();
		if(commandLine.hasOption("o")){
			PrintWriter os = new PrintWriter(new FileWriter(commandLine.getOptionValue("o")));
			 for(String s: ems.getIds()){
				 os.println(ems.get(s));				 
			 }
			 os.close();
		}else{
			 for(String s: ems.getIds()){
				 System.out.println(ems.get(s));				 
			 }
			 System.out.close();
		}
		if(commandLine.hasOption("e")){
			ems.writeToFile(commandLine.getOptionValue("e"));
		}
		if(commandLine.hasOption("l")){
			PrintWriter os = new PrintWriter(new FileWriter(commandLine.getOptionValue("l")));
			os.print(BE.getLog());
			os.close();
		}
	}
	
}
