package eem;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.*;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

import sun.reflect.Reflection;
import utility.*;

public class ParameterOptimizingCoherenceBasedEEMsearch extends AbstractParameterOptimizingEEMsearch{
	private List <Double> relativeRadiuses;
	//private boolean bothCorrelation;
	public ParameterOptimizingCoherenceBasedEEMsearch(MyMat E){
		super();
		EEMsearchTemplate = new CoherenceBasedEEMsearch(E);
		relativeRadiuses = new ArrayList<Double>();
		relativeRadiuses.add(0.05);
		relativeRadiuses.add(0.1);
		relativeRadiuses.add(0.15);
		//bothCorrelation = false;
	}
	
	public void setRelativeRadiuses(List <Double> r){
		relativeRadiuses = new ArrayList<Double>(r);
	}
	
	public void setCoreGeneSize(int i){
		((CoherenceBasedEEMsearch) EEMsearchTemplate).setCoreGeneSize(i);
	}
	/*public void useAbsoluteCorrelation(){
		((CoherenceBasedEEMsearch) EEMsearchTemplate).useAbsoluteCorrelation();
	}
	public void useBothCorrelation(){
		((CoherenceBasedEEMsearch) EEMsearchTemplate).useNomalCorrelation();
	     bothCorrelation = true;
	}*/
	
	public void setCor(Dist D){
		((CoherenceBasedEEMsearch) EEMsearchTemplate).setCor(D);
	}
	protected void initializeEEMsearches(){
		if(((CoherenceBasedEEMsearch) EEMsearchTemplate) != null & ((CoherenceBasedEEMsearch) EEMsearchTemplate).Cor == null){
			((CoherenceBasedEEMsearch) EEMsearchTemplate).calculateCor();
		}
		/*if(bothCorrelation){
			EEMsearch EEMsearchTemplateDeep = new  CoherenceBasedEEMsearch((CoherenceBasedEEMsearch) EEMsearchTemplate,true);
			((CoherenceBasedEEMsearch) EEMsearchTemplateDeep ).useAbsoluteCorrelation();
			for(int i = 0; i < relativeRadiuses.size(); i++){
				CoherenceBasedEEMsearch eem  = new CoherenceBasedEEMsearch((CoherenceBasedEEMsearch) EEMsearchTemplate,false);
				eem.setRelativeRadius(relativeRadiuses.get(i));
				EEMsearches.add(eem);
				CoherenceBasedEEMsearch eem2  = new CoherenceBasedEEMsearch((CoherenceBasedEEMsearch) EEMsearchTemplateDeep,false);
				eem2.setRelativeRadius(relativeRadiuses.get(i));
				EEMsearches.add(eem2);
				
			}
		}else{*/
			for(int i = 0; i < relativeRadiuses.size(); i++){
				CoherenceBasedEEMsearch eem  = new CoherenceBasedEEMsearch((CoherenceBasedEEMsearch) EEMsearchTemplate,false);
				eem.setRelativeRadius(relativeRadiuses.get(i));
				EEMsearches.add(eem);
			}
		//}
		
	}
	public String toString(){
		List <String> tmp = new ArrayList<String>();
		tmp.add("CoherenceBasedEEM");
		//tmp.add("absoluteCorrelation= " + ((CoherenceBasedEEMsearch) EEMsearchTemplate).absoluteCorrelation);
		tmp.add("coreGeneSize=" +  ((CoherenceBasedEEMsearch) EEMsearchTemplate).coreGeneSize);
		tmp.add("relativeRadius= " + relativeRadiuses);
		tmp.add("maxGeneSetSize= " + ((CoherenceBasedEEMsearch) EEMsearchTemplate).maxSeedGeneSize);
		tmp.add("minGeneSetSize= " + ((CoherenceBasedEEMsearch) EEMsearchTemplate).minSeedGeneSize);	
		tmp.add("itrForPvalue2Calculation= " +  ((CoherenceBasedEEMsearch) EEMsearchTemplate).itrForPvalue2Calculation);
		tmp.add("Pvalue1Cutoff= " +  ((CoherenceBasedEEMsearch) EEMsearchTemplate).Pvalue1Cutoff);
		return MyFunc.join(" ", tmp); 
	}

	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("R", "relrad", true,  "relative radiuses (delimited by ':')");
		options.addOption("C", "cor", true,  "correlation file");
		options.addOption("p", "pcut", true,  "first P value cutoff");
		options.addOption("i", "itr", true,  "itration for second P value calculation");
		options.addOption("a", "absolute", false, "use absolute correlation");
		options.addOption("b", "both", false, "use both correlation");
		options.addOption("o", "outfile", true, "out file");
		options.addOption("m", "mingeneset", true, "min geneset size");
		options.addOption("M", "maxgeneset", true, "max geneset size");
		options.addOption("e", "expmodset", true, "expression module set file");
		options.addOption("l", "log", true, "log file");
		options.addOption("r", "recnull", false, "recycle null distribution");
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
		ParameterOptimizingCoherenceBasedEEMsearch  CE = new ParameterOptimizingCoherenceBasedEEMsearch (new MyMat(argList.get(0)));	
		CE.setGeneSets(MyFunc.readGeneSetFromGmtFile(argList.get(1)));
		if(commandLine.hasOption("C")){
			CE.setCor(Dist.readFromBinary(commandLine.getOptionValue("C")));
		}
		if(commandLine.hasOption("c")){
			CE.setCoreGeneSize(Integer.valueOf(commandLine.getOptionValue("c")));
		}
		if(commandLine.hasOption("p")){
			CE.setPvalue1Cutoff(Double.valueOf(commandLine.getOptionValue("p")));
		}
		if(commandLine.hasOption("i")){
			CE.setItrForPvalue2Calculation(Integer.valueOf(commandLine.getOptionValue("i")));
		}
		/*if(commandLine.hasOption("a")){
			CE.useAbsoluteCorrelation();
		}
		if(commandLine.hasOption("b")){
			CE.useBothCorrelation();
		}*/
		if(commandLine.hasOption("m")){
			CE.setMinGeneSetSize(Integer.valueOf(commandLine.getOptionValue("m")));
		}
		if(commandLine.hasOption("M")){
			CE.setMaxGeneSetSize(Integer.valueOf(commandLine.getOptionValue("M")));
		}
		if(commandLine.hasOption("R")){
			CE.setRelativeRadiuses(MyFunc.toDouble(MyFunc.split(":", (commandLine.getOptionValue("R")))));
		}
		if(commandLine.hasOption("r")){
			CE.recycleNullDistribution();
		}
		CE.perform();
		ExpressionModuleSet ems = CE.getExpressionModuleSet();
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
			os.print(CE.getLog());
			os.close();
		}
	}
	
}
