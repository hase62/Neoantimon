package eem;

import java.util.*;
import java.io.*;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

import sun.print.resources.serviceui;
import sun.reflect.Reflection;
import utility.Dist;
import utility.MyException;
import utility.MyFunc;
import utility.MyMat;
import utility.StopWatch;


public class CoherenceBasedEEMsearch extends AbstractEEMsearch{
	MyMat Exp;
	Dist Cor;
	//boolean absoluteCorrelation;
	double absoluteRadius;
	double relativeRadius;
	int coreGeneSize;
	DistConverter distConverter;
	
	interface DistFunc {
		double get(List <Double> a, List <Double> b);
		double get(String a, String b);
	}
	private class Correlation implements DistFunc{
		public double get(List <Double> a, List <Double> b){
			return 1 - MyFunc.pearsonCorrelationForNormarizedList(a,b);
		}
		public double get(String a, String b){
			return 1 - Cor.get(a,b);
		}
	}
	private class AbsoluteCorrelation implements DistFunc{
		public double get(List <Double> a, List <Double> b){
			return 1 - Math.abs(MyFunc.pearsonCorrelationForNormarizedList(a,b));
		}
		public double get(String a, String b){
			return 1 - Math.abs(Cor.get(a,b));
		}
		
		
	}
	DistFunc distfunc; 
	
	public CoherenceBasedEEMsearch(CoherenceBasedEEMsearch C, boolean deep){
		super(C);
		//absoluteCorrelation = C.absoluteCorrelation;
		absoluteRadius = C.absoluteRadius;
		relativeRadius = C.relativeRadius;
		coreGeneSize = C.coreGeneSize;
		Exp = C.Exp;
		Cor = C.Cor;
		if(deep){
			//distfunc = absoluteCorrelation?new AbsoluteCorrelation():new Correlation();
			distfunc = new Correlation();
			distConverter = new DistConverter(this);
		}else{
			distfunc = C.distfunc;
			distConverter = C.distConverter;
		}
	}
	
	public CoherenceBasedEEMsearch(MyMat E){
		super();
		itrForPvalue2Calculation = 300;
		Pvalue1Cutoff = -Math.log10(0.05);
		//absoluteCorrelation = false;
		absoluteRadius = 0.0;
		relativeRadius = 0.0;
		coreGeneSize = 10;
		originalExp = E;
		Exp = new MyMat (E);
		Exp.normalizeRows();
		Cor = null;
		allGenes = Exp.getRowNames(); 
		distfunc =new Correlation();
		distConverter =  new DistConverter(this);
	};
	
	
	protected void setEEM(){
		if(seedGeneSets == null || seedGeneSets.isEmpty()){
			throw new MyException("seedGeneSets must be set!");
		}
		if(Cor == null){
			calculateCor();
		}
		if(absoluteRadius == 0.0){
			setRelativeRadius(0.05);
		}
		
		for(Map.Entry<String,List<String>> s: seedGeneSets.entrySet()){
			EEM e =  new CoherenceBasedEEM(this, s.getValue());
			eems.put(s.getKey(), e);
			seeds.add(s.getKey());
		}
		
	
	}
	public void setCoreGeneSize(int i){
		if(i>3){
			coreGeneSize = i;
		}
	}
	
	
	//public void useAbsoluteCorrelation(){
		//absoluteCorrelation = true;
		//distfunc =new AbsoluteCorrelation();
	//}
	
	//public void useNomalCorrelation(){
		//absoluteCorrelation = false;
		//distfunc =new Correlation();
	//}
	
	public void setCor(Dist D){
		if(!D.getNames().containsAll(Exp.getRowNames())){
			throw new MyException("Correlation data must include data for all genes in expression data");
		}
		Cor = D;
	}
	public void calculateCor(){
		System.err.println("Calculating correlations...");
		Cor = new Dist(Exp, 'C');	
	}
	
	public void setAbsoluteRadius(double r){
		if(Cor == null){
			calculateCor();
		}
		absoluteRadius = r;
		relativeRadius = distConverter.convertAbsolute2relativeDist(r);
	}
	public void setRelativeRadius(double r){
		if(Cor == null){
			calculateCor();
		}
		relativeRadius = r;
		absoluteRadius = distConverter.convertRelative2absoluteDist(r);
	}
	
	
	
	public String toString(){
		List <String> tmp = new ArrayList<String>();
		tmp.add("CoherenceBasedEEM");
		//tmp.add("absoluteCorrelation= " + absoluteCorrelation);
		tmp.add("coreGeneSize=" + coreGeneSize);
		tmp.add("absoluteRadius= " + absoluteRadius);
		tmp.add("relativeRadius= " + relativeRadius);
		tmp.add("maxGeneSetSize= " + maxSeedGeneSize);
		tmp.add("minGeneSetSize= " + minSeedGeneSize);	
		tmp.add("itrForPvalue2Calculation= " + itrForPvalue2Calculation);
		tmp.add("Pvalue1Cutoff= " + Pvalue1Cutoff);
		tmp.add("seeds=" + seeds.size());
		tmp.add("candidates=" + candidates.size());
		return MyFunc.join(" ", tmp); 
	}
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("A", "absrad", true,  "absolute radius");
		options.addOption("R", "relrad", true,  "relative radius");
		options.addOption("C", "cor", true,  "correlation file");
		options.addOption("p", "pcut", true,  "first P value cutoff");
		options.addOption("i", "itr", true,  "itration for second P value calculation");
		options.addOption("a", "absolute", false, "use absolute correlation");
		options.addOption("o", "outfile", true, "output file");
		options.addOption("m", "mingeneset", true, "min geneset size");
		options.addOption("M", "maxgeneset", true, "max geneset size");
		options.addOption("e", "expmodset", true, "expression module set file");
		options.addOption("l", "log", true, "logFile");
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
		CoherenceBasedEEMsearch CE = new CoherenceBasedEEMsearch(new MyMat(argList.get(0)));	
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
		//if(commandLine.hasOption("a")){
		//	CE.useAbsoluteCorrelation();
		//}
		if(commandLine.hasOption("m")){
			CE.setMinGeneSetSize(Integer.valueOf(commandLine.getOptionValue("m")));
		}
		if(commandLine.hasOption("M")){
			CE.setMaxGeneSetSize(Integer.valueOf(commandLine.getOptionValue("M")));
		}
		if(commandLine.hasOption("A")){
			CE.setAbsoluteRadius(Double.valueOf(commandLine.getOptionValue("A")));
		}	
		if(commandLine.hasOption("R")){
			CE.setRelativeRadius(Double.valueOf(commandLine.getOptionValue("R")));
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
	
	
