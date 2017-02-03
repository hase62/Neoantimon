package eem;

import java.util.*;
import java.util.zip.DataFormatException;
import java.io.*;

import utility.LinearRegression;
import utility.MyException;
import utility.MyFunc;
import utility.MyMat;
import org.apache.commons.cli.*;

public class BiclusterBasedEEMsearch extends AbstractEEMsearch {
	MyMat Eg; //expression normalized for gene
	MyMat Ec; //expression normalized for conditions/samples
	double Tg; //threshold for genes
	double Tc; // threshold for conditions
	int maxItrForISA;
	Double extrapolationCutoff = null;
	CoherenceBasedEEMsearch coEEM;
	public static class BiclusterType implements Serializable{
		private static final long serialVersionUID = 3446392147178993187L;
		private final String name;
		private BiclusterType(String name){this.name = name;}
		public BiclusterType(BiclusterType b){this.name = b.name;}
		public String toString(){return name;}
		public boolean equals(BiclusterType b){return name.equals(b.toString());};
		public static final BiclusterType UP = new BiclusterType("Up");
		public static final BiclusterType DOWN = new BiclusterType("Down");
		public static final BiclusterType ABSOLUTE = new BiclusterType("Absolute");
		public static final BiclusterType gUPcDOWN  = new BiclusterType("GeneUpButConditionDown");
	}
	BiclusterType biclusterType; 
	
	public BiclusterBasedEEMsearch(BiclusterBasedEEMsearch B, boolean deep){
		super(B);
		if(deep){
			Eg = new MyMat(B.Eg);
			Ec = new MyMat(B.Ec);
			Tg = B.Tg;
			Tc = B.Tc;
			maxItrForISA = B.maxItrForISA;
			biclusterType = new BiclusterType(B.biclusterType);
		}else{
			Eg = B.Eg;
			Ec = B.Ec;
			Tg = B.Tg;
			Tc = B.Tc;
			maxItrForISA = B.maxItrForISA;
			biclusterType = B.biclusterType;
		}
	}
	
	public  BiclusterBasedEEMsearch(MyMat Exp){
		super();
		itrForPvalue2Calculation = 10000;
		Pvalue1Cutoff = 2.0;
		this.originalExp = Exp;
		Eg = new MyMat(Exp);
		Ec = new MyMat(Exp);
		Eg.normalizeCols();
		Ec.normalizeRows();
		allGenes = new ArrayList<String>(Exp.getRowNames());
		Tg = 0.1;
		Tc = 0.1;
		maxItrForISA = 30;
		biclusterType = BiclusterType.UP;
	}
	
	protected void initializedCoherenceBasedEEM(){
		coEEM = new CoherenceBasedEEMsearch(Ec);
		coEEM.calculateCor();
		coEEM.setRelativeRadius(0.15);
	}
	
	protected void setEEM(){
		if(seedGeneSets == null || seedGeneSets.isEmpty()){
			throw new MyException("seedGeneSets must be set!");
		}
	
		for(Map.Entry<String,List<String>> s: seedGeneSets.entrySet()){
			EEM e = new BiclusterBasedEEM(this, s.getValue());
			eems.put(s.getKey(), e);
			seeds.add(s.getKey());
		}
		
		
	}
	public void setBiclusterType(BiclusterType type){
		if(type.equals(BiclusterType.UP)){
			setBiclusterType2Up();
			return;	
		}
		if(type.equals(BiclusterType.gUPcDOWN)){
			setBiclusterType2GeneUpButConditionDown();
			return;
		}	
		if(type.equals(BiclusterType.DOWN)){
			setBiclusterType2Down();
			return;	
		}
		if(type.equals(BiclusterType.ABSOLUTE)){
			setBiclusterType2Absolute();
			return;
		}			
	}
	public void setBiclusterType2Up(){
		Eg = new MyMat(originalExp);
		Ec = new MyMat(originalExp);
		Eg.normalizeCols();
		Ec.normalizeRows();
		biclusterType = BiclusterType.UP;
	}
    public void setBiclusterType2Down(){
    	if(!biclusterType.equals(BiclusterType.UP)){
			setBiclusterType2Up();
		}	
		int i,j;
		for(i=0;i<Ec.rowSize();i++){
			for(j=0;j<Ec.colSize();j++){
				Ec.set(i, j, -Ec.get(i, j));
			}	
		}
		for(i=0;i<Eg.rowSize();i++){
			for(j=0;j<Eg.colSize();j++){
				Eg.set(i, j, -Eg.get(i, j));
			}	
		}
		biclusterType = BiclusterType.DOWN;
	}
    public void setBiclusterType2Absolute(){
		if(!biclusterType.equals(BiclusterType.UP)){
			setBiclusterType2Up();
		}	
		int i,j;
		for(i=0;i<Ec.rowSize();i++){
			for(j=0;j<Ec.colSize();j++){
				Ec.set(i, j, Math.abs(Ec.get(i, j)));
			}	
		}
		for(i=0;i<Eg.rowSize();i++){
			for(j=0;j<Eg.colSize();j++){
				Eg.set(i, j, Math.abs(Eg.get(i, j)));
			}	
		}
		biclusterType = BiclusterType.ABSOLUTE;
	}
	
    public void setBiclusterType2GeneUpButConditionDown(){
		if(!biclusterType.equals(BiclusterType.UP)){
			setBiclusterType2Up();
		}	
		int i,j;
		for(i=0;i<Ec.rowSize();i++){
			for(j=0;j<Ec.colSize();j++){
				Ec.set(i, j, -Ec.get(i, j));
			}	
		}
		biclusterType = BiclusterType.gUPcDOWN;		
	}
	
	
	
	String errLog = ""; 
	
	public void setMaxItrForISA(int i){
		if(i > 30){
			maxItrForISA = i;
		}
	}
	
	
	public void setTc(double d){
		if(d > 1){
			Tc = d/100;
		}else{
			Tc = d;
		}	
	}
	
	public void setTg(double d){
		if(d > 1){
			Tg = d/100;
		}else{
			Tg = d;
		}
		
	}
	
	
	public void setExtrapolationCutoff(double d){
		extrapolationCutoff =d;
	}
	
	public void extrapolateSmallPvalues(){
		System.err.println("Extrapolating small Pvalues...");
		List<Double> trainingX = new ArrayList<Double>();
		List <Double>  trainingY = new ArrayList<Double>() ;
		for(String s: candidates){
			if(eems.get(s).getPvalue1() > extrapolationCutoff && eems.get(s).getPvalue2() != null && eems.get(s).getPvalue2()!=Double.MAX_VALUE){
				trainingX.add(eems.get(s).getPvalue1());
				trainingY.add(eems.get(s).getPvalue2());
			}
		}
		
		if(trainingX.size() < 10){
			throw new MyException("training data size ("+ trainingX.size() +")is too small!");	
		}
		
		LinearRegression L = new LinearRegression(trainingX, trainingY);	
		L.learn();
		System.err.println(L);
		
		for(String s: candidates){
			if(eems.get(s).getPvalue2()==Double.MAX_VALUE){
				double p = L.predict(eems.get(s).getPvalue1());
				eems.get(s).setPvalue2(p);
				eems.get(s).setPvalue(p);
			}
		}
	}
	
	public void perform(){
		try{
			setEEM();
			findModuleGenes();
			calculatePvalue1();
			findCandidates();
			calculatePvalue2();
			if(extrapolationCutoff != null){
				try{	
					extrapolateSmallPvalues();
				} catch (Exception e2) {
					System.err.println(e2.getMessage());
				}
			}
		} catch (Exception e) {
			System.err.println(e.getMessage());
			e.printStackTrace();
		}
		finally{
			stopWatch.stop();
			writeTimeLog();
		}
	}
	
	public String toString(){
		List <String> tmp = new ArrayList<String>();
		tmp.add("biclusterType= " + biclusterType);
		tmp.add("Tg= " + Tg);
		tmp.add("Tc= " + Tc);
		tmp.add("maxGeneSetSize= " + maxSeedGeneSize);
		tmp.add("minGeneSetSize= " + minSeedGeneSize);	
		tmp.add("itrForPvalue2Calculation= " + itrForPvalue2Calculation);
		tmp.add("maxItrForISA= " + maxItrForISA);
		tmp.add("Pvalue1Cutoff= " + Pvalue1Cutoff);
		tmp.add("seeds=" + seeds.size());
		tmp.add("candidates=" + candidates.size());
		return MyFunc.join(" ", tmp); 
	}
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("g", "tg", true,  "Tg parameter");
		options.addOption("c", "tc", true,  "Tc parameter");
		options.addOption("p", "cut", true,  "first P value cutoff");
		options.addOption("i", "itr", true,  "itration for second P value calculation");
		options.addOption("a", "absolute", false, "absolute bicluster type");
		options.addOption("d", "down", false, "down bicluster type");
		options.addOption("G", "gupcdown", false, "GeneUpButConditionDown bicluster type");
		options.addOption("o", "outfile", true, "output file");
		options.addOption("m", "mingeneset", true, "min geneset size");
		options.addOption("M", "maxgeneset", true, "max geneset size");
		options.addOption("e", "expmodset", true, "expression module set file");
		options.addOption("E", "extrapolation", true, "set a pvalue cutoff used for extrapolation");
		options.addOption("l", "log", true, "logfile");
		options.addOption("r", "recnull", false, "recycle nulldistribution");
		options.addOption("S", "maxisa", true,  "maximum itration for ISA");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp("BiclusterBasedEEMsearch [options] expfile gmtfile", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 2){
			formatter.printHelp("BiclusterBasedEEMsearch [options] expfile gmtfile", options);
			return;
		}
		BiclusterBasedEEMsearch BE = new BiclusterBasedEEMsearch(new MyMat(argList.get(0)));	
		BE.setGeneSets(MyFunc.readGeneSetFromGmtFile(argList.get(1)));
		
		if(commandLine.hasOption("g")){
			BE.setTg(Double.valueOf(commandLine.getOptionValue("g")));
		}	
		if(commandLine.hasOption("c")){
			BE.setTc(Double.valueOf(commandLine.getOptionValue("c")));
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
		
		BE.initializedCoherenceBasedEEM();
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
