package snp.old;

import java.io.*;
import java.util.*;

import org.apache.commons.cli.*;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.*;

import snp.ProbeInfo;
import snp.SegmentContainerMap;
import sun.reflect.Reflection;
import utility.*;

public class SAIRIC2 extends SAIRIC {

	protected List <Double> binomialParameterP;
	
	protected List <Distribution> nullDistributions;
	
	protected double normalApproxCutoff = 5;
	
	protected int nullDistSize = 100000;
	
	protected MyMat zscores;
	
	protected boolean useZscore = false;
	
	public SAIRIC2(SegmentContainerMap baf,SegmentContainerMap logr, ProbeInfo pi)throws IOException {
		super(baf, logr, pi);
	}

	
	private void calculateBinomialParameterP(){
		binomialParameterP = new ArrayList<Double>();
		for(int i = 0; i < scoreTypes.size();i++){
			binomialParameterP.add(0.0);
		}
		for(int i = 0; i < nullDistSize; i++){
			String p = generateRandomProbes();
			String s = generateRandomSamples();
			int chr = probeInfo.chr(p);
			int pos = probeInfo.pos(p);
			double baf = BAF.get(s).get(chr, pos).value();
			double logr = LogR.get(s).get(chr, pos).value();
			//amp
			if(logr > ampCutoff){
				binomialParameterP.set(0,binomialParameterP.get(0)+1);
			}
			//del
			if(logr < delCutoff){
				binomialParameterP.set(1,binomialParameterP.get(1)+1);
			}
			//neut
			if(logr >= delCutoff & logr <= ampCutoff){
				binomialParameterP.set(2,binomialParameterP.get(2)+1);
			}	
			//AI
			if(baf>= AICutoff){
				binomialParameterP.set(3,binomialParameterP.get(3)+1);
			}
			//ampAI
			if(baf>= AICutoff & logr > ampCutoff){
				binomialParameterP.set(4,binomialParameterP.get(4)+1);
			}
			//delAI
			if(baf>= AICutoff & logr < delCutoff){
				binomialParameterP.set(5,binomialParameterP.get(5)+1);
			}
			//neutAI
			if(baf>= AICutoff & logr >= delCutoff & logr <= ampCutoff){
				binomialParameterP.set(6,binomialParameterP.get(6)+1);
			}
		}
		for(int i = 0; i < scoreTypes.size();i++){
			binomialParameterP.set(i,binomialParameterP.get(i)/ nullDistSize);
		}	
	}

	protected void prepearNullDistributions(){
		nullDistributions = new ArrayList<Distribution>();
		calculateBinomialParameterP();
		for(int i = 0; i < scoreTypes.size();i++){
			double p = binomialParameterP.get(i);
			int n = sampleList.size();
			Distribution D;
			if(n*p > normalApproxCutoff & n*(1-p)> normalApproxCutoff){
				D = new NormalDistributionImpl(n*p, Math.pow(n*p*(1-p),0.5));
			}else{
				D = new BinomialDistributionImpl(n,p);
			}
			nullDistributions.add(D);
		}
	}
	
	protected  void calculatePvalues(){
		pvalues =  new MyMat(probeList, scoreTypes);
		for(int i = 0; i < scoreTypes.size(); i++){
			for(String p: probeList){
				double score = scores.get(p, scoreTypes.get(i));
				double pvalue = 1;
			    try {
					pvalue = 1- nullDistributions.get(i).cumulativeProbability(score);
				} catch (MathException e) {
						e.printStackTrace();
				} 
				pvalues.set(p, scoreTypes.get(i), pvalue);
			}
		}
	}
	
	
	protected  void calculateZscore(){
		zscores =  new MyMat(probeList, scoreTypes);
		calculateBinomialParameterP();
		for(int i = 0; i < scoreTypes.size();i++){
			double p = binomialParameterP.get(i);
			int n = sampleList.size();
			double mean = n*p;
			double sd = Math.pow(n*p*(1-p),0.5);
			for(String probe: probeList){
				double score = scores.get(probe, scoreTypes.get(i));
				double z = (score-mean)/sd;
				zscores.set(probe, scoreTypes.get(i), Math.abs(z));
			}
		}
	}
	
	public void perform(){
		System.err.println("calculate scores....");
		calculateScores();
		System.err.println("prepare null distributions....");
		if(useZscore){
			System.err.println("calculate zscores....");
			calculateZscore();
		}else{
			prepearNullDistributions();
			System.err.println("calculate pvalues....");
			calculatePvalues();
			System.err.println("calculate qvalues....");
			calculateQvalues();
		}
	}
	
	public MyMat getZscores(){
		return zscores;
	}
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("a", "ampcutoff", true, "cutoff for amplification call");
		options.addOption("d", "delcutoff", true, "cutoff for deletion call");
		options.addOption("A", "aicutoff", true, "cutoff for AI call");
		options.addOption("n", "ndsize", true, "size of null distributions");
		options.addOption("o", "outfile", true, "output file name");
		options.addOption("t", "thinpi", true, "thin down probe info");
		options.addOption("i", "interval", true, "interval for psuedo probe info");
		options.addOption("N", "normap", true, "cutoff for normal approximation");
		options.addOption("z", "zscore", false, "use zscore");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] BAFSegFile LogRSegFile probeTsvFile", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(!(argList.size() == 3 | argList.size() == 2)){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] BAFSegFile LogRSegFile probeTsvFile", options);
			return;
		}
		SegmentContainerMap BAF = new SegmentContainerMap(argList.get(0));
		SegmentContainerMap LogR = new SegmentContainerMap(argList.get(1));
		ProbeInfo PI;
		List<String> sample = new ArrayList<String>(BAF.keySet());
		if(argList.size() == 3 ){
			PI = 	ProbeInfo.getProbeInfoFromTsvFile(argList.get(2));
			if(commandLine.hasOption("t")){
				PI.thinDown(Integer.valueOf(commandLine.getOptionValue("t")));
			}
			PI.filter(BAF.get(sample.get(0)));
		}else{
			int interval = 10000000;
			if(commandLine.hasOption("i")){
				interval = Integer.valueOf(commandLine.getOptionValue("i"));
			}
			PI = ProbeInfo.generatePsuedoProbeInfo(interval);
			PI.filter(BAF.get(sample.get(0)));
			}
		
		SAIRIC2 SAIRIC = new SAIRIC2(BAF,LogR,PI);
		if(commandLine.hasOption("a")){
			SAIRIC.ampCutoff = Double.valueOf(commandLine.getOptionValue("a"));
		}
		if(commandLine.hasOption("d")){
			SAIRIC.delCutoff = Double.valueOf(commandLine.getOptionValue("d"));
		}
		if(commandLine.hasOption("A")){
			SAIRIC.AICutoff = Double.valueOf(commandLine.getOptionValue("A"));
		}
		if(commandLine.hasOption("p")){
			SAIRIC.nullDistSize = Integer.valueOf(commandLine.getOptionValue("p"));
		}
		if(commandLine.hasOption("N")){
			SAIRIC.normalApproxCutoff = Double.valueOf(commandLine.getOptionValue("N"));
		}
		if(commandLine.hasOption("z")){
			SAIRIC.useZscore = true;
		}
		
		SAIRIC.perform();
		MyMat tmp;
		if(SAIRIC.useZscore){
			tmp = SAIRIC.getZscores(); 
		}else{
			tmp = SAIRIC.getMinusLogQvalues();
		}
		PrintWriter os;
		if(commandLine.hasOption("o")){
			 os = new PrintWriter(new FileWriter(commandLine.getOptionValue("o")));
		}else{
			os = new PrintWriter(System.out);
		}
		os.println("Probe" + "\t" +  "Chrom" + "\t" +  "BasePair" + "\t" + MyFunc.join("\t", SAIRIC.scoreTypes));
		for(String s: PI.getProbeSet()){
			os.println(s + "\t" +  PI.chr(s) + "\t" +  PI.pos(s) + "\t" + MyFunc.join("\t", tmp.getRow(s)));			
		}
		os.flush();
	}
	
	
	
	
	
	
}
