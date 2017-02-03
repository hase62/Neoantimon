package snp.old;

import java.util.*;
import java.io.*;
import java.lang.*;

import org.apache.commons.cli.*;

import utility.MyFunc;
import utility.MyMat;

import snp.ProbeInfo;
import snp.Segment;
import snp.SegmentContainer;
import snp.SegmentContainerMap;
import snp.SegmentContainer.SegmentIterator;
import sun.reflect.Reflection;


// Selective Allelic Imbalance Region Identification in Cancer 

public class SAIRIC {
		
	protected SegmentContainerMap  BAF; 
	protected SegmentContainerMap  LogR;
	protected ProbeInfo probeInfo;
	
	protected List<String>probeList;
	protected List<String>sampleList; 
	
	protected MyMat scores;
	protected MyMat pvalues;
	protected MyMat qvalues;
	protected ArrayList<String> scoreTypes; // amp, del, neut, AI, ampAI, delAI, neutAI;
	
	protected MyMat nullScores;
	
	protected double AICutoff =  0.55;
	protected double ampCutoff = 0.1;
	protected double delCutoff = -0.1;
	
	protected int nullDistSize = 100000;

	
	public SAIRIC(SegmentContainerMap baf, SegmentContainerMap logr, ProbeInfo pi) throws IOException{
		BAF = baf;
		LogR = logr;
		probeInfo = pi; 
		scoreTypes = new ArrayList<String>();
		scoreTypes.add("amp");
		scoreTypes.add("del");
		scoreTypes.add("neut");
		scoreTypes.add("AI");
		scoreTypes.add("ampAI");
		scoreTypes.add("delAI");
		scoreTypes.add("neutAI");
		probeList = new ArrayList<String>(probeSet());
		sampleList = new ArrayList<String>(sampleSet());
	}
	
	protected Set<String> sampleSet(){
		return BAF.keySet();
	}
	
	protected Set<String> probeSet(){
		return probeInfo.getProbeSet();
	}
		
	protected void calculateScores(){
		scores =  new MyMat(probeList, scoreTypes);
		for(String s: sampleSet()){
			SegmentContainer bafSC = BAF.get(s);
			SegmentContainer logrSC = LogR.get(s);
			SegmentContainer.SegmentIterator bafSCitr  = bafSC.iterator();
			SegmentContainer.SegmentIterator logrSCitr  = logrSC.iterator();
			Segment S = bafSCitr.next();
			Segment S2 = logrSCitr.next();
			for(String p: probeSet()){
				int chr = probeInfo.chr(p);
				int pos = probeInfo.pos(p);
				while(chr > S.chr()){
					S = bafSCitr.next();
				}
				while(pos > S.end()){
					S = bafSCitr.next();
				}
				
				while(chr > S2.chr()){
					S2 = bafSCitr.next();
				}
				while(pos > S2.end()){
					S2 = bafSCitr.next();
				}
				
				double baf = S.value();
				double logr = S2.value();
				
				
				//amp
				if(logr > ampCutoff){
					scores.set(p, scoreTypes.get(0), scores.get(p, scoreTypes.get(0))+1);
				}
				//del
				if(logr < delCutoff){
					scores.set(p, scoreTypes.get(1), scores.get(p, scoreTypes.get(1))+1);
				}
				//neut
				if(logr >= delCutoff & logr <= ampCutoff){
					scores.set(p, scoreTypes.get(2), scores.get(p, scoreTypes.get(2))+1);
				}	
				//AI
				if(baf>= AICutoff){
					scores.set(p, scoreTypes.get(3), scores.get(p, scoreTypes.get(3))+1);
				}
				//ampAI
				if(baf>= AICutoff & logr > ampCutoff){
					scores.set(p, scoreTypes.get(4), scores.get(p, scoreTypes.get(4))+1);
				}
				//delAI
				if(baf>= AICutoff & logr < delCutoff){
					scores.set(p, scoreTypes.get(5), scores.get(p, scoreTypes.get(5))+1);
				}
				//neutAI
				if(baf>= AICutoff & logr >= delCutoff & logr <= ampCutoff){
					scores.set(p, scoreTypes.get(6), scores.get(p, scoreTypes.get(6))+1);
				}
			}
		}
	}
	
	
	protected String generateRandomProbes(){
		int i = (int) Math.floor(Math.random()*probeList.size());
		return probeList.get(i);
	}
	
	protected String generateRandomSamples(){
		int i = (int) Math.floor(Math.random()*sampleList.size());
		return sampleList.get(i);
	}
	
	protected  void calculateNullScores(){
		List <String> tmp = new ArrayList<String>();
		for(int i = 0; i<nullDistSize; i++){
			tmp.add("null" + i);
		}
		nullScores =  new MyMat(tmp, scoreTypes);
		for(String s: sampleSet()){
			SegmentContainer bafSC = BAF.get(s);
			SegmentContainer logrSC = LogR.get(s);
			for(String p:tmp){
				String randProbe = generateRandomProbes();
				int chr = probeInfo.chr(randProbe);
				int pos = probeInfo.pos(randProbe);
				double baf = bafSC.get(chr, pos).value();
				double logr = logrSC.get(chr, pos).value();
				//amp
				if(logr > ampCutoff){
					nullScores.set(p, scoreTypes.get(0), nullScores.get(p, scoreTypes.get(0))+1);
				}
				//del
				if(logr < delCutoff){
					nullScores.set(p, scoreTypes.get(1), nullScores.get(p, scoreTypes.get(1))+1);
				}
				//neut
				if(logr >= delCutoff & logr <= ampCutoff){
					nullScores.set(p, scoreTypes.get(2), nullScores.get(p, scoreTypes.get(2))+1);
				}	
				//AI
				if(baf>= AICutoff){
					nullScores.set(p, scoreTypes.get(3), nullScores.get(p, scoreTypes.get(3))+1);
				}
				//ampAI
				if(baf>= AICutoff & logr > ampCutoff){
					nullScores.set(p, scoreTypes.get(4), nullScores.get(p, scoreTypes.get(4))+1);
				}
				//delAI
				if(baf>= AICutoff & logr < delCutoff){
					nullScores.set(p, scoreTypes.get(5), nullScores.get(p, scoreTypes.get(5))+1);
				}
				//neutAI
				if(baf>= AICutoff & logr >= delCutoff & logr <= ampCutoff){
					nullScores.set(p, scoreTypes.get(6), nullScores.get(p, scoreTypes.get(6))+1);
				}
			}
		}
	}
		
	
	protected  void calculatePvalues(){
		pvalues =  new MyMat(probeList, scoreTypes);
		for(String s: scoreTypes){
			List <Double> nullScoresTmp = nullScores.getCol(s);
			Collections.sort(nullScoresTmp);
			Collections.reverse(nullScoresTmp);
			for(String p: probeList){
				double score = scores.get(p, s);
				double i;
				for(i=0;i<nullScoresTmp.size() && nullScoresTmp.get((int)i) > score;i++){}
				if(i==0){
					i=1;
				}
				double P = i / nullScoresTmp.size();
				pvalues.set(p, s, P);		
			}
		}
	}
	
	protected  void calculateQvalues(){
		qvalues =  new MyMat(probeList, scoreTypes);
		for(String s: scoreTypes){
			Map <String, Double> pmap = pvalues.getColMap(s);
			Map <String, Double> qmap = MyFunc.calculateQvalue(pmap);
			for(String p: qmap.keySet()){ 
				qvalues.set(p, s, qmap.get(p)>1?1:qmap.get(p));
			}
		}
	}
	
	public void perform(){
		System.err.println("calculate scores....");
		calculateScores();
		System.err.println("calculate null scores....");
		calculateNullScores();
		System.err.println("calculate pvalues....");
		calculatePvalues();
		System.err.println("calculate qvalues....");
		calculateQvalues();
	}
	
	public MyMat getNullScore(){
		return nullScores;
	}
	
	public MyMat getQvalues(){
		return qvalues;
	}
	
	public MyMat getPvalues(){
		return pvalues;
	}
	
	public MyMat getScores(){
		return scores;
	}
	public MyMat getMinusLogQvalues(){
		MyMat minusLogQvalues = new MyMat(probeList, scoreTypes);
		for(String s: scoreTypes){
			for(String p: probeList){
				double tmp = qvalues.get(p, s);
				if(tmp==1){
					tmp=0;
				}else{
					tmp = - Math.log10(tmp);
				}
				minusLogQvalues.set(p, s, tmp);
			}
		}
		return minusLogQvalues;
	}
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("a", "ampcutoff", true, "cutoff for amplification call");
		options.addOption("d", "delcutoff", true, "cutoff for deletion call");
		options.addOption("A", "aicutoff", true, "cutoff for AI call");
		options.addOption("n", "ndsize", true, "size of null distributions");
		options.addOption("o", "outfile", true, "output file name");
		options.addOption("l", "logp", false, "get minus log pvalues");
		options.addOption("t", "thinpi", true, "thin down probe info");
		options.addOption("i", "interval", true, "interval for psuedo probe info");
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
		
		SAIRIC SAIRIC = new SAIRIC(BAF,LogR,PI);
		if(commandLine.hasOption("a")){
			SAIRIC.ampCutoff = Double.valueOf(commandLine.getOptionValue("a"));
		}
		if(commandLine.hasOption("d")){
			SAIRIC.delCutoff = Double.valueOf(commandLine.getOptionValue("d"));
		}
		if(commandLine.hasOption("A")){
			SAIRIC.AICutoff = Double.valueOf(commandLine.getOptionValue("A"));
		}
		if(commandLine.hasOption("n")){
			SAIRIC.nullDistSize = Integer.valueOf(commandLine.getOptionValue("n"));
		}
		
		SAIRIC.perform();
		MyMat tmp;
		if(commandLine.hasOption("l")){
			tmp = SAIRIC.getMinusLogQvalues();
		}else{
			tmp = SAIRIC.getQvalues();
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
