package snp.old;

import java.util.*;
import java.io.*;
import org.apache.commons.cli.*;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.*;

import snp.ProbeInfo;
import snp.Segment;
import snp.SegmentContainer;
import snp.SegmentContainerMap;
import snp.SegmentContainer.SegmentIterator;
import sun.reflect.Reflection;
import utility.MyFunc;
import utility.MyMat;

// this overestimates significance.
public class SAIRICcomparative0 {
	protected SegmentContainerMap  BAF;
	protected SegmentContainerMap  BAF2;

	protected ProbeInfo probeInfo;
	protected List<String>probeList;
	protected List<String>sampleList; 
	
	
	protected Map<String, Double> scores;
	protected Map<String, Double> scores2;
	protected Map<String, Double> pvalues;
	protected Map<String, Double> qvalues;
	
	protected double AICutoff =  0.55;
	

	public SAIRICcomparative0(SegmentContainerMap baf, SegmentContainerMap baf2,  ProbeInfo pi) throws IOException{
		BAF = baf;
		BAF2 = baf2;
		probeInfo = pi; 
		probeList = new ArrayList<String>(probeSet());
		sampleList = new ArrayList<String>(sampleSet());
	}
	
	protected Set<String> sampleSet(){
		return BAF.keySet();
	}
	protected Set<String> sampleSet2(){
		return BAF2.keySet();
	}
	
	protected Set<String> probeSet(){
		return probeInfo.getProbeSet();
	}
	
	protected void calculateScores(){
		scores =  new LinkedHashMap<String, Double>();
		for(String p: probeSet()){
			scores.put(p,0.0);
		}
		
		for(String s: sampleSet()){
			SegmentContainer bafSC = BAF.get(s);
			SegmentContainer.SegmentIterator bafSCitr  = bafSC.iterator();
			Segment S = bafSCitr.next();
			for(String p: probeSet()){
				int chr = probeInfo.chr(p);
				int pos = probeInfo.pos(p);
				while(chr > S.chr()){
					S = bafSCitr.next();
				}
				while(pos > S.end()){
					S = bafSCitr.next();
				}
				
				double baf = S.value();
				if(baf>= AICutoff){
					scores.put(p, scores.get(p)+1);
				
				}
			}
		}
		
		scores2 =  new LinkedHashMap<String, Double>();
		for(String p: probeSet()){
			scores2.put(p,0.0);
		}
		
		for(String s: sampleSet2()){
			SegmentContainer bafSC = BAF2.get(s);
			SegmentContainer.SegmentIterator bafSCitr  = bafSC.iterator();
			Segment S = bafSCitr.next();
			for(String p: probeSet()){
				int chr = probeInfo.chr(p);
				int pos = probeInfo.pos(p);
				while(chr > S.chr()){
					S = bafSCitr.next();
				}
				while(pos > S.end()){
					S = bafSCitr.next();
				}
				
				double baf = S.value();
				if(baf>= AICutoff){
					scores2.put(p, scores2.get(p)+1);
				
				}
			}
		}
		
	}
	
	protected  void calculatePvalues(){
		pvalues =  new LinkedHashMap<String, Double>();
		for(String p: probeSet()){
			int s = (int)(double)scores.get(p);
			int s2 = (int)(double)scores2.get(p);
			double pv =  MyFunc.calculateFisheExactPvalue(s, s2, sampleSet().size()-s, sampleSet2().size()-s2);
			pvalues.put(p, pv);
			 System.err.println(s + " " + s2 + " " + (sampleSet().size()-s) + " " + (sampleSet2().size()-s2) + " " + pv);
		}
	}
	
	protected  void calculateQvalues(){
		qvalues = MyFunc.calculateQvalue(pvalues);
		for(String p: qvalues.keySet()){ 
			qvalues.put(p, qvalues.get(p)>1?1:qvalues.get(p));
		}
	}
	public void perform(){
		System.err.println("calculate scores....");
		calculateScores();
		
		System.err.println("calculate pvalues....");
		calculatePvalues();
		
		System.err.println("calculate qvalues....");
		calculateQvalues();
	}
	
	public Map <String, Double> getQvalues(){
		return qvalues;
	}
	
	public Map <String, Double> getPvalues(){
		return pvalues;
	}
	
	public Map <String, Double> getScores(){
		return scores;
	}
	
	public Map <String, Double> getMinusLogQvalues(){
		Map <String, Double> minusLogQvalues = new LinkedHashMap<String, Double>();
		for(String p: probeList){
		double tmp = qvalues.get(p);
		if(tmp==1){
			tmp=0;
		}else{
			tmp = - Math.log10(tmp);
		}
		minusLogQvalues.put(p, tmp);
		}
		return minusLogQvalues;
	}
	public Map <String, Double> getMinusLogPvalues(){
		Map <String, Double> minusLogPvalues = new LinkedHashMap<String, Double>();
		for(String p: probeList){
		double tmp = pvalues.get(p);
		if(tmp==1){
			tmp=0;
		}else{
			tmp = - Math.log10(tmp);
		}
		minusLogPvalues.put(p, tmp);
		}
		return minusLogPvalues;
	}
	

	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("c", "cutoff", true, "cutoff for AI call");
		options.addOption("o", "outfile", true, "output file name");
		options.addOption("t", "thinpi", true, "thin down probe info");
		options.addOption("i", "interval", true, "interval for psuedo probe info");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] SegFile SegFile2 probeTsvFile", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(!(argList.size() == 2 | argList.size() == 3)){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] SegFile SegFile2 probeTsvFile", options);
			return;
		}
		SegmentContainerMap BAF = new SegmentContainerMap(argList.get(0));
		SegmentContainerMap BAF2 = new SegmentContainerMap(argList.get(1));
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
		
		SAIRICcomparative0 SAIRIC = new SAIRICcomparative0(BAF,BAF2,PI);
		if(commandLine.hasOption("c")){
			SAIRIC.AICutoff = Double.valueOf(commandLine.getOptionValue("c"));
		}
		
		
		SAIRIC.perform();
		
		Writer os;
		if(commandLine.hasOption("o")){
			 os = new BufferedWriter(new FileWriter(commandLine.getOptionValue("o")));	 
		}else{
			os = new PrintWriter(System.out);
		}
			os.write("Probe" + "\t" +  "Chrom" + "\t" +  "BasePair"  + "\t" + "pvalue" + "\t" + "qvalue" + "\n");
		for(String s: PI.getProbeSet()){
			os.write(s + "\t" +  PI.chr(s) + "\t" +  PI.pos(s) + "\t" + SAIRIC.getMinusLogPvalues().get(s) + "\t" + SAIRIC.getMinusLogQvalues().get(s) + "\n");			
		}
		os.flush();
	}
	
	
	
	

}
