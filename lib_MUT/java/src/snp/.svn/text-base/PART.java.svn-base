package snp;

import java.io.*;
import java.util.*;
import java.util.zip.DataFormatException;

import utility.*;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.*;
import org.apache.commons.cli.*;
import sun.reflect.Reflection;


public class PART {
	
	protected SegmentContainerMap  BAF; 
	protected SegmentContainerMap  LogR;
	protected ProbeInfo probeInfo;
	
	protected List<String>probeList;
	protected List<String>sampleList;
	protected Set<String>sampleSet; 
	
	protected MyMat scores;
	protected MyMat pvalues;
	protected MyMat qvalues;
	
	protected String[] scoreTypesArray = {"amp", "del", "neut", "AI", "ampAI", "delAI", "neutAI", "LOH"};
	protected ArrayList<String> scoreTypes; // amp, del, neut, AI, ampAI, delAI, neutAI, LOH;
	protected int N = scoreTypesArray.length;
	
	protected MyMat nullScores;
	
	protected double AICutoff =  0.65;
	protected double ampCutoff = 0.1;
	protected double delCutoff = -0.1;
	
	protected MyMat cnStatus;
	
	protected MyMat poissonBinomialP; //[amp, del, neut, AI, ampAI, delAI, neutAI] x sample 
	
	protected List<Condition> ConditionFunctions; 
	
	interface Condition {
		boolean check(double logr, double baf);
	}
	
	private void initializeConditionFunctions(){
		ConditionFunctions = new ArrayList<Condition>();
		//amp
		ConditionFunctions.add(
				new Condition(){
					public boolean check(double logr, double baf){
						return logr > ampCutoff;
					}
				}
			);
		//del
		ConditionFunctions.add(
				new Condition(){
					public boolean check(double logr, double baf){
						return logr < delCutoff;
					}
				}
			);
		//neut
		ConditionFunctions.add(
				new Condition(){
					public boolean check(double logr, double baf){
						return logr >= delCutoff & logr <= ampCutoff;
					}
				}
			);
		//AI
		ConditionFunctions.add(
				new Condition(){
					public boolean check(double logr, double baf){
						return baf>= AICutoff;
					}
				}
			);
		//ampAI
		ConditionFunctions.add(
				new Condition(){
					public boolean check(double logr, double baf){
						return baf>= AICutoff & logr > ampCutoff;
					}
				}
			);
		//delAI
		ConditionFunctions.add(
				new Condition(){
					public boolean check(double logr, double baf){
						return baf>= AICutoff & logr < delCutoff;
					}
				}
			);
		//neutAI
		ConditionFunctions.add(
				new Condition(){
					public boolean check(double logr, double baf){
						return baf>= AICutoff & logr >= delCutoff & logr <= ampCutoff;
					}
				}
			);
		//LOH
		ConditionFunctions.add(
				new Condition(){
					public boolean check(double logr, double baf){
						return baf>= AICutoff & logr <= ampCutoff;
					}
				}
			);
	}
	
	public PART(SegmentContainerMap baf,SegmentContainerMap logr, ProbeInfo pi)throws IOException, DataFormatException {
		BAF = baf;
		LogR = logr;
		if(!BAF.checkConsistency()){
			System.err.println("ERR: in SegmentContainerMap consistensiy");
			throw new DataFormatException("");
		}
		if(!LogR.checkConsistency()){
			System.err.println("ERR: in SegmentContainerMap consistensiy");
			throw new DataFormatException("");
		}
		probeInfo = pi; 
		scoreTypes = new ArrayList<String>();
		for(int i = 0; i < N; i++){
			scoreTypes.add(scoreTypesArray[i]);
		}
		sampleSet = new HashSet<String>();
		for(String s: LogR.keySet()){
			if(BAF.containsKey(s) & LogR.containsKey(s)){
				sampleSet.add(s);
			}
		}
		probeList = new ArrayList<String>(probeSet());
		sampleList = new ArrayList<String>(sampleSet());
		initializeConditionFunctions();
	}
	
	
	protected Set<String> sampleSet(){
		return sampleSet;
	}
	
	protected Set<String> probeSet(){
		return probeInfo.getProbeSet();
	}
	
	
	
	protected void calculateScores(){
		scores =  new MyMat(probeList, scoreTypes);
		MyMat CountInSample = new MyMat(new ArrayList<String>(sampleSet()), scoreTypes);
		
		cnStatus = new MyMat(new ArrayList<String>(probeSet()),new ArrayList<String>(sampleSet()));
		
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
					S2 = logrSCitr.next();
				}
				while(pos > S2.end()){
					S2 = logrSCitr.next();
				}
				
				double baf = S.value();
				double logr = S2.value();
				
				for(int i = 0; i < N; i++){
					if(ConditionFunctions.get(i).check(logr, baf)){
						scores.set(p, scoreTypes.get(i), scores.get(p, scoreTypes.get(i))+1);
						CountInSample.set(s, scoreTypes.get(i), CountInSample.get(s, scoreTypes.get(i))+1);
						cnStatus.set(p, s, i);
					}
					}
			}
		}
		
		poissonBinomialP = new MyMat(scoreTypes,sampleList);
		
		for(String s1: scoreTypes){
			for(String s2: sampleSet()){
				poissonBinomialP.set(s1, s2, CountInSample.get(s2,s1)/probeSet().size());
			}
		}			
	}
	
	
	protected void showAberrationRatio(){
		for(String s: scoreTypes){
			double m  = MyFunc.mean(poissonBinomialP.getRow(s));
			System.err.println(s + ": " + m);
		}
	}
	

	public  void calculatePvalues(){
		pvalues =  new MyMat(probeList, scoreTypes);
		
		for(String s: scoreTypes){
			List <Double> tmp = scores.getCol(s);
			double M = MyFunc.max(tmp);
			double m = MyFunc.min(tmp);
			List <Integer> kk = new ArrayList<Integer>();
			for(int i=(int) m;i<=(int)M;i++){
				kk.add(i);
			}
			List<Double> pp = poissonBinomialP.getRow(s);
			List<Double> pv = getPoissonBinomialPvalue(kk,pp);
			
			Map<Integer, Double> kk2pv = new HashMap<Integer, Double>();
			for(int i=0; i<kk.size();i++){
				kk2pv.put(kk.get(i),pv.get(i));
			}
			
			for(String p: probeList){
				int score = (int)scores.get(p, s);
				double pvalue = kk2pv.get(score);
				pvalues.set(p, s, pvalue);
			}
		}
	}
	
	
	public static List<Double> getPoissonBinomialPvalue(List<Integer> kk, List<Double> pp){
		   String tmpFile = "tmp" + Math.round(Math.random()*100000000); 
		   try {
			   PrintWriter os;
			   os = new PrintWriter(new FileWriter(tmpFile + ".kk"));
			   for(int i = 0, n = kk.size(); i < n; i++){
				   os.println(kk.get(i).toString());
			   }
			   os.flush();
			   os.close();
			   os = new PrintWriter(new FileWriter(tmpFile + ".pp"));
			   for(int i = 0, n = pp.size(); i < n; i++){
				   os.println(pp.get(i).toString());
			   }
			   os.flush();
			   os.close();
		   
			   os = new PrintWriter(new FileWriter(tmpFile + ".R"));
			   os.println("kk<-scan(\"" +  tmpFile + ".kk" + "\")");
			   os.println("pp<-scan(\"" +  tmpFile + ".pp" + "\")");
			   os.println("library(poibin)");
			   os.println("pv<-1-ppoibin(kk=kk, pp=pp, method = \"RF\")");
			   os.println("pv[pv<=0]<-min(pv[pv>0])");
			   os.println("write(pv, ncolumns=1, file=\"" + tmpFile + ".pv" + "\")");
			   os.flush();
			   os.close();
			   
			   MyFunc.runRscript(tmpFile + ".R");
			   List<Double> p = MyFunc.readDoubleList(tmpFile + ".pv");
		   
			   File f = new File(tmpFile + ".kk");
			   if(!f.delete()){
				   f.deleteOnExit();
		   		}
			   f = new File(tmpFile + ".pp");
			   if(!f.delete()){
				   f.deleteOnExit();
			   	}
			   f = new File(tmpFile + ".R");
			   if(!f.delete()){
				   f.deleteOnExit();
			   	}
			   f = new File(tmpFile + ".pv");
			   if(!f.delete()){
				   f.deleteOnExit();
			   }
			   return p;
		   } catch (Exception e) {
			   System.err.println(kk);
			   System.err.print(pp);
			   	throw  new ArithmeticException();
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
	
	protected String generateRandomProbes(){
		int i = (int) Math.floor(Math.random()*probeList.size());
		return probeList.get(i);
	}
	
	protected String generateRandomSamples(){
		int i = (int) Math.floor(Math.random()*sampleList.size());
		return sampleList.get(i);
	}
	protected void setLogRcutoffByPercentile(double LogRcutoffByPercentile){
		List <Double> nullLogR = new ArrayList<Double>();
		for(int i = 0; i < 500; i++){
			String p = generateRandomProbes();
			String s = generateRandomSamples();
			int chr = probeInfo.chr(p);
			int pos = probeInfo.pos(p);
			double logr = LogR.get(s).get(chr, pos).value();
			nullLogR.add(logr);
		}
		if(LogRcutoffByPercentile > 1){
			LogRcutoffByPercentile /= 100;
		}
		ampCutoff = MyFunc.percentile(nullLogR,1-LogRcutoffByPercentile);
		delCutoff = MyFunc.percentile(nullLogR,LogRcutoffByPercentile);
		System.err.println("ampCutoff=" + ampCutoff);
		System.err.println("delCutoff=" + delCutoff);
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
	

	public MyMat getQvalues(){
		return qvalues;
	}
	
	public MyMat getPvalues(){
		return pvalues;
	}
	
	public MyMat getScores(){
		return scores;
	}
	
	public void perform(){
		System.err.println("calculate scores....");
		calculateScores();
		//System.err.println(scores);
		showAberrationRatio();
		
		System.err.println("calculate pvalues....");
		calculatePvalues();
		
		System.err.println("calculate qvalues....");
		calculateQvalues();
	}
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("a", "ampcutoff", true, "cutoff for amplification call");
		options.addOption("d", "delcutoff", true, "cutoff for deletion call");
		options.addOption("A", "aicutoff", true, "cutoff for AI call");
		options.addOption("n", "ndsize", true, "size of null distributions");
		options.addOption("o", "outfile", true, "output file name");
		options.addOption("t", "thinpi", true, "thin down probe info");
		options.addOption("n", "nump", true, "# of psuedo probes");
		options.addOption("c", "cutbyper", true, "amp and del cutoff by percentile");
		options.addOption("s", "score", true, "score file name");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine = null;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] BAFSegFile LogRSegFile probeTsvFile", options);
			System.exit(1);
		}
		List <String> argList = commandLine.getArgList();
		if(!(argList.size() == 3 | argList.size() == 2)){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] BAFSegFile LogRSegFile probeTsvFile", options);
			System.exit(1);
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
			int n = 10000;
			if(commandLine.hasOption("n")){
				n = Integer.valueOf(commandLine.getOptionValue("n"));
			}
			PI = BAF.generatePsuedoProbeInfo(n);
		}
		
		PART PART = new PART(BAF,LogR,PI);
		if(commandLine.hasOption("a")){
			PART.ampCutoff = Double.valueOf(commandLine.getOptionValue("a"));
		}
		if(commandLine.hasOption("d")){
			PART.delCutoff = Double.valueOf(commandLine.getOptionValue("d"));
		}
		if(commandLine.hasOption("A")){
			PART.AICutoff = Double.valueOf(commandLine.getOptionValue("A"));
		}
		
		if(commandLine.hasOption("c")){
			PART.setLogRcutoffByPercentile(Double.valueOf(commandLine.getOptionValue("c")));
		}
		
		PART.perform();
		MyMat tmp =  PART.getMinusLogQvalues();
		
		Writer os;
		if(commandLine.hasOption("o")){
			 os = new BufferedWriter(new FileWriter(commandLine.getOptionValue("o")));
		}else{
			os = new PrintWriter(System.out);
		}
		os.write("Probe" + "\t" +  "Chrom" + "\t" +  "BasePair" + "\t" + MyFunc.join("\t", PART.scoreTypes) + "\n");
		for(String s: PI.getProbeSet()){
			os.write(s + "\t" +  PI.chr(s) + "\t" +  PI.pos(s) + "\t" + MyFunc.join("\t", tmp.getRow(s))+ "\n");			
		}
		os.flush();
		if(commandLine.hasOption("s")){
			 os = new BufferedWriter(new FileWriter(commandLine.getOptionValue("s")));
			 tmp = PART.getScores();
			 os.write("Probe" + "\t" +  "Chrom" + "\t" +  "BasePair" + "\t" + MyFunc.join("\t", PART.scoreTypes)+ "\n");
			 for(String s: PI.getProbeSet()){
				 os.write(s + "\t" +  PI.chr(s) + "\t" +  PI.pos(s) + "\t" + MyFunc.join("\t", tmp.getRow(s))+ "\n");			
			 }
			 os.flush();
		}
	}

}
