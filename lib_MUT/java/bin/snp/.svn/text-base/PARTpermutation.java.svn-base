package snp;

import java.util.*;
import java.io.*;
import org.apache.commons.cli.*;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.*;

import sun.reflect.Reflection;
import utility.MyFunc;
import utility.MyMat;


public class PARTpermutation extends PARTsimple {

	
	
	protected int nullDistSize = 100000;
	protected List <Integer> nullDist;
	
	
	public PARTpermutation(SegmentContainerMap baf,  ProbeInfo pi) throws IOException{
		super(baf,pi);
		
	}
	
	protected String generateRandomProbes(){
		int i = (int) Math.floor(Math.random()*probeList.size());
		return probeList.get(i);
	}
	protected String generateRandomSamples(){
		int i = (int) Math.floor(Math.random()*sampleList.size());
		return sampleList.get(i);
	}
	
	protected List<String> shufflesProbes(){
		List <String> tmp = probeList;
		Collections.shuffle(tmp);
		return tmp;
	}
	
	
	/*protected void generateNullDist(){
		nullDist = new ArrayList<Integer>();
		for(int i = 0; i <  probeList.size(); i++){
			nullDist.add(0);
		}
		Map<String, Double> poissonBinomialPnull = new HashMap<String, Double>();
		for(String s: sampleSet()){
			poissonBinomialPnull.put(s,0.0);
		}
		for(String s: sampleSet()){
			SegmentContainer bafSC = BAF.get(s);
			int i = 0;
			for(String p: shufflesProbes()){
				int chr = probeInfo.chr(p);
				int pos = probeInfo.pos(p);
				double baf = bafSC.get(chr, pos).value();
				if((countLess & baf<= AICutoff)| (!countLess & baf>= AICutoff)){
					nullDist.set(i,nullDist.get(i)+1);
					poissonBinomialPnull.put(s, poissonBinomialPnull.get(s) + 1);
				}
				i++;
			}
		}
		
		for(String s: sampleSet()){
			poissonBinomialPnull.put(s, poissonBinomialPnull.get(s)/probeSet().size());
		}
		for(String s: sampleSet()){
		System.err.println(s + " " + poissonBinomialP.get(s) + " " + poissonBinomialPnull.get(s));
		}
	}*/
	
	protected void generateNullDist() {
		nullDist = new ArrayList<Integer>();
		for(int i = 0; i <  nullDistSize; i++){
			nullDist.add(0);
		}
		for(String s: sampleSet()){
			SegmentContainer bafSC = BAF.get(s);
		
			for(int i = 0; i <  nullDistSize; i++){
				String randProbe = generateRandomProbes();
				int chr = probeInfo.chr(randProbe);
				int pos = probeInfo.pos(randProbe);
				
				double baf = 0;
				//try{
					baf = bafSC.get(chr, pos).value();
				//}catch(Exception e){
					//continue;
					//System.err.println(chr + " " + pos);
				//}
				if((countLess & baf<= AICutoff)| (!countLess & baf>= AICutoff)){
					nullDist.set(i,nullDist.get(i)+1);
				}
				
				//if(Math.random() <= poissonBinomialP.get(s)){
				//	nullDist.set(i,nullDist.get(i)+1);
				//}
			}
		}
	}
	
	protected  void calculatePvalues(){
		generateNullDist();
		pvalues =  new LinkedHashMap<String,Double>();;
		for(String p: probeList){
				int score = (int)(double)scores.get(p);
				double pvalue = 0;
				for(int i = 0; i < nullDist.size(); i++){
					if(score<nullDist.get(i)){
						pvalue++;
					}
				}
				pvalue = (pvalue==0)?1:pvalue;
				pvalue/= nullDist.size();
				pvalues.put(p, pvalue);
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
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("c", "cutoff", true, "cutoff for AI call");
		options.addOption("o", "outfile", true, "output file name");
		options.addOption("t", "thinpi", true, "thin down probe info");
		options.addOption("n", "nump", true, "# of psuedo probes");
		options.addOption("O", "optcutoff", true, "optimize cutoff with [lower:upper:numberOfCutoffs]");
		options.addOption("l", "less", false, "count less than cutoff");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine = null;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] SegFile  probeTsvFile", options);
			System.exit(1);
		}
		List <String> argList = commandLine.getArgList();
		if(!(argList.size() == 1 | argList.size() == 2)){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] SegFile probeTsvFile", options);
			System.exit(1);
		}
		SegmentContainerMap BAF = new SegmentContainerMap(argList.get(0));
		ProbeInfo PI;
		List<String> sample = new ArrayList<String>(BAF.keySet());
		if(argList.size() == 2 ){
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
		
		PARTpermutation PART = new PARTpermutation(BAF,PI);
		if(commandLine.hasOption("c")){
			PART.AICutoff = Double.valueOf(commandLine.getOptionValue("c"));
		}
		if(commandLine.hasOption("l")){
			PART.countLess = true;
		}
		if(commandLine.hasOption("O")){
			List <String>tmp = MyFunc.split(":",commandLine.getOptionValue("O"));
			if(tmp.size()==3){
				PART.AIUpperCutoff = Double.valueOf(tmp.get(1));
				PART.AILowerCutoff = Double.valueOf(tmp.get(0));
				PART.numberOfTriedCutoffs = Integer.valueOf(tmp.get(2));
				PART.optimizeCutoffs();
			}else{
				formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] SegFile  probeTsvFile", options);
				System.exit(1);
			}
		}
		
		PART.perform();
		
		Writer os;
		if(commandLine.hasOption("o")){
			 os = new BufferedWriter(new FileWriter(commandLine.getOptionValue("o")));	 
		}else{
			os = new PrintWriter(System.out);
		}
			os.write("Probe" + "\t" +  "Chrom" + "\t" +  "BasePair" + "\t" + "score" + "\t" + "pvalue" + "\t" + "qvalue" + "\n");
		for(String s: PI.getProbeSet()){
			os.write(s + "\t" +  PI.chr(s) + "\t" +  PI.pos(s) + "\t" + PART.scores.get(s) + "\t" + PART.getMinusLogPvalues().get(s) + "\t" + PART.getMinusLogQvalues().get(s) + "\n");			
		}
		os.flush();
	}
		

	
	
}
