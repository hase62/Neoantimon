package snp.old;

import java.io.*;
import java.util.*;

import org.apache.commons.cli.*;

import snp.ProbeInfo;
import snp.Segment;
import snp.SegmentContainer;
import snp.SegmentContainerMap;
import snp.SegmentContainer.SegmentIterator;
import sun.reflect.Reflection;
import utility.MyFunc;

public class SAIRICsimplePaired extends SAIRICsimple {

	protected  SegmentContainerMap BAF0;
	protected  Set <String> sampleSet;
	public SAIRICsimplePaired(){};
	public SAIRICsimplePaired(SegmentContainerMap baf, SegmentContainerMap baf0, ProbeInfo pi) throws IOException {
		BAF = baf;
		probeInfo = pi; 
		probeList = new ArrayList<String>(probeSet());
		BAF0 = baf0;
		sampleSet = new HashSet<String>();
		for(String s:BAF.keySet()){
			if(BAF0.containsKey(s)){
				sampleSet.add(s);
			}
		}
	}
	
	protected Set<String> sampleSet(){
		return sampleSet;
	}
	
	protected void calculateScores(){
		scores =  new LinkedHashMap<String, Double>();
		for(String p: probeSet()){
			scores.put(p,0.0);
		}
		poissonBinomialP = new HashMap<String, Double>();
		for(String s: sampleSet()){
			poissonBinomialP.put(s,0.0);
		}
		for(String s: sampleSet()){
			SegmentContainer bafSC = BAF.get(s);
			SegmentContainer.SegmentIterator bafSCitr  = bafSC.iterator();
			SegmentContainer bafSC0 = BAF0.get(s);
			SegmentContainer.SegmentIterator bafSC0itr  = bafSC0.iterator();
			Segment S = bafSCitr.next();
			Segment S0 = bafSC0itr.next();
			for(String p: probeSet()){
				int chr = probeInfo.chr(p);
				int pos = probeInfo.pos(p);
				
				while(chr > S.chr()){
					S = bafSCitr.next();
				}
				while(pos > S.end()){
					S = bafSCitr.next();
				}
				while(chr > S0.chr()){
					S0 = bafSCitr.next();
				}
				while(pos > S0.end()){
					S0 = bafSCitr.next();
				}
				
				double baf = S.value() - S0.value();
				
				if((countLess & baf<= AICutoff)| (!countLess & baf>= AICutoff)){
					scores.put(p, scores.get(p)+1);
					poissonBinomialP.put(s, poissonBinomialP.get(s) + 1);
				}
			}
		}
		
		for(String s: sampleSet()){
			poissonBinomialP.put(s, poissonBinomialP.get(s)/probeSet().size());
		}
		System.err.println("Ratio of velues above cutoff: " + MyFunc.mean(new ArrayList<Double>(poissonBinomialP.values())));
		
	}
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("c", "cutoff", true, "cutoff for AI call");
		options.addOption("o", "outfile", true, "output file name");
		options.addOption("t", "thinpi", true, "thin down probe info");
		options.addOption("i", "interval", true, "interval for psuedo probe info");
		options.addOption("O", "optcutoff", true, "optimize cutoff with [lower:upper:numberOfCutoffs]");
		options.addOption("l", "less", false, "count less than cutoff");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine = null;
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
			int interval = 10000;
			if(commandLine.hasOption("i")){
				interval = Integer.valueOf(commandLine.getOptionValue("i"));
			}
			PI = ProbeInfo.generatePsuedoProbeInfo(interval);
			PI.filter(BAF.get(sample.get(0)));
			}
		
		SAIRICsimplePaired SAIRIC = new SAIRICsimplePaired(BAF,BAF2,PI);
		if(commandLine.hasOption("c")){
			SAIRIC.AICutoff = Double.valueOf(commandLine.getOptionValue("c"));
		}
		if(commandLine.hasOption("l")){
			SAIRIC.countLess = true;
		}
		if(commandLine.hasOption("O")){
			List <String>tmp = MyFunc.split(":",commandLine.getOptionValue("O"));
			if(tmp.size()==3){
				SAIRIC.AIUpperCutoff = Double.valueOf(tmp.get(1));
				SAIRIC.AILowerCutoff = Double.valueOf(tmp.get(0));
				SAIRIC.numberOfTriedCutoffs = Integer.valueOf(tmp.get(2));
				SAIRIC.optimizeCutoffs();
			}else{
				formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] SegFile  probeTsvFile", options);
				System.exit(1);
			}
		}
		
		SAIRIC.perform();
		
		Writer os;
		if(commandLine.hasOption("o")){
			 os = new BufferedWriter(new FileWriter(commandLine.getOptionValue("o")));	 
		}else{
			os = new PrintWriter(System.out);
		}
			os.write("Probe" + "\t" +  "Chrom" + "\t" +  "BasePair" + "\t" + "score" + "\t" + "pvalue" + "\t" + "qvalue" + "\n");
		for(String s: PI.getProbeSet()){
			os.write(s + "\t" +  PI.chr(s) + "\t" +  PI.pos(s) + "\t" + SAIRIC.scores.get(s) + "\t" + SAIRIC.getMinusLogPvalues().get(s) + "\t" + SAIRIC.getMinusLogQvalues().get(s) + "\n");			
		}
		os.flush();
	}

}
