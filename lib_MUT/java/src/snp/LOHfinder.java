package snp;

import java.io.*;
import java.util.*;
import java.util.zip.DataFormatException;

import utility.*;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.*;
import org.apache.commons.cli.*;
import sun.reflect.Reflection;


public class LOHfinder {
	
	protected SegmentContainerMap  BAF; 
	protected SegmentContainerMap  LogR;
	protected SegmentContainerMap  LOH;
	
	protected List<String>sampleList;
	protected Set<String>sampleSet; 
	
	protected double AICutoff =  0.65;
	protected double ampCutoff = 0.1;
	protected double delCutoff = -0.1;
	
	protected List<Condition> ConditionFunctions; 
	
	interface Condition {
		boolean check(double logr, double baf);
	}
	
	private void initializeConditionFunctions(){
		ConditionFunctions = new ArrayList<Condition>();

		//LOH
		ConditionFunctions.add(
				new Condition(){
					public boolean check(double logr, double baf){
						return baf>= AICutoff & logr <= ampCutoff;
					}
				}
			);
		//LOHbD
		ConditionFunctions.add(
				new Condition(){
					public boolean check(double logr, double baf){
						return baf>= AICutoff & logr < delCutoff;
					}
				}
			);
		//UPD
		ConditionFunctions.add(
				new Condition(){
					public boolean check(double logr, double baf){
						return baf>= AICutoff & logr >= delCutoff & logr <= ampCutoff;
					}
				}
			);
		
	}
	
	int conditionIndex = 0;
	
	public LOHfinder(SegmentContainerMap baf,SegmentContainerMap logr)throws IOException, DataFormatException {
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
		
		sampleSet = new HashSet<String>();
		for(String s: LogR.keySet()){
			if(BAF.containsKey(s) & LogR.containsKey(s)){
				sampleSet.add(s);
			}
		}
		sampleList = new ArrayList<String>(sampleSet);
		
		LOH= new SegmentContainerMap();
		initializeConditionFunctions();
	}
	
	public void find(){
		for(String s: sampleSet){
			SegmentContainer bafSC = BAF.get(s);
			SegmentContainer logrSC = LogR.get(s);
			
			Set<Integer> tmp = new LinkedHashSet<Integer>(bafSC.chrList());
			Set<Integer> tmp2 = new LinkedHashSet<Integer>(logrSC.chrList());
			tmp.retainAll(tmp2);
			List <Integer> chrList = new ArrayList<Integer>(tmp);
			
			SegmentContainer lohSC = new SegmentContainer();
			
			for(Integer I: chrList){
				List <Segment> bafSL = bafSC.getByChr(I);
				List <Segment> logrSL = logrSC.getByChr(I);
				
				int i=0,j=0;	
				while(i < bafSL.size() & j < logrSL.size()){
					Segment b = bafSL.get(i);
					Segment l = logrSL.get(j);
					
					if(b.end() < l.start()){
						i++;
						continue;
					}else if (l.end() < b.start()){
						j++;
						continue;
					}
					
					int start = (b.start()>l.start())?b.start():l.start(); 
					int end = (b.end()<l.end())?b.end():l.end();
					int value = ConditionFunctions.get(conditionIndex).check(l.value(), b.value())?1:0;
					int count = (end -start)/10000;
					Segment newSeg = new Segment(I, start, end, count, value);
					lohSC.add(newSeg);
					if(b.end()<l.end()){
						i++;
					}else{
						j++;
					}	
				}	
			}	
			LOH.put(s,lohSC);
		}
		LOH.mergeSegments();
	}
	
	public SegmentContainerMap getLOH(){
		return LOH;
		
	}
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("a", "ampcut", true, "cutoff for amplification call");
		options.addOption("d", "delcut", true, "cutoff for deletion call");
		options.addOption("A", "aicut", true, "cutoff for AI call");
		options.addOption("u", "upd", false, "call upd");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine = null;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] BAFSegFile LogRSegFile", options);
			System.exit(1);
		}
		List <String> argList = commandLine.getArgList();
		if(!(argList.size() == 2)){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] BAFSegFile LogRSegFile", options);
			System.exit(1);
		}
		SegmentContainerMap BAF = new SegmentContainerMap(argList.get(0));
		SegmentContainerMap LogR = new SegmentContainerMap(argList.get(1));
		
		LOHfinder LF = new LOHfinder(BAF, LogR);
		
		if(commandLine.hasOption("c")){
			LF.ampCutoff = Double.valueOf(commandLine.getOptionValue("a"));
		}
		if(commandLine.hasOption("d")){
			LF.delCutoff = Double.valueOf(commandLine.getOptionValue("d"));
		}
		if(commandLine.hasOption("A")){
			LF.AICutoff = Double.valueOf(commandLine.getOptionValue("A"));
		}
		if(commandLine.hasOption("u")){
			LF.conditionIndex = 1;
		}
		
		LF.find();
		LF.getLOH().print();
	}

}
