package mutation;

import java.io.IOException;
import java.util.*;
import java.util.zip.DataFormatException;

import network.Link;
import network.NullLinkGenerator;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.math3.distribution.BinomialDistribution;

import sun.reflect.Reflection;
import utility.MyFunc;




public class HotNeighborTest {
	private	Link link;
	private List <String> genes;
	private	List <String> allGenes;
	private Map <String,Integer> mutationCount;
	private Map <String,Integer> mutationSumInNeighbor;
	private Map <String,Integer> neighborCount;
	private Map <String, List<String>> hotNeighbor;
	private Map <String,Double> pvalue;
	private int totalMutationCount;
	private double pcutoff = 3;
	private double pmax = 20;
	private int itrq = 20;
	private Map <String,Double> qvalue;
	
	public HotNeighborTest() {}
	
	public void setLink(Link link){
		this.link = link;
		allGenes = link.getNodeName();
	}
	public void setLink(String infile) throws IOException, DataFormatException{
		link = new Link(infile); 
		allGenes = link.getNodeName();
	}
	
	public void setPcutoff(double d){
		pcutoff = d;
	}
	
	public void setItrQ(int d){
		itrq = d;
	}
	
	public void applyCeiling2MutationCount(int i){
		for(String s: mutationCount.keySet()){
			if(mutationCount.get(s) > i){
				mutationCount.put(s,i);
			}
		}
		setTotalMutationCount();
	}
	
	public void cutSmallMutationCount(int i){
		Set <String> tmp = new HashSet<String>(mutationCount.keySet());
		for(String s: tmp){
			if(mutationCount.get(s) < i){
				mutationCount.remove(s);
			}
		}
		setTotalMutationCount();
		genes = new ArrayList<String>(mutationCount.keySet());
	}
	
	private void setTotalMutationCount(){
		mutationCount = new HashMap<String, Integer>();
		totalMutationCount = 0;
		for(String k: mutationCount.keySet()){
			totalMutationCount += Integer.valueOf(mutationCount.get(k));
		}
	}
	
	
	public void setMutationCount(String infile) throws IOException, DataFormatException{
		mutationCount = new HashMap<String, Integer>();
		Map<String, String> tmp = MyFunc.readStringStringMap(infile);
		for(String k: tmp.keySet()){
			if(link.containsNode(k)){
				mutationCount.put(k, Integer.valueOf(tmp.get(k)));
			}
		}
		setTotalMutationCount();
		genes = new ArrayList<String>(mutationCount.keySet());
	
	}
	
	public void setMutationCount(Map<String, Integer> count){
		mutationCount = new HashMap<String, Integer>();
		for(String k: count.keySet()){
			if(link.containsNode(k)){
				mutationCount.put(k, count.get(k));
			}
		}
		setTotalMutationCount();
		genes = new ArrayList<String>(mutationCount.keySet());
	}
	
	private void getMutationSumInNeighbor(){
		mutationSumInNeighbor = new HashMap<String, Integer>();
		neighborCount = new HashMap<String, Integer>();
		hotNeighbor = new HashMap<String, List<String>>();
		for(String gene: genes){
			List <String> neighbors = link.getNeighbors(gene);
			List <String> hot = new ArrayList<String>();
			int i = 0;
			int j = 0;
 			for(String neighbor: neighbors){
				if(mutationCount.containsKey(neighbor)){
					hot.add(neighbor);
					i += mutationCount.get(neighbor);
				}
				j++;
			}
 			hotNeighbor.put(gene, hot);
 			mutationSumInNeighbor.put(gene, i);
 			neighborCount.put(gene,j);
		}
	}

	private void calculatePvalues(){
		pvalue = new HashMap<String, Double>();
		for(String gene: genes){
			pvalue.put(gene, calculatePvalue(gene));
		}
	}
	
	private  double calculatePvalue(String gene){
		int n =  totalMutationCount;
		int k = mutationSumInNeighbor.get(gene);
		double pvalue = 1;
		if(n>0 & k > 0){
			double p = (double)(neighborCount.get(gene))/allGenes.size();
			BinomialDistribution BN = new BinomialDistribution(n, p); 
			try {
				pvalue = 1- BN.cumulativeProbability(k); 
			}catch (Exception e){
				e.printStackTrace(); 
			}
		}
		pvalue = -Math.log10(pvalue);
		if(pvalue > pmax){
			pvalue = pmax;
		}
		return pvalue;
	}
	
	public void printResult(){
		System.out.println("gene\tpvalue" + (qvalue==null?"":"\tqvalue") + "\tmutationSumInNeighbor\tneighborCount\tmutatedNeighbor");
		for(String gene: MyFunc.sortKeysByDescendingOrderOfValues(pvalue)){
			if(pvalue.get(gene) >= pcutoff){
				System.out.print(gene + "[" +  mutationCount.get(gene) + "]" + "\t"  + pvalue.get(gene) + (qvalue==null?"":("\t"+qvalue.get(gene))) + "\t"+  mutationSumInNeighbor.get(gene) + "\t" + neighborCount.get(gene)  + "\t");
				List <String> tmp = hotNeighbor.get(gene);
				java.util.Collections.sort(tmp, new geneComp());
				List <String> tmp2 = new ArrayList<String>(); 
				for(String hn: tmp){
					tmp2.add(hn + "[" +  mutationCount.get(hn) + "]");
				}
				System.out.println(tmp2.isEmpty()?"":MyFunc.join(" ", tmp2));
			}
		}
	}
	
	public void printResultSimple(){
		for(String gene: MyFunc.sortKeysByDescendingOrderOfValues(pvalue)){
			if(pvalue.get(gene) >= pcutoff){
				System.out.println(gene  + "\t"  + pvalue.get(gene) + (qvalue==null?"":("\t"+qvalue.get(gene)))); 
			}
		}
	}
	
	
	
	
	public class geneComp implements Comparator<String> {
        public int compare(String g1, String g2) {
        	if(mutationCount.get(g1) < mutationCount.get(g2)){
        		return 1;
        	}else if (mutationCount.get(g1) > mutationCount.get(g2)){
        		return -1;
        	}else{
				return 0;
        	}
        }
	}
	
	public void perform(){
		getMutationSumInNeighbor();
		calculatePvalues();
		if(itrq > 0){
			calculateQvalues();
		}
	}
	
	private List <Double> getNullPvalues(){
		List <Double> nullp = new ArrayList<Double>();
		NullLinkGenerator NLG = new NullLinkGenerator(link);
		for(int i = 0; i < itrq; i++){
			System.err.println((i+1) + "-th iteration....."); 
			Link nullLink = NLG.getRondomNetwork();
			HotNeighborTest HNA = new  HotNeighborTest();
			HNA.setLink(nullLink);
			HNA.setMutationCount(mutationCount);
			HNA.getMutationSumInNeighbor();
			HNA.calculatePvalues();
			nullp.addAll(HNA.pvalue.values());
		}
		return nullp;
	}
	
	private void calculateQvalues(){
		qvalue = new HashMap <String, Double>();
		List <Double> tmp =  getNullPvalues();
		for(String s: pvalue.keySet()){
			double p = pvalue.get(s);
			double q = 0;
			for(double d: tmp){
				if(d > p){
					q++;
				}
			}
			q /= tmp.size();
			qvalue.put(s,q);
		}	
	}
	
	public static void main(String [] args) throws Exception{	
		Options options = new Options();
		options.addOption("p", "pcut", true,  "pvalue cutoff");
		options.addOption("s", "simple", false,  "print simple result");
		options.addOption("i", "itr", true,  "iteration number for qvalue caluculation");
		options.addOption("c", "ceiling", true,  "mutation count ceiling");
		options.addOption("C", "cut", true,  "mutation count cutoff ");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options]  linkFile geneCountFile", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 2){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] linkFile geneCountFile", options);
			return;
		}
		HotNeighborTest H = new HotNeighborTest();
		H.setLink(argList.get(0));
		H.setMutationCount(argList.get(1));
		if(commandLine.hasOption("p")){
			H.setPcutoff(Integer.valueOf(commandLine.getOptionValue("p")));
		}
		if(commandLine.hasOption("i")){
			H.setItrQ(Integer.valueOf(commandLine.getOptionValue("i")));
		}
		if(commandLine.hasOption("c")){
			H.applyCeiling2MutationCount(Integer.valueOf(commandLine.getOptionValue("c")));
		}
		if(commandLine.hasOption("C")){
			H.cutSmallMutationCount(Integer.valueOf(commandLine.getOptionValue("C")));
		}
		H.perform();
		if(commandLine.hasOption("s")){
			H.printResultSimple();
		}else{
			H.printResult();
		}
	}
}
	
	
	
	
	
	
	
	
	
