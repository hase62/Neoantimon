package mutation;

import java.io.*;
import java.util.*;
import java.util.zip.DataFormatException;

import org.apache.commons.cli.*;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.*;

import sun.reflect.Reflection;
import utility.*;

public class NSmutationEnrichmentTest {
	
	private Set <String> genes;
	private Map <String, Integer> total; 
	private Map <String, Integer> NS;
	
	private Map <String, List<String>> geneset;
	private Set <String> targetGenes;	
	
	private Double NSratio;	
	
	private Map <String, Double> pvalues;
	private Map <String, Double> qvalues;
	private Map <String, Integer> total2;
	private Map <String, Integer> NS2;
	
	private Map <String, List<String>> geneset2;
	
	private int n = 0; // target genes with more than n recurrent mutations
	private int m = 10000; //  calculate NSratio using  genes with less than m recurrent mutations
	
	public NSmutationEnrichmentTest(){};
	
	public void readCountFile(String infile) throws IOException{
		BufferedReader inputStream = new BufferedReader(new FileReader(infile));
		String line;
		genes = new LinkedHashSet<String>();
		total = new HashMap<String, Integer>();
		NS = new HashMap<String, Integer>();
		while((line = inputStream.readLine()) != null){
			List<String> tmp = Arrays.asList((line.split("\t")));
			genes.add(tmp.get(0));
			total.put(tmp.get(0), Integer.valueOf(tmp.get(1)));
			NS.put(tmp.get(0), Integer.valueOf(tmp.get(2)));	
		}
	}
	
	public void readGenesetFile(String infile) throws IOException, DataFormatException{
		geneset = MyFunc.readGeneSetFromGmtFile(infile);
	}
	
	public void readGeneSetList(String infile) throws IOException{
		geneset = new HashMap<String, List<String>>();
		geneset.put(infile, MyFunc.readStringList2(infile));
	}
	
	public void readTargetGenetList(String infile) throws IOException{
		targetGenes = new HashSet <String> (MyFunc.readStringList2(infile));
	}
	
	public void calculateNSratio(){
		double t = 0;
		double ns = 0;
		for(String s: genes){
			if(total.get(s) < m){
				t += total.get(s);		
				ns += NS.get(s);
			}
		}
		NSratio = ns/t;
		System.err.println("NSratio = " + NSratio);
	}
	
	public void calculatePvalues(){
		pvalues = new HashMap<String, Double>();
		total2 = new HashMap<String, Integer>();
		NS2 = new HashMap<String, Integer>();
		geneset2 = new HashMap<String, List<String>>();
		for(String gsid: geneset.keySet()){
			List <String> gs =  geneset.get(gsid);
			int t = 0;
			int ns = 0;
			List<String> tmp = new ArrayList<String>();
			for(String g: gs){
				if(targetGenes != null){
					if(!targetGenes.contains(g)){
						continue;
					}
				}
				if(!genes.contains(g)){
					continue;
				}
				if(!(total.get(g) < m)){
					continue;
				}
				if(total.get(g) > n){
					t += total.get(g);
					ns += NS.get(g);
					tmp.add(g);
				}
			}
			BinomialDistribution BD = new BinomialDistributionImpl(t, NSratio);
			double p = 1;
			try {
				p = 1 - BD.cumulativeProbability(ns-1);
			} catch (MathException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} 
			pvalues.put(gsid, p);
			NS2.put(gsid, ns);
			total2.put(gsid, t);
			geneset2.put(gsid, tmp);
		}
		if(pvalues.size()>10){
				qvalues= MyFunc.calculateStoreyQvalue(pvalues);
		}
	}
	
	public void print(){
		for(String s: MyFunc.sortKeysByAscendingOrderOfValues(pvalues)){
			if(qvalues == null){
				System.out.println(s + "\t" + pvalues.get(s) + "\t" + total2.get(s) + "\t" + NS2.get(s) + "\t" + geneset2.get(s));
			}else{
				System.out.println(s + "\t" + pvalues.get(s) + "\t" +  qvalues.get(s) + "\t" + total2.get(s) + "\t" + NS2.get(s) + "\t" + geneset2.get(s));
			}
		}
	}
	
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("L", "list", false, "gene set from list file");
		options.addOption("r", "rmrec", true, "remove recurrent mutations");
		options.addOption("t", "target", true, "read a target file");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] NSmutCountFile gmtFile", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 2){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] tNSmutCountFile gmtFile", options);
			return;
		}
		
		NSmutationEnrichmentTest NSMET = new  NSmutationEnrichmentTest();
		NSMET.readCountFile(argList.get(0));
		if(commandLine.hasOption("L")){
			NSMET.readGeneSetList(argList.get(1));	
		}else{
			NSMET.readGenesetFile(argList.get(1));	
		}
		if(commandLine.hasOption("r")){
			NSMET.m = Integer.valueOf(commandLine.getOptionValue("r"));
		}
		if(commandLine.hasOption("t")){
			NSMET.readTargetGenetList(commandLine.getOptionValue("t"));
		}
		NSMET.calculateNSratio();
		NSMET.calculatePvalues();
		NSMET.print();
	}
		
}
