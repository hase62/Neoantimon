package clone;

import java.io.IOException;
import java.util.*;
import java.util.zip.DataFormatException;

import mutation.MutualExclusivityTest;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.math3.distribution.BinomialDistribution;

import sun.reflect.Reflection;
import utility.MyFunc;
import utility.MyMat;

public class CCFP {
	Map <String, Integer> nA;
	Map <String, Integer> nB;
	Map <String, Integer> mac;
	Map <String, Integer> depth;
	Map <String, Double> purity;
	List <String> mutID;
	
	Map <String, List<Double>> ccfProbability;
	Map <String, Double> mean;
	Map <String, Double> upperBound;
	Map <String, Double> lowerBound;
	
	boolean mutBeforeCNalteration = false;

	
	public CCFP(String infile) throws IOException, DataFormatException{
		MyMat M = new MyMat(infile);
		mutID = M.getRowNames();
		if(M.containsColName("nA") ){
			nA = new HashMap<String, Integer>();
			for(String s:mutID){
				nA.put(s, (int)M.get(s, "nA"));	
			}
		}else{
			throw new DataFormatException();	
		}
		if(M.containsColName("nB") ){
			nB = new HashMap<String, Integer>();
			for(String s:mutID){
				nB.put(s, (int)M.get(s, "nB"));	
			}
		}else{
			throw new DataFormatException();	
		}
		if(M.containsColName("depth") ){
			depth = new HashMap<String, Integer>();
			for(String s:mutID){
				depth.put(s, (int)M.get(s, "depth"));	
			}
		}else{
			throw new DataFormatException();	
		}
		if(M.containsColName("purity") ){
			purity = new HashMap<String, Double>();
			for(String s:mutID){
				purity.put(s, M.get(s, "purity"));	
			}
		}else{
			throw new DataFormatException();	
		}
		if(M.containsColName("maf")){
			mac = new HashMap<String, Integer>();
			for(String s:mutID){
				mac.put(s, (int) ((M.get(s, "maf") <= 1)? Math.round(M.get(s, "maf")*depth.get(s)) : M.get(s, "maf")));	
			}
		}else{
			throw new DataFormatException();	
		}
		ccfProbability = new HashMap<String,List<Double>>();
		mean = new HashMap<String,Double>();
		upperBound = new HashMap<String,Double>();
		lowerBound = new HashMap<String,Double>();		
	}
	
	
	
	void calculateFccProbability(){
		List <String> mutID2 = new ArrayList<String>();
		for(String s:mutID){
			if(nA.get(s) == 0 &  nB.get(s) == 0){
				continue;
			}
			System.err.println(s);
			List <Double > ccfP =  new ArrayList<Double>();
			for(int i = 0; i< 100; i++){
				ccfP.add(0.0); 
			}
			if(nA.get(s) == nB.get(s)){
				if( nA.get(s) ==1){ // cn neutral
					for(int i = 0; i< 100; i++){
						ccfP.set( i, ccfP.get(i) + LikelihoodFunction(s,0.01 * (i +1), 1)); 
					}
				}else{
					if(!mutBeforeCNalteration){
						for(int i = 0; i< 100; i++){
							ccfP.set( i, ccfP.get(i) + LikelihoodFunction(s,0.01 * (i +1), 1)); 
							ccfP.set( i, ccfP.get(i) + LikelihoodFunction(s,0.01 * (i +1), nA.get(s))); 
						}
					}else{
						for(int i = 0; i< 100; i++){
							ccfP.set( i, ccfP.get(i) + LikelihoodFunction(s,0.01 * (i +1), nA.get(s))); 
						}
					}
				}
			}else if(nA.get(s) * nB.get(s) == 0 ){ //LOH
				if(mutBeforeCNalteration){
					int q = nA.get(s) + nB.get(s);
					for(int i = 0; i< 100; i++){
						ccfP.set( i, ccfP.get(i) + LikelihoodFunction(s,0.01 * (i +1), q));
					}
				}else{
					for(int i = 0; i< 100; i++){
						ccfP.set( i, ccfP.get(i) + LikelihoodFunction(s,0.01 * (i +1), 1)); 
					}
					int q = nA.get(s) + nB.get(s);
					if(q >1){
						for(int i = 0; i< 100; i++){
							ccfP.set( i, ccfP.get(i) + LikelihoodFunction(s,0.01 * (i +1),q)); 
						}
					}
				}
			}else{
				if(!mutBeforeCNalteration){
					for(int i = 0; i< 100; i++){
						ccfP.set( i, ccfP.get(i) + LikelihoodFunction(s,0.01 * (i +1), 1)*2); 
					}
				}
				if(nA.get(s) >1){
					for(int i = 0; i< 100; i++){
						ccfP.set( i, ccfP.get(i) + LikelihoodFunction(s,0.01 * (i +1),nA.get(s))); 
					}
				}
				if(nB.get(s) >1){
					for(int i = 0; i< 100; i++){
						ccfP.set( i, ccfP.get(i) + LikelihoodFunction(s,0.01 * (i +1),nB.get(s))); 
					}
				}
			}
			double tmp = MyFunc.sum(ccfP);
			for(int i = 0; i< 100; i++){
				ccfP.set( i,ccfP.get(i)/tmp);
			}
			ccfProbability.put(s, ccfP);
			mutID2.add(s);
		}
		mutID = mutID2;
	}
	
	private double LikelihoodFunction(String s, double c, double p){ 
		// s= mutID, c = ccf, p = mutated allele copy number
		double q = nA.get(s) + nB.get(s); //q = total copy number 
		double alpha = purity.get(s); // alpha = purity
		// f = mutated allele fraction
		double f = alpha * c * p;  
		f /= 2*(1 - alpha)  +  alpha * q;
		return (new BinomialDistribution(depth.get(s), f)).probability(mac.get(s));
	}
	
	void calculateStatistics(){
		for(String s:mutID){
			double tmp = 0;	
			for(int i = 0; i< 100; i++){
				tmp += ccfProbability.get(s).get(i);
				if(tmp > 0.05 & !lowerBound.containsKey(s) ){
					lowerBound.put(s, (i+1)*0.01);
				}
				if(tmp > 0.5 & !mean.containsKey(s) ){
					mean.put(s, (i+1)*0.01);
				}
				if(tmp > 0.95 & !upperBound.containsKey(s) ){
					upperBound.put(s, (i+1)*0.01);
				}
			}
		}
		
	}
	
	void printResult(){
		System.out.println("\t" + "mean" + "\t" + "lowerBound" + "\t" + "upperBound");
		for(String s:mutID){
			System.out.println(s + "\t" + mean.get(s) + "\t" + lowerBound.get(s) + "\t" + upperBound.get(s));	
		}
	}
	
	MyMat getCcfProbabilityMatrix(){
		List <String> L = new ArrayList <String>();
		for(int i = 0; i< 100; i++){
			String tmp = "ccf" + (i+1)*0.01;
			if(tmp.length() > 7){
				tmp = tmp.substring(0,7);
			}
			L.add(tmp);
		}
		MyMat M = new MyMat(mutID,L);
		for(int i = 0; i < mutID.size(); i++){
			for(int j = 0; j< 100; j++){
				M.set(i,j, ccfProbability.get(mutID.get(i)).get(j));
			}
		}
		return M;
	}
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("m", "mat", true, "get probabirity matrix");
		options.addOption("b", "bf", false, "assume before copy number alteration");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] tabfile (containing 'nA', 'nB', 'maf', 'depth' and 'purity' columns)", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 1){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] tabfile (containing 'nA', 'nB', 'maf', 'depth' and 'purity' columns)", options);
			return;
		}
			
		CCFP CCFP = new CCFP(argList.get(0));
		if(commandLine.hasOption("b")){
			CCFP. mutBeforeCNalteration = true;
		}
		CCFP.calculateFccProbability();
		CCFP.calculateStatistics();
		CCFP.printResult();
		if(commandLine.hasOption("m")){
			CCFP.getCcfProbabilityMatrix().print(commandLine.getOptionValue("m"));
		}
	}

}
