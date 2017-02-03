package snp;

import java.util.*;

import sun.reflect.Reflection;
import utility.*;

import org.apache.commons.cli.*;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.*;

public class AllelicImbalanceTest{
	private MyMat AI; // 0, No AI; 1, AI to A; 2, AI to B
	
	private List <Double> binomialP;
	private List <Double> correctedP;
	private List <Integer> A;
	private List <Integer> B;
	
	private List <Double> BratioForEachArray;
	
	private Double Alabel = 1.0;
	private Double Blabel = 2.0;
	 
	private MyMat result;
	
	private BinomialDistributionImpl BD;
	
	private List <Double> topPermutationP;
	private int numberOfPermutations = 100;
	
	public AllelicImbalanceTest(MyMat AI){
		System.err.println("Setuping the test.....");
		this.AI = AI;	
		BD = new BinomialDistributionImpl(0, 0.5);
	}
	
	public void perform(){
		calculateBinomialP();
		calculateBratioForEachArray();
		perfomPermutation();
		calculateCorrectedP();
		prepareResultMatrix();
	}
	
	public MyMat getResult(){
		return result;
	}
	
	private void calculateBinomialP(){
		System.err.println("Caluluating Binomial Pvalues.....");
		binomialP  = new ArrayList<Double>();
		A = new ArrayList<Integer>();
		B = new ArrayList<Integer>();
		int j,i,ncol,nrow;
		nrow = AI.rowSize();
		ncol = AI.colSize();
		for(i=0; i < nrow; i++){
			int a = 0;
			int b = 0;
			for(j=0; j < ncol; j++ ){
				if(AI.get(i, j)==Alabel){
					a++;
				}else if(AI.get(i,j)==Blabel){
					b++;
				}
			}
			A.add(a);
			B.add(b);
			double p=1;
			if(a+b>0){
				int m = Math.min(a,b);
				BD.setNumberOfTrials(a+b);
				try {
					p = 2*BD.cumulativeProbability(m);
				} catch (MathException e) {
					e.printStackTrace();
				}
				if(p>1){
					p=1;
				}
			}
			binomialP.add(p);
		}
	}

	private void calculateBratioForEachArray(){
		System.err.println("Caluluating B ratio for each aray.....");
		BratioForEachArray = new ArrayList<Double>();
		int j,i,ncol,nrow;
		nrow = AI.rowSize();
		ncol = AI.colSize();
		for(j=0; j < ncol; j++ ){
			double a = 0;
			double b = 0;
			for(i=0; i < nrow; i++){
				if(AI.get(i, j)==Alabel){
					a++;
				}else if(AI.get(i,j)==Blabel){
					b++;
				}
			}
			if(a+b==0){
				System.err.println("AI err: " +  j + "-th col has no value!");
				BratioForEachArray.add(Double.NaN);
			}else{
				BratioForEachArray.add(b/(a+b));
			}
		}
	}
	
	
	private void perfomPermutation(){
		System.err.println("Perform Permutaion.....");
		int i;
		topPermutationP = new ArrayList<Double>();
		for(i = 0; i < numberOfPermutations; i++){
			System.err.println(i + "/" + numberOfPermutations);
			MyMat permAI = permutateAI();
			List <Double> P = calculateBinomialP(permAI);
			topPermutationP.add(MyFunc.min(P));
		}
	}
	
	private void calculateCorrectedP(){
		System.err.println("Calculate corrected Pvalues.....");
		int i,j;
		int n = binomialP.size();
		int m = topPermutationP.size();
		correctedP = new ArrayList<Double>();
		for(i = 0; i < n; i++){
			double k = 0;
			for(j = 0; j < m; j++){
				if(binomialP.get(i) > topPermutationP.get(j)){
					k++;
				}
			}
			double p= k/m;
			correctedP.add(p);	
		}
	}
	
	private void prepareResultMatrix(){
		System.err.println("Prepare a result matrix.....");
		List <String> field = new ArrayList<String>();
		field.add("A");
		field.add("B");
		field.add("binomialP");
		field.add("correctedP");
		List <String> id = AI.getRowNames();
		result = new MyMat(id,field);
		int i;
		for(i=0;i<result.rowSize();i++){
			result.set(i, 0, A.get(i));
			result.set(i, 1, B.get(i));
			result.set(i, 2, binomialP.get(i));
			result.set(i, 3, correctedP.get(i));
		}
		Map <String,Double> tmp = result.getColMap(2);
		id = MyFunc.sortKeysByAscendingOrderOfValues(tmp);
		result.reorderRows(id);
	}
	
	
	private List <Double>  calculateBinomialP(MyMat AI){
		List <Double> P = new ArrayList<Double>();
		int j,i,ncol,nrow;
		nrow = AI.rowSize();
		ncol = AI.colSize();
		for(i=0; i < nrow; i++){
			int a = 0;
			int b = 0;
			for(j=0; j < ncol; j++ ){
				if(AI.get(i, j)==Alabel){
					a++;
				}else if(AI.get(i,j)==Blabel){
					b++;
				}
			}
			double p=1;
			if(a+b>0){
				int m = Math.min(a,b);
				BD.setNumberOfTrials(a+b);
				try {
					p = 2*BD.cumulativeProbability(m);
				} catch (MathException e) {
					e.printStackTrace();
				}
				if(p>1){
					p=1;
				}
			}
			P.add(p);
		}
		return P;
	}
	
	private MyMat permutateAI(){
		MyMat permAI = new MyMat(AI.getRowNames(), AI.getColNames());
		int j,i,ncol,nrow;
		nrow = AI.rowSize();
		ncol = AI.colSize();
		Random rnd = new Random();
		for(i=0; i < nrow; i++){
			for(j=0; j < ncol; j++ ){
				if(AI.get(i, j)==Alabel||AI.get(i, j)==Blabel){
					if(rnd.nextDouble() > BratioForEachArray.get(j)){
						permAI.set(i,j,Alabel);
					}else{
						permAI.set(i,j,Blabel);
					}
				}else{
					permAI.set(i,j,Double.NaN);
				}
			}
		}
		return permAI;
	}
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("a", "alabel", true, "value of a label");
		options.addOption("b", "blabel", true, "value of b label");
		options.addOption("p", "perm", true, "number of permutaions");
		options.addOption("o", "outfile", true, "output file name");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] AIfile", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 1){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] AIfile", options);
			return;
		}
		MyMat AI = new MyMat(argList.get(0));
		AllelicImbalanceTest AIT = new AllelicImbalanceTest(AI);
		if(commandLine.hasOption("a")){
			AIT.Alabel = Double.valueOf(commandLine.getOptionValue("a"));
		}
		if(commandLine.hasOption("b")){
			AIT.Blabel = Double.valueOf(commandLine.getOptionValue("b"));
		}
		if(commandLine.hasOption("p")){
			AIT.numberOfPermutations = Integer.valueOf(commandLine.getOptionValue("p"));
		}
		AIT.perform();
		if(commandLine.hasOption("o")){
			AIT.getResult().print(commandLine.getOptionValue("o"));
		}else{
			AIT.getResult().print();
		}
	}
}