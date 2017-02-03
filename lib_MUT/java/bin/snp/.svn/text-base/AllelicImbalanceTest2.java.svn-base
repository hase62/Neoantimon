package snp;

import java.util.*;

import org.apache.commons.math.stat.inference.TTestImpl;
import org.apache.commons.cli.*;
import sun.reflect.Reflection;
import utility.*;

public class AllelicImbalanceTest2 {
	private MyMat CN; 
	private MyMat GT; // AA,0; AB,1; BB,2; NA,-1
	
	private List <Double> ttestP;
	private List <Double> correctedP;
	private List <Integer> Acount;
	private List <Double> Amean;
	private List <Double> Asd;
	private List <Integer> Bcount;
	private List <Double> Bmean;
	private List <Double> Bsd;
	
	private Double Alabel = 0.0;
	private Double Blabel = 2.0;
	
	private MyMat result;
	
	private List <Double> topPermutationP;
	private int numberOfPermutations = 100;
	
	private TTestImpl TT;
	
	public AllelicImbalanceTest2(MyMat CN, MyMat GT){
		System.err.println("Setuping the test.....");
		this.CN = CN;
		this.GT = GT;
		TT = new TTestImpl();
	}
	
	public void perform(){
		calculateTtestP();
		perfomPermutation();
		calculateCorrectedP();
		prepareResultMatrix();
	}
	
	public MyMat getResult(){
		return result;
	}
	
	private void calculateTtestP(){
		System.err.println("Caluluating t-test Pvalues.....");
		ttestP  = new ArrayList<Double>();
		Acount = new ArrayList<Integer>();
		Bcount = new ArrayList<Integer>();
		Amean = new ArrayList<Double>();
		Bmean = new ArrayList<Double>();
		Asd = new ArrayList<Double>();
		Bsd = new ArrayList<Double>();
		int j,i,ncol,nrow;
		nrow = CN.rowSize();
		ncol = CN.colSize();
		for(i=0; i < nrow; i++){
			List <Double> A = new ArrayList<Double>();
			List <Double> B = new ArrayList<Double>();
			for(j=0; j < ncol; j++ ){
				if(GT.get(i, j)==Alabel){
					A.add(CN.get(i, j));
				}else if(GT.get(i,j)==Blabel){
					B.add(CN.get(i, j));
				}
			}
			Acount.add(A.size());
			if(A.size()==0){
				Amean.add(Double.NaN);
				Asd.add(Double.NaN);
			}else{
				Amean.add(MyFunc.mean(A));
				Asd.add(MyFunc.sd(A));
			}
			Bcount.add(B.size());
			if(B.size()==0){
				Bmean.add(Double.NaN);
				Bsd.add(Double.NaN);
			}else{
				Bmean.add(MyFunc.mean(B));
				Bsd.add(MyFunc.sd(B));
			}
			if(A.size()<2 || B.size() <2){
				ttestP.add(1.0);
				continue;
			}
			double A2[] = new double[A.size()]; 
			double B2[] = new double[B.size()];
			for(j=0;j<A.size();j++){
				A2[j] = A.get(j);
			}
			for(j=0;j<B.size();j++){
				B2[j] = B.get(j);
			}
			double p = 1;
			try {
				p = TT.homoscedasticTTest(A2, B2);
			} catch (Exception e) {
				e.printStackTrace();
			} 
			if(Double.isNaN(p)){
				ttestP.add(1.0);
			}else{
				ttestP.add(p);
			}
		}
	}
	
	private List <Double> calculateTtestP(MyMat GT){
		List <Double> ttestP  = new ArrayList<Double>();
		int j,i,ncol,nrow;
		nrow = CN.rowSize();
		ncol = CN.colSize();
		for(i=0; i < nrow; i++){
			List <Double> A = new ArrayList<Double>();
			List <Double> B = new ArrayList<Double>();
			for(j=0; j < ncol; j++ ){
				if(GT.get(i, j)==Alabel){
					A.add(CN.get(i, j));
				}else if(GT.get(i,j)==Blabel){
					B.add(CN.get(i, j));
				}
			}
			if(A.size()<2 || B.size() <2){
				ttestP.add(1.0);
				continue;
			}
			double A2[] = new double[A.size()]; 
			double B2[] = new double[B.size()];
			for(j=0;j<A.size();j++){
				A2[j] = A.get(j);
			}
			for(j=0;j<B.size();j++){
				B2[j] = B.get(j);
			}
			double p = 1;
			try {
				p = TT.homoscedasticTTest(A2, B2);
			} catch (Exception e) {
				e.printStackTrace();
			} 
			if(Double.isNaN(p)){
				ttestP.add(1.0);
			}else{
				ttestP.add(p);
			}
		}
		return ttestP;	
	}
	
	private MyMat permutateGT(){
		MyMat permGT = new MyMat(GT.getRowNames(), GT.getColNames());
		int j,i,ncol,nrow;
		nrow = GT.rowSize();
		ncol = GT.colSize();
		for(i = 0; i < nrow; i++){
			List <Double> gt = GT.getRow(i);
			Collections.shuffle(gt);
			for(j = 0; j < ncol; j++){
				permGT.set(i, j, gt.get(j));
			}
		}
		return permGT;
	}
		
	private void perfomPermutation(){
		System.err.println("Perform Permutaion.....");
		int i;
		topPermutationP = new ArrayList<Double>();
		for(i = 0; i < numberOfPermutations; i++){
			System.err.println(i + "/" + numberOfPermutations);
			MyMat permGT = permutateGT();
			List <Double> P = calculateTtestP(permGT);
			topPermutationP.add(MyFunc.min(P));
		}
	}
	
	private void calculateCorrectedP(){
		System.err.println("Calculate corrected Pvalues.....");
		int i,j;
		int n = ttestP.size();
		int m = topPermutationP.size();
		correctedP = new ArrayList<Double>();
		for(i = 0; i < n; i++){
			double k = 0;
			for(j = 0; j < m; j++){
				if(ttestP.get(i) > topPermutationP.get(j)){
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
		field.add("Acount");
		field.add("Bcount");
		field.add("Amean");
		field.add("Bmean");
		field.add("Asd");
		field.add("Bsd");
		field.add("ttestP");
		field.add("correctedP");
		List <String> id = GT.getRowNames();
		result = new MyMat(id,field);
		int i;
		for(i=0;i<result.rowSize();i++){
			result.set(i, 0, Acount.get(i));
			result.set(i, 1, Bcount.get(i));
			result.set(i, 2, Amean.get(i));
			result.set(i, 3, Bmean.get(i));
			result.set(i, 4, Asd.get(i));
			result.set(i, 5, Bsd.get(i));
			result.set(i, 6, ttestP.get(i));
			result.set(i, 7, correctedP.get(i));
		}
		Map <String,Double> tmp = result.getColMap(6);
		id = MyFunc.sortKeysByAscendingOrderOfValues(tmp);
		result.reorderRows(id);
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
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] CNfile GTfile", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 2){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] CNfile GTfile", options);
			return;
		}
		MyMat CN = new MyMat(argList.get(0));
		MyMat GT = new MyMat(argList.get(1));
		AllelicImbalanceTest2 AIT = new AllelicImbalanceTest2(CN,GT);
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
