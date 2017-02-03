package eem;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

import sun.reflect.Reflection;
import utility.ExpRegression;
import utility.MyException;
import utility.MyFunc;
import utility.MyMat;

public class SMIEAtest extends MatrixInformationEnrichmentAnalysis {

	
	int geneSetSize = 100;
	int n = 100;
	
	
	public SMIEAtest(MyMat E) {
		super(E);
	}

	
	public void perform(){
		
		if(geneSetSize >= Exp.rowSize()){
			n = 1;
			geneSetSize = Exp.rowSize();
		}
	
		List <Double> IC = new ArrayList<Double>();
		for(int i = 0; i < n ; i++){
			List <String> gs = MyFunc.sample(Exp.getRowNames(), geneSetSize);
			try {
				IC.add(getInformationContent(gs));
			} catch (NotConvergedException e) {
				e.printStackTrace();
			}
		}
		
	
		
		
		double mean  = MyFunc.mean(IC);
		double sd = MyFunc.sd(IC);
		System.out.println(geneSetSize + "\t" + mean + "\t" + sd);
		
	}
	
	
	public void perform3(){
		
		if(geneSetSize >= Exp.rowSize()){
			n = 1;
			geneSetSize = Exp.rowSize();
		}
	
		List <Double> IC = new ArrayList<Double>();
		for(int i = 0; i < n ; i++){
			List <String> gs = MyFunc.sample(Exp.getRowNames(), geneSetSize);
			try {
				IC.add(getInformationContent(gs));
			} catch (NotConvergedException e) {
				e.printStackTrace();
			}
		}
		
		
		List <Double> nullDist = IC;
		MyFunc.Density D = new MyFunc.Density(nullDist); 
		double minX = MyFunc.percentile(nullDist, 0.6);
		double maxX = MyFunc.percentile(nullDist, 0.95);
		double d = (maxX - minX) / 500;
		List <Double> X = new ArrayList<Double>();
		double x;
		for(x = minX; x <= maxX; x += d ){
			X.add(x);
		}
		List <Double> Y = D.estimate(X);
		ExpRegression ER  = new ExpRegression();
		ER.setTrainingData(X, Y);
		List <Double> powers = new ArrayList<Double>();
		double minPower = 1.0;
		double maxPower = 3.0;
		for( double p = minPower; p <= maxPower; p += 0.1){
			powers.add(p);
		}
		try {
			ER.optimizePowerOfX(powers);
		} catch (MyException e) {
			e.printStackTrace();
		}
		
		//System.err.println(ER); 
		System.err.println(ER.getPowerOfX()); 
		double upperLimit = MyFunc.mean(nullDist)*50;
		
		for(double p = 0.9; p < 0.990001;  p += 0.01){
			double  t   = MyFunc.percentile(nullDist, p);
			double P = ER.integrate(t, upperLimit);
			System.err.println(1-p + "\t" + P + "\t" + t + "\t" + upperLimit);
		}
		System.out.println(MyFunc.join("\n", IC));	
	}
	
public void perform2(){
		
		if(geneSetSize >= Exp.rowSize()){
			n = 1;
			geneSetSize = Exp.rowSize();
		}
	
		List <Double> IC = new ArrayList<Double>();
		for(int i = 0; i < n ; i++){
			List <String> gs = MyFunc.sample(Exp.getRowNames(), geneSetSize);
			
			
			try {
				IC.add(getInformationContent(gs));
			} catch (NotConvergedException e) {
				e.printStackTrace();
			}
		}
		//System.out.println(MyFunc.shapiroTest(IC));
		System.out.println(MyFunc.join("\n", IC));	
	}



public void perform4(){
	
	if(geneSetSize >= Exp.rowSize()){
		n = 1;
		geneSetSize = Exp.rowSize();
	}

	List <Double> IC = new ArrayList<Double>();
	
	Map <String, String> out = new HashMap<String, String>();
	Map <String, Double> p = new HashMap<String, Double>();
	for(int i = 0; i < n ; i++){
		List <String> gs = MyFunc.sample(Exp.getRowNames(), geneSetSize);
		String id = "gs" + i;
		
		
		try {
			double ic = getInformationContent(gs);
			p.put(id, ic);
			out.put(id, id + "\t" + ic + "\t" + MyFunc.join("\t", gs));
		} catch (NotConvergedException e) {
			e.printStackTrace();
		}
	}
	for(String s: MyFunc.sortKeysByDescendingOrderOfValues(p)){
		System.out.println(out.get(s));
	}	
}
	
public ArrayList<Double> getSV(List <String>  geneset) throws NotConvergedException{
	
	 DenseMatrix M;
	 if(useCovMatrix){
		 M = Exp.getSubMatrix(geneset, geneset).getMTJDenseMatrix();			 
	 }else if (useSampleCovMatrix){
		 M = Exp.getSubMatByRow(geneset).getCovMatForCol().getMTJDenseMatrix();
	 }else{
		 M = Exp.getSubMatByRow(geneset).getMTJDenseMatrix();
	 }
	 SVD svd = SVD.factorize(M);
	 double[] sv =  svd.getS();
	 
	 ArrayList<Double> SV = new ArrayList<Double>();
	 for(int i = 0; i < sv.length; i++){
		 SV.add(sv[i]);
	 }
	 return  SV;
}

public void perform5(){
	
	if(geneSetSize >= Exp.rowSize()){
		n = 1;
		geneSetSize = Exp.rowSize();
	}

	List <Double> IC = new ArrayList<Double>();
	
	Map <String, String> out = new HashMap<String, String>();
	Map <String, Double> p = new HashMap<String, Double>();
	for(int i = 0; i < n ; i++){
		List <String> gs = MyFunc.sample(Exp.getRowNames(), geneSetSize);
		String id = "gs" + i;
		
		
		try {
			double ic = getInformationContent(gs);
			p.put(id, ic);
			out.put(id, ic + "\t" + MyFunc.join("\t", getSV(gs)));
		} catch (NotConvergedException e) {
			e.printStackTrace();
		}
	}
	for(String s: MyFunc.sortKeysByDescendingOrderOfValues(p)){
		System.out.println(out.get(s));
	}	
}
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("n", "number", true,  "sample number");
		options.addOption("g", "gs", true, "gene set size");
		options.addOption("c", "covmat", false, "useCovMatrix");
		options.addOption("a", "a", false, "useCovMatrix");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] expfile", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 1){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] expfile", options);
			return;
		}
		
	
		SMIEAtest S = new SMIEAtest(new MyMat(argList.get(0)));
		
		S.normalizeExpRows();
		if(commandLine.hasOption("c")){
			S.useCovMatrix();
		}
		if(commandLine.hasOption("n")){
			S.n = Integer.valueOf(commandLine.getOptionValue("n"));
		}
		if(commandLine.hasOption("g")){
			S.geneSetSize = Integer.valueOf(commandLine.getOptionValue("g"));
		}
		if(commandLine.hasOption("a")){
			S.perform4();
		}else{
		S.perform5();
		}
	}
	

	
	
	
}
