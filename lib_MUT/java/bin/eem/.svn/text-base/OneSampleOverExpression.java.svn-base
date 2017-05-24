package eem;


import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.*;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ChiSquaredDistribution;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;


import sun.reflect.Reflection;
import utility.*;

public class OneSampleOverExpression extends AbstractGeneSetAnalysis{
	MyMat Exp;
	double overExpressionCutoff;
	int itrNullPdistCalculation = 0;
	double covSum = 0; 
	Map <String, Double> Pvalue;
	MyMat PvalueForSample;
	
	boolean printPvalueForSample = false;
	boolean absolute = false;
	boolean underExpression = false;
	boolean takeBestPvalue = false;
	
	public OneSampleOverExpression(MyMat E){
		Exp = new MyMat(E);
		Exp.normalizeRows();
		allGenes = Exp.getRowNames();
		overExpressionCutoff = 0.05;
		Pvalue = new HashMap<String, Double>();
	}
	
	
	public void perform(){
		calculateNullPCov();
		calculatePvalue();
	}
	public void setOverExpressionCutoff(double d){
		if(d < 0){
			overExpressionCutoff = d;
		}else{
			overExpressionCutoff = d/100;
		}
		
	}
	public void findAbsoluteOverExpression(){
		absolute = true;
	}
	public void findUnderExpressionAlso(){
		underExpression = true;
	}
	
	
	public void calculateNullPCov(){
		if(itrNullPdistCalculation> 0){
			System.err.println("calculating null P value covatiance ...");	
			int colsize =  underExpression?Exp.colSize()*2:Exp.colSize();
			MyMat NullP = new MyMat(itrNullPdistCalculation, colsize);
			List <String> geneSetIDs = MyFunc.sampleWithReplacement(new ArrayList<String>(seedGeneSets.keySet()), itrNullPdistCalculation);
			for(int i=0; i < itrNullPdistCalculation; i++){
				System.err.println( i+1 + "/"  + itrNullPdistCalculation );	
				List <String> nullGeneSet = MyFunc.sample(Exp.getRowNames(), seedGeneSets.get(geneSetIDs.get(i)).size());	
				if(underExpression){
					for(int j = 0; j < Exp.colSize(); j++){
						NullP.set(i,j*2, getOverexpressPvalue(nullGeneSet, Exp.getColNames().get(j)));
						NullP.set(i,j*2+1,getUnderexpressPvalue(nullGeneSet, Exp.getColNames().get(j)));
					}
				}else{
					for(int j = 0; j < Exp.colSize(); j++){
						NullP.set(i,j, getOverexpressPvalue(nullGeneSet, Exp.getColNames().get(j)));
					}	
				}
			}	
			MyMat CovMat = NullP.getCovMatForCol();
			for(int i = 0; i < CovMat.colSize(); i++){
				for(int j = 0; j < CovMat.rowSize(); j++){
					covSum += CovMat.get(i, j);
				}
			}
			covSum *= Math.pow(Math.log(10), 2);
		}		
	}
	
	
	private void  calculatePvalue(){
		if(absolute){
			for(int i=0;i<Exp.rowSize();i++){
				for(int j=0;j<Exp.colSize();j++){
					Exp.set(i, j, Math.abs(Exp.get(i, j)));
				}	
			}
		}
		int n = seedGeneSets.size();
		int k=0;
		
		if(!underExpression){
			PvalueForSample = new MyMat(new ArrayList<String>(seedGeneSets.keySet()),Exp.getColNames());
		}else{
			List <String> tmp = new ArrayList<String>();	
			for(String colname: Exp.getColNames()){
				tmp.add(colname + "_up");
				tmp.add(colname + "_down");
			}
			PvalueForSample = new MyMat(new ArrayList<String>(seedGeneSets.keySet()),tmp);
		}
		
		for(String geneset: seedGeneSets.keySet()){
			k++;
			List <Double> P = new ArrayList<Double>();
			for(String colname: Exp.getColNames()){
				
				
				P.add(getOverexpressPvalue(new ArrayList<String>(seedGeneSets.get(geneset)), colname));
				
				double upP = getOverexpressPvalue(new ArrayList<String>(seedGeneSets.get(geneset)), colname);
				P.add(upP);
				if(underExpression){
					double downP = getUnderexpressPvalue(new ArrayList<String>(seedGeneSets.get(geneset)), colname);
					P.add(downP);
					
					PvalueForSample.set(geneset, colname + "_up",  upP==1?0:-Math.log10(upP));
					PvalueForSample.set(geneset, colname + "_down", downP==1?0:-Math.log10(downP));
				}else{
					PvalueForSample.set(geneset, colname, upP==1?0:-Math.log10(upP));
				}		
			}
			
			if(printPvalueForSample){ 
				System.err.println(geneset +  "(" + k + "/"  + n + "):"  + "\t" + PvalueForSample.getRow(geneset));	
			}else{
				double p;
				if(takeBestPvalue){
					p = MyFunc.min(P);
				}else{
					p = combinePvalues(P, covSum);
				}
				p = p==1?0:-Math.log10(p);
				p = Double.isInfinite(p) ? Double.MAX_VALUE : p;
				Pvalue.put(geneset, p);
				System.err.println(geneset +  "(" + k + "/"  + n + "):"  + "\t" + p);	
			}
		}
	}
	
	 public static double combinePvalues(List <Double> Pvalue, double covSum){
		   double  chi = 0;
		   double df;
		   for(Double p : Pvalue){
			   chi += -Math.log(p);
		   }
		   if(covSum == 0){ 
			   chi *=2;
			   df = 2*Pvalue.size();
		   }else{
			   chi *= 2*Pvalue.size()/(covSum);
			   df = 2*Math.pow(Pvalue.size(), 2)/(covSum);
		   }
		   ChiSquaredDistribution chidist =  new ChiSquaredDistributionImpl(df);
		   double p = 0;
		   try {
			   p = 1 - chidist.cumulativeProbability(chi);
			   return p;
			   } catch (MathException e) {
				   e.printStackTrace();   
			   }
			   return p;
	   }
	

	public double getOverexpressPvalue(List <String> seedGenes, String sample){
		Map <String, Double> exp = Exp.getColMap(sample);
		double percentile  = MyFunc.percentile(new ArrayList<Double>(exp.values()), 1-overExpressionCutoff);
		List <String> overExpressedGenes = new ArrayList<String>((MyFunc.getSubMapWithValueAboveCutoff(exp, percentile)).keySet());
		int isect = MyFunc.isect(seedGenes, overExpressedGenes).size();
		double p = isect==0?1:MyFunc.calculatePvalueForSetOverlap(allGenes.size(), overExpressedGenes.size(), seedGenes.size(), isect);
		return p;
		
	}
	
	public double getUnderexpressPvalue(List <String> seedGenes,  String sample){
		Map <String, Double> exp = Exp.getColMap(sample);
		double percentile  = MyFunc.percentile(new ArrayList<Double>(exp.values()), overExpressionCutoff);
		List <String> underExpressedGenes = new ArrayList<String>((MyFunc.getSubMapWithValueBelowCutoff(exp, percentile)).keySet());
		int isect = MyFunc.isect(seedGenes, underExpressedGenes).size();
		double p = isect==0?1:MyFunc.calculatePvalueForSetOverlap(allGenes.size(), underExpressedGenes.size(), seedGenes.size(), isect);
		return p;
	}
	
	
	
	public Map <String, Double> getPvalues(){
		return Pvalue;
	}
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("o", "outfile", true, "outFile");
		options.addOption("m", "mingeneset", true, "minGeneSetSize");
		options.addOption("M", "maxgeneset", true, "maxGeneSetSize");
		options.addOption("a", "absolute", false, "find absolute overexpression");
		options.addOption("u", "underexpression", false, "find underexpression also");
		options.addOption("c", "cutoff", true, "overexpression cutoff");
		options.addOption("p", "pcovsum", true, "unll pvalue covariance sum");
		options.addOption("i", "itr", true, "# of null pvalue to be calculated");
		options.addOption("s", "sample", false, "print p-value for each sample");
		options.addOption("b", "bestp", false, "take the best p-value");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] expfile gmtfile", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 2){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] expfile gmtfile", options);
			return;
		}
		 OneSampleOverExpression GSA = new  OneSampleOverExpression(new MyMat(argList.get(0)));	
		GSA.setGeneSets(MyFunc.readGeneSetFromGmtFile(argList.get(1)));
		if(commandLine.hasOption("a")){
			GSA.findAbsoluteOverExpression();
		}
		if(commandLine.hasOption("u")){
			GSA.findUnderExpressionAlso();
		}
		if(commandLine.hasOption("m")){
			GSA.setMinGeneSetSize(Integer.valueOf(commandLine.getOptionValue("m")));
		}
		if(commandLine.hasOption("M")){
			GSA.setMaxGeneSetSize(Integer.valueOf(commandLine.getOptionValue("M")));
		}
		if(commandLine.hasOption("c")){
			GSA.setOverExpressionCutoff(Double.valueOf(commandLine.getOptionValue("c")));
		}
		if(commandLine.hasOption("i")){
			GSA.itrNullPdistCalculation = Integer.valueOf(commandLine.getOptionValue("i"));
		}
		if(commandLine.hasOption("p")){
			GSA.covSum = Double.valueOf(commandLine.getOptionValue("p"));
			GSA.itrNullPdistCalculation = 0;
		}
		if(commandLine.hasOption("s")){
			GSA.printPvalueForSample = true;
			GSA.itrNullPdistCalculation = 0;
		}
		if(commandLine.hasOption("b")){
			GSA.takeBestPvalue = true;
			GSA.itrNullPdistCalculation = 0;
		}
		GSA.perform();
		Map <String, Double> Pvalues = GSA.getPvalues();
		if(commandLine.hasOption("o")){
			PrintWriter os = new PrintWriter(new FileWriter(commandLine.getOptionValue("o")));
			if(GSA.printPvalueForSample){
				os.print(GSA.PvalueForSample);	
			}else{
			for(String s: MyFunc.sortKeysByDescendingOrderOfValues(Pvalues)){
				 os.println(s + "\t" + Pvalues.get(s));				 
			 }
			}
			 os.close();
		}else{
			if(GSA.printPvalueForSample){
				 System.out.print(GSA.PvalueForSample);	
			}else{
			 for(String s: MyFunc.sortKeysByDescendingOrderOfValues(Pvalues)){
				 System.out.println(s + "\t" + Pvalues.get(s));				 
			 }
			}
			 System.out.close();
		}
		
	}
	
	
	
}
