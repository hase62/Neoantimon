package eem;


import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.*;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;


import sun.reflect.Reflection;
import utility.*;

public class GeneSetExpressionCoherenceAnalysis  extends AbstractGeneSetAnalysis{
	MyMat originalExp;
	MyMat Exp;
	Dist Cor;
	Map <String, Double> Pvalue;
	int itrForPvalueCalculation = 1000;
	boolean absoluteCorrelation = false;
	boolean itrAdjustment = false;
	boolean useExpReg = false;
	
	interface DistFunc {
		double get(String a, String b);
	}
	private class Correlation implements DistFunc{
		public double get(String a, String b){
			return  Cor.get(a,b);
		}
	}
	private class AbsoluteCorrelation implements DistFunc{
		public double get(String a, String b){
			return Math.abs(Cor.get(a,b));
		}
	}
	DistFunc distfunc; 
	public GeneSetExpressionCoherenceAnalysis(MyMat E){
		Exp = new MyMat(E);
		Exp.normalizeRows();
		Cor = null;
		allGenes = Exp.getRowNames();
		Pvalue = new HashMap<String, Double>();
		distfunc =new Correlation();
	}
	
	public void useAbsoluteCorrelation(){
		absoluteCorrelation = true;
		distfunc =new AbsoluteCorrelation();
	}
	
	
	public void calculateCor(){
		Cor = new Dist(Exp, 'C');	
	}
	public void setCor(Dist D){
		if(D.getNames().containsAll(Exp.getRowNames())){
			throw new MyException("Correlation data must include data for all genes in expression data");
		}
		Cor = D;
	}
	
	public void setItrForPvalueCalculation(int i){
		itrForPvalueCalculation = i;
	}
	
	public void supressItrAdjustment(){
		itrAdjustment = false;
	}
	public double getCorMean(List <String>  geneset){
		List <Double> v = new ArrayList<Double>();
		int i,j, n;
		for(i=0, n = geneset.size();i<n;i++){
			for(j=0;j<i;j++){
				v.add(distfunc.get(geneset.get(i), geneset.get(j)));
			}
		}	
		return  MyFunc.mean(v);
	}
	
	
	private void  calculatePvalue(){
		int n = seedGeneSets.size();
		int k=0;
		L:for(String genesetId: seedGeneSets.keySet()){
			k++;
			List <String> geneset = seedGeneSets.get(genesetId);
			double corMean = getCorMean(geneset);
			List <Double> nullDist = new ArrayList<Double>();
			int i=0;
			
			if(itrAdjustment){
				int d = 1;	
				double LogItrForPvalueCalculation = Math.log10(itrForPvalueCalculation);
				
				
				for(int I = 1 + d; I < LogItrForPvalueCalculation; I++){
					for(int m = (int) Math.pow(10, I); i<m;i++){
						nullDist.add(getCorMean(MyFunc.sample(allGenes, geneset.size())));
					}
					double  p = 0;
					for(Double D: nullDist){
						if(D > corMean){
							p++;
							}
					}
					if (p == 0){
						p = 1;
					}
					p /= nullDist.size();
					p = -Math.log10(p);
						
						
					if(p  + d < I){
						Pvalue.put(genesetId, p);
						System.err.println(genesetId +  "(" + k + "/"  + n + "):"  + "\t" + p + "\t" + "( based on a null distribution of size " + nullDist.size() + " )");
						continue L;
					}
			
				}
			}
			for(; i< itrForPvalueCalculation; i++){ 
				nullDist.add(getCorMean(MyFunc.sample(allGenes, geneset.size())));
			}	
			double  p = 0;
			for(Double d: nullDist){
				if(d > corMean){
					p++;
				}
			}
			p /= nullDist.size();
			if (p == 0){
				if(useExpReg){
					ExtrapolatePvalueUsingExpRegression E = new ExtrapolatePvalueUsingExpRegression(nullDist);
					double upperLimit = MyFunc.mean(nullDist)*50;
					 E.setXrange(0.75, 0.95);
					 E.setPoweRange(1, 3);
					 E.setUpperLimit(upperLimit);
					 E.regress();
					 p = E.estimatePvalue(corMean);
					 if(p > 1.0/nullDist.size() ){
						 System.err.println("not good:\t" + p + "\t" +  1/nullDist.size() );
					 }
				}else{
					p = 1.0/nullDist.size();
				}	
			}
			
			double m = MyFunc.mean(nullDist);
			double sd = MyFunc.sd(nullDist);
			
			
			p = -Math.log10(p);
			Pvalue.put(genesetId, p);
			//System.err.println(genesetId +  "(" + k + "/"  + n + "):"  + "\t" + p );
			System.err.println(genesetId +  "(" + k + "/"  + n + "):"  + "\t" + p  + "\t" +  corMean + "\t" +  m  +"\t" +  sd);	
			
		}
	}
	

	public void perform(){
		if(Cor == null){
			calculateCor();
		}
		calculatePvalue();
	}
	
	public Map <String, Double> getPvalues(){
		return Pvalue;
	}

	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("c", "cor", true,  "correlationFile");
		options.addOption("i", "itr", true,  "itrForPvalueCalculation");
		options.addOption("I", "itrad", false,  "doItrAdjustment");
		options.addOption("a", "absolute", false, "useAbsoluteCorrelation");
		options.addOption("o", "outfile", true, "outFile");
		options.addOption("m", "mingeneset", true, "minGeneSetSize");
		options.addOption("M", "maxgeneset", true, "maxGeneSetSize");
		options.addOption("e", "expreg", false, "useExpReg");
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
		
		GeneSetExpressionCoherenceAnalysis GSA = new  GeneSetExpressionCoherenceAnalysis(new MyMat(argList.get(0)));	
		GSA.setGeneSets(MyFunc.readGeneSetFromGmtFile(argList.get(1)));
		if(commandLine.hasOption("c")){
			GSA.setCor(new Dist(commandLine.getOptionValue("c")));
		}
		if(commandLine.hasOption("i")){
			GSA.setItrForPvalueCalculation(Integer.valueOf(commandLine.getOptionValue("i")));
		}
		if(commandLine.hasOption("a")){
			GSA.useAbsoluteCorrelation();
		}
		if(commandLine.hasOption("m")){
			GSA.setMinGeneSetSize(Integer.valueOf(commandLine.getOptionValue("m")));
		}
		if(commandLine.hasOption("M")){
			GSA.setMaxGeneSetSize(Integer.valueOf(commandLine.getOptionValue("M")));
		}
		if(commandLine.hasOption("I")){
			GSA.itrAdjustment = true;
		}
		if(commandLine.hasOption("e")){
			GSA.useExpReg = true;
		}
		GSA.perform();
		Map <String, Double> Pvalues = GSA.getPvalues();
		if(commandLine.hasOption("o")){
			PrintWriter os = new PrintWriter(new FileWriter(commandLine.getOptionValue("o")));
			 for(String s: MyFunc.sortKeysByDescendingOrderOfValues(Pvalues)){
				 os.println(s + "\t" + Pvalues.get(s) );				 
			 }
			 os.close();
		}else{
			 for(String s: MyFunc.sortKeysByDescendingOrderOfValues(Pvalues)){
				 System.out.println(s + "\t" + Pvalues.get(s) + "\t" );				 
			 }
			 System.out.close();
		}
	}
	
	
	
}
