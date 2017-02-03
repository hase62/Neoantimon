package eem;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.*;

import no.uib.cipr.matrix.*;

import org.apache.commons.cli.*;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistributionImpl;

import utility.Dist;
import utility.MyException;
import utility.MyFunc;
import utility.MyMat;

public class SingularValueEntoropyAnalysis  implements GeneSetAnalysis{
	MyMat Exp;
	List <String> allGenes; 
	Map<String, List<String>> seedGeneSets;
	int maxSeedGeneSize = 2000;
	int minSeedGeneSize = 10;
	Map <String, Double> Pvalue;
	Map <String, Double> informationContent; 
	int itrForPvalueCalculation = 300;
	int itrAdjustment = 1;
	
	NormalDistributionImpl ND = new NormalDistributionImpl();

	public SingularValueEntoropyAnalysis (MyMat E){
		Exp = E;
		allGenes = Exp.getRowNames();
		Pvalue = new HashMap<String, Double>();
		informationContent = new HashMap<String, Double>();
	}
	public void setGeneSets(Map<String, List<String>> geneSets){
		seedGeneSets = new HashMap<String, List<String>>();
		for(Map.Entry<String, List<String>> e: geneSets.entrySet() ){
			List <String> tmp = MyFunc.isect(e.getValue(), allGenes);
			if(tmp.size() < minSeedGeneSize || tmp.size() > maxSeedGeneSize){
				continue; 
			}
			seedGeneSets.put(e.getKey(), tmp);	
		}
	}
	
	public void setMaxGeneSetSize(int i){
		maxSeedGeneSize = i;
	}
	public void setMinGeneSetSize(int i){
		minSeedGeneSize = i;
	}
	public void setItrForPvalueCalculation(int i){
		itrForPvalueCalculation = i;
	}
	
	public void supressItrAdjustment(){
		itrAdjustment = 0;
	}
	
	
	private  static double shannonEntropy(double[] v){
		   double e = 0;
		   double s = 0;
		   List <Double> v2 = new ArrayList<Double>();
		   for(int i = 0, n = v.length; i < n; i++){
			   //double tmp = Math.pow(v[i], 2);	
			   //s  +=  tmp;
			   //v2.add(tmp);
			   s  +=  v[i];
			   v2.add(v[i]);
			   
			   
		   }
		   for(Double d: v2){
			   if(d != 0){
			    double p = d/s;
			    e += p*Math.log(p);
			   }
		   }
		   return - e/Math.log(v.length);   
	 }
	
	private double calculatePvalue(double score, List<Double> nullDist){
		double p = 0;
		for(Double d: nullDist){
			if(d > score){
				p++;
			}
		}
		if (p == 0){
			p = 1;
		}
		p /= nullDist.size();
		return -Math.log10(p);
	}
	
	private  double calculatePvalue2(double score, List<Double> nullDist) throws MathException{
		double nullMean = MyFunc.mean(nullDist);
		double nullSd = MyFunc.sd(nullDist);
		
		double z = (score - nullMean) /nullSd; 

		return -Math.log10(1-ND.cumulativeProbability(z));
		
		
	}
	
	 
	public double getInformationContent(List <String>  geneset) throws NotConvergedException{
		 DenseMatrix M = Exp.getSubMatByRow(geneset).getMTJDenseMatrix();
		 SVD svd = SVD.factorize(M);
		 double[] sv =  svd.getS();
		return  1-shannonEntropy(sv);
	}
	public void perform () {
		L:for(String genesetId: seedGeneSets.keySet()){
		try{
			List <String> geneset = seedGeneSets.get(genesetId);
			double infoCont = getInformationContent(geneset);
			List <Double> nullDist = new ArrayList<Double>();
			int i=0;
			informationContent.put(genesetId, infoCont);
			double  p = 0;
			int d = itrAdjustment;	
			double LogItrForPvalueCalculation = Math.log10(itrForPvalueCalculation);
			
			
			if(d > 0){
				for(int I = 1 + d; I < LogItrForPvalueCalculation; I++){
					for(int n = (int) Math.pow(10, I); i<n;i++){
						nullDist.add(getInformationContent(MyFunc.sample(allGenes, geneset.size())));
					}
					p = calculatePvalue(infoCont, nullDist);
					if(p  + d < I){
						Pvalue.put(genesetId, p);
						System.err.println(genesetId + "\t" + p + "\t" + infoCont +  "\t" + Math.log10(nullDist.size()));
						continue L;
					}
				}
			}
			for(; i< itrForPvalueCalculation; i++){ 
				nullDist.add(getInformationContent(MyFunc.sample(allGenes, geneset.size())));
			}
			p = calculatePvalue2(infoCont, nullDist);
			Pvalue.put(genesetId, p);
			System.err.println(genesetId + "\t" + p + "\t" + infoCont  +  "\t" + Math.log10(nullDist.size()));	
		}catch (Exception e) {
			System.err.println(genesetId + ": failed!");
			e.printStackTrace();
			informationContent.put(genesetId, -1.0);
			Pvalue.put(genesetId, -1.0);
		}
		}	
	}
	

	
	public Map <String, Double> getPvalues(){
		return Pvalue;
	}

	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("i", "pitr", true,  "itrForPvalueCalculation (log10 scale)");
		options.addOption("o", "outfile", true, "outFile");
		options.addOption("m", "mingeneset", true, "minGeneSetSize");
		options.addOption("M", "maxgeneset", true, "maxGeneSetSize");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp("CoherenceBasedEEMsearch [options] expfile gmtfile", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 2){
			formatter.printHelp("CoherenceBasedEEMsearch [options] expfile gmtfile", options);
			return;
		}
		SingularValueEntoropyAnalysis GSECA = new  SingularValueEntoropyAnalysis(new MyMat(argList.get(0)));	
		GSECA.setGeneSets(MyFunc.readGeneSetFromGmtFile(argList.get(1)));
		
		if(commandLine.hasOption("i")){
			GSECA.setItrForPvalueCalculation(Integer.valueOf(commandLine.getOptionValue("i")));
		}
		
		if(commandLine.hasOption("m")){
			GSECA.setMinGeneSetSize(Integer.valueOf(commandLine.getOptionValue("m")));
		}
		if(commandLine.hasOption("M")){
			GSECA.setMaxGeneSetSize(Integer.valueOf(commandLine.getOptionValue("M")));
		}
		GSECA.perform();
		Map <String, Double> Pvalues = GSECA.getPvalues();
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
