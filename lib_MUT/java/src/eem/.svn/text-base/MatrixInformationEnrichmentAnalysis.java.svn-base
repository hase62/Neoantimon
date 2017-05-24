package eem;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.*;

import no.uib.cipr.matrix.*;

import org.apache.commons.cli.*;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistributionImpl;

import sun.reflect.Reflection;
import utility.Dist;
import utility.ExtrapolatePvalueUsingExpRegression;
import utility.MyException;
import utility.MyFunc;
import utility.MyMat;

public class MatrixInformationEnrichmentAnalysis extends AbstractGeneSetAnalysis{
	MyMat Exp;
	Map <String, Double> Pvalue;
	Map <String, Double> informationContent; 
	int itrForPvalueCalculation = 1000;
	
	boolean itrAdjustment = false;
	boolean useExpReg = false;
	boolean gaussian  = false;
	
	NormalDistributionImpl ND = new NormalDistributionImpl();

	boolean useCovMatrix = false;
	boolean useSampleCovMatrix = false;
	
	public MatrixInformationEnrichmentAnalysis(MyMat E){
		Exp = E;
		allGenes = Exp.getRowNames();
		Pvalue = new HashMap<String, Double>();
		informationContent = new HashMap<String, Double>();
	}

	public void setItrForPvalueCalculation(int i){
		itrForPvalueCalculation = i;
	}
	
	//public void supressItrAdjustment(){
	//	itrAdjustment = 0;
	//}
	
	public void normalizeExpRows(){
		 Exp.normalizeRows();
	}
	
	
	public void useCovMatrix(){
		Exp = Exp.getCovMatForRow();
		useCovMatrix = true;
	}
	
	public void useSampleCovMatrix(){
		useSampleCovMatrix = true;
	}
	
	public  static double shannonEntropy(double[] v){
		   double e = 0;
		   double s = 0;
		   List <Double> v2 = new ArrayList<Double>();
		   for(int i = 0, n = v.length; i < n; i++){
			   //double tmp = Math.pow(v[i], 2);	
			   //s  +=  tmp;
			   //v2.add(tmp);
			   s  +=  v[i];
			   v2.add(v[i]);
			   //s  +=  Math.sqrt(v[i]);
			   //v2.add( Math.sqrt(v[i]));
			   
		   }
		   for(Double d: v2){
			   if(d > 0){
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
		  
		 double ic = 1-shannonEntropy(sv);
		 return  ic;
	}
	public void perform0 () {
		int n = seedGeneSets.size();
		int k=0;
		
		L:for(String genesetId: seedGeneSets.keySet()){
		k++;	
		try{
			List <String> geneset = seedGeneSets.get(genesetId);
			double infoCont = getInformationContent(geneset);
			
			if(Double.isNaN(infoCont)){
				System.err.println(genesetId +  "(" + k + "/"  + n + "):" + " failed! (infoCont is NaN)");
				informationContent.put(genesetId, -1.0);
				Pvalue.put(genesetId, -1.0);
				continue L;
			}
			
			List <Double> nullDist = new ArrayList<Double>();
			int i=0;
			informationContent.put(genesetId, infoCont);
			double  p = 0;
			while(nullDist.size() < itrForPvalueCalculation){ 
				double nullInfoCont = getInformationContent(MyFunc.sample(allGenes, geneset.size()));
				if(!Double.isNaN(nullInfoCont)){
					nullDist.add((nullInfoCont));
				}else{
					i++;
					if(i > itrForPvalueCalculation){
						System.err.println(genesetId +  "(" + k + "/"  + n + "):" + " failed! (unable to  generate null samples) ");
						informationContent.put(genesetId, -1.0);
						Pvalue.put(genesetId, -1.0);
						continue L;
					}				
				}
			}
			
			
			p = calculatePvalue2(infoCont, nullDist);
			
			
			if(Double.isInfinite(p) || Double.isNaN(p)){
				p = Double.MAX_VALUE;
			}
			
			double mean = MyFunc.mean(nullDist);
			double sd = MyFunc.sd(nullDist);
			System.err.println(nullDist);
			System.err.println(genesetId +  "(" + k + "/"  + n + "):"  + "\t" + p  + "\t" + infoCont + "\t" + mean + "\t" + sd);	
			
			Pvalue.put(genesetId, p);
			//System.err.println(genesetId +  "(" + k + "/"  + n + "):"  + "\t" + p);	
		}catch (Exception e) {
			System.err.println(genesetId +  "(" + k + "/"  + n + "):" + " failed!");
			e.printStackTrace();
			informationContent.put(genesetId, -1.0);
			Pvalue.put(genesetId, -1.0);
		}
		}	
	}
	

	public void perform1() {
		int n = seedGeneSets.size();
		int k=0;
		
		L:for(String genesetId: seedGeneSets.keySet()){
		k++;	
		try{
			List <String> geneset = seedGeneSets.get(genesetId);
			double infoCont = getInformationContent(geneset);
			
			if(Double.isNaN(infoCont)){
				System.err.println(genesetId +  "(" + k + "/"  + n + "):" + " failed! (infoCont is NaN)");
				informationContent.put(genesetId, -1.0);
				Pvalue.put(genesetId, -1.0);
				continue L;
			}
			
			
			
			
			List <Double> nullDist = new ArrayList<Double>();
			int i=0;
			int j = 0;
			informationContent.put(genesetId, infoCont);
			
			if(itrAdjustment){
				int d = 1;	
				double LogItrForPvalueCalculation = Math.log10(itrForPvalueCalculation);
					
				for(int I = 1 + d; I < LogItrForPvalueCalculation; I++){
					for(int m = (int) Math.pow(10, I); i<m;i++){
						double nullInfoCont = getInformationContent(MyFunc.sample(allGenes, geneset.size()));
						if(!Double.isNaN(nullInfoCont)){
							nullDist.add((nullInfoCont));
						}else{
							j++;
							if(j > itrForPvalueCalculation){
								System.err.println(genesetId +  "(" + k + "/"  + n + "):" + " failed! (unable to  generate null samples) ");
								informationContent.put(genesetId, -1.0);
								Pvalue.put(genesetId, -1.0);
								continue L;
							}				
						}
					}	
					double  p = 0;
					for(Double D: nullDist){
						if(D > infoCont){
							p++;
							}
					}
					if (p == 0){
						p = 1;
					}
					p /= nullDist.size();
					p = (p==1)?0:(-Math.log10(p));
						
						
					if(p  + d < I){
						Pvalue.put(genesetId, p);
						System.err.println(genesetId +  "(" + k + "/"  + n + "):"  + "\t" + p + "\t" + "( based on a null distribution of size " + nullDist.size() + " )");
						continue L;
					}
			
				}
			}
					
					
			for(; i< itrForPvalueCalculation; i++){ 
				double nullInfoCont = getInformationContent(MyFunc.sample(allGenes, geneset.size()));
				if(!Double.isNaN(nullInfoCont)){
					nullDist.add((nullInfoCont));
				}else{
					j++;
					if(j > itrForPvalueCalculation){
						System.err.println(genesetId +  "(" + k + "/"  + n + "):" + " failed! (unable to  generate null samples) ");
						informationContent.put(genesetId, -1.0);
						Pvalue.put(genesetId, -1.0);
						continue L;
					}				
				}
			}
			
			double  p = 0;
			for(Double d: nullDist){
				if(d > infoCont){
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
					 p = E.estimatePvalue(infoCont);
					 if(p > 1.0/nullDist.size() ){
						 System.err.println("not good:\t" + p + "\t" +  1/nullDist.size() );
					 }
				}else{
					p = 1.0/nullDist.size();
				}	
			}
			
			double m = MyFunc.mean(nullDist);
			double sd = MyFunc.sd(nullDist);
			
			p = (p==1)?0:(-Math.log10(p));
			if(Double.isInfinite(p)){
				p = Double.MAX_VALUE;
			}
			Pvalue.put(genesetId, p);
			//System.err.println(genesetId +  "(" + k + "/"  + n + "):"  + "\t" + p);	
			System.err.println(genesetId +  "(" + k + "/"  + n + "):"  + "\t" + p  + "\t" +  infoCont + "\t" +  m  +"\t" +  sd);	
			
		}catch (Exception e) {
			System.err.println(genesetId +  "(" + k + "/"  + n + "):" + " failed!");
			e.printStackTrace();
			informationContent.put(genesetId, -1.0);
			Pvalue.put(genesetId, -1.0);
		}
		}	
	}
	
	
	public void perform() {
		int n = seedGeneSets.size();
		int k=0;
		
		L:for(String genesetId: seedGeneSets.keySet()){
		k++;	
		try{
			List <String> geneset = seedGeneSets.get(genesetId);
			double infoCont = getInformationContent(geneset);
			
			if(Double.isNaN(infoCont)){
				System.err.println(genesetId +  "(" + k + "/"  + n + "):" + " failed! (infoCont is NaN)");
				informationContent.put(genesetId, -1.0);
				Pvalue.put(genesetId, -1.0);
				continue L;
			}
			
			informationContent.put(genesetId, infoCont);
			
			
			List <Double> nullDist = new ArrayList<Double>();
			int i=0;
			int j = 0;
			
			
			if(gaussian){
				double shapiroTestCutoff = 0.01;
				while(nullDist.size() < 1000){ 
					double nullInfoCont = getInformationContent(MyFunc.sample(allGenes, geneset.size()));
					if(!Double.isNaN(nullInfoCont)){
						nullDist.add((nullInfoCont));
						i++;
					}else{
						j++;
						if(j > itrForPvalueCalculation){
							System.err.println(genesetId +  "(" + k + "/"  + n + "):" + " failed! (unable to  generate null samples) ");
							informationContent.put(genesetId, -1.0);
							Pvalue.put(genesetId, -1.0);
							continue L;
						}				
					}
				}
				double shapiroP = MyFunc.shapiroTest(nullDist);
				if(shapiroP > shapiroTestCutoff){
					double p = calculatePvalue2(infoCont, nullDist);
					if(Double.isInfinite(p) || Double.isNaN(p)){
						p = Double.MAX_VALUE;
					}
					Pvalue.put(genesetId, p);
					System.err.println(genesetId +  "(" + k + "/"  + n + "):"  + "\t" + p + "\t" + "( by approximating nulldist as gaussian [shapiro p = " + shapiroP + "])");
					continue L;
				}
			}
			
			
			
			
			if(itrAdjustment){
				int d = 1;	
				double LogItrForPvalueCalculation = Math.log10(itrForPvalueCalculation);
					
				for(int I = 1 + d; I < LogItrForPvalueCalculation; I++){
					for(int m = (int) Math.pow(10, I); i<m;i++){
						double nullInfoCont = getInformationContent(MyFunc.sample(allGenes, geneset.size()));
						if(!Double.isNaN(nullInfoCont)){
							nullDist.add((nullInfoCont));
						}else{
							j++;
							if(j > itrForPvalueCalculation){
								System.err.println(genesetId +  "(" + k + "/"  + n + "):" + " failed! (unable to  generate null samples) ");
								informationContent.put(genesetId, -1.0);
								Pvalue.put(genesetId, -1.0);
								continue L;
							}				
						}
					}	
					double  p = 0;
					for(Double D: nullDist){
						if(D > infoCont){
							p++;
							}
					}
					if (p == 0){
						p = 1;
					}
					p /= nullDist.size();
					p = (p==1)?0:(-Math.log10(p));
						
						
					if(p  + d < I){
						Pvalue.put(genesetId, p);
						System.err.println(genesetId +  "(" + k + "/"  + n + "):"  + "\t" + p + "\t" + "( based on a null distribution of size " + nullDist.size() + " )");
						continue L;
					}
			
				}
			}
					
					
			for(; i< itrForPvalueCalculation; i++){ 
				double nullInfoCont = getInformationContent(MyFunc.sample(allGenes, geneset.size()));
				if(!Double.isNaN(nullInfoCont)){
					nullDist.add((nullInfoCont));
				}else{
					j++;
					if(j > itrForPvalueCalculation){
						System.err.println(genesetId +  "(" + k + "/"  + n + "):" + " failed! (unable to  generate null samples) ");
						informationContent.put(genesetId, -1.0);
						Pvalue.put(genesetId, -1.0);
						continue L;
					}				
				}
			}
			
			double  p = 0;
			for(Double d: nullDist){
				if(d > infoCont){
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
					 p = E.estimatePvalue(infoCont);
					 if(p > 1.0/nullDist.size() ){
						 System.err.println("not good:\t" + p + "\t" +  1/nullDist.size() );
					 }
				}else{
					p = 1.0/nullDist.size();
				}	
			}
			
			double m = MyFunc.mean(nullDist);
			double sd = MyFunc.sd(nullDist);
			
			p = (p==1)?0:(-Math.log10(p));
			if(Double.isInfinite(p)){
				p = Double.MAX_VALUE;
			}
			Pvalue.put(genesetId, p);
			//System.err.println(genesetId +  "(" + k + "/"  + n + "):"  + "\t" + p);	
			System.err.println(genesetId +  "(" + k + "/"  + n + "):"  + "\t" + p  + "\t" +  infoCont + "\t" +  m  +"\t" +  sd);	
			
		}catch (Exception e) {
			System.err.println(genesetId +  "(" + k + "/"  + n + "):" + " failed!");
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
		options.addOption("i", "pitr", true,  "itrForPvalueCalculation");
		options.addOption("o", "outfile", true, "outFile");
		options.addOption("m", "mingeneset", true, "minGeneSetSize");
		options.addOption("M", "maxgeneset", true, "maxGeneSetSize");
		options.addOption("c", "covmat", false, "useCovMatrix");
		options.addOption("s", "smpcovmat", false, "useSampleCovMatrix");
		options.addOption("n", "rownorm", false, "normalizeExpRows");
		options.addOption("I", "itrad", false,  "doItrAdjustment");
		options.addOption("e", "expreg", false, "useExpReg");
		options.addOption("G", "gaussian", false, "approximate null dist as gaussian");
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
		MatrixInformationEnrichmentAnalysis  GSA = new  MatrixInformationEnrichmentAnalysis (new MyMat(argList.get(0)));	
		GSA.setGeneSets(MyFunc.readGeneSetFromGmtFile(argList.get(1)));
		if(commandLine.hasOption("G")){
			GSA.gaussian = true;
		}
		if(commandLine.hasOption("n")){
			GSA.normalizeExpRows();
		}
		if(commandLine.hasOption("c")){
			GSA.useCovMatrix();
		}
		if(commandLine.hasOption("s")){
			GSA.useSampleCovMatrix();
		}
		if(commandLine.hasOption("i")){
			GSA.setItrForPvalueCalculation(Integer.valueOf(commandLine.getOptionValue("i")));
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
				 System.out.println(s + "\t" + Pvalues.get(s) );				 
			 }
			 System.out.close();
		}
	}
	
	
	
}
