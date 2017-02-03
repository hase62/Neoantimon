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
import org.apache.commons.math.distribution.NormalDistributionImpl;


import sun.reflect.Reflection;
import utility.MyFunc;
import utility.MyMat;

public class OneSampleOverExpression2 extends AbstractGeneSetAnalysis{
	MyMat Exp;
	MyMat ScoreMat;	
	Map <String, Double> ScoreMap;
	List <String> sample;
	int itrForZvalueCalculation = 100;
	NormalDistributionImpl ND = new NormalDistributionImpl();
	double covSum = 0;
	double ceilingScore  = 30;
	boolean useZscore = false;
	boolean printScoresForEachSample = false;
	public OneSampleOverExpression2(MyMat E){
		Exp  = new MyMat(E);	
		Exp.normalizeCols();
		Exp.normalizeRows();
		allGenes = Exp.getRowNames();
		ScoreMap = new HashMap<String, Double>();
		sample = Exp.getColNames();
	}
	
	public void perform(){
		calculateZvalues();
	}
	
	private void calculateZvalues() {
		ScoreMat = new MyMat( new ArrayList<String>(seedGeneSets.keySet()),sample);	
		int k=0;
		int n = seedGeneSets.size();
		for(String geneSet: seedGeneSets.keySet()){
			k++;
			List <String> genes = seedGeneSets.get(geneSet);
			List <Double> exp  = getMeanExp(genes);
			List <List <Double>> nullMeanAndSd  = getNullMeanAndSd(genes.size());
			List <Double> vec  = new ArrayList<Double>();
			for(int i = 0; i < exp.size(); i++){
				double nullMean = nullMeanAndSd.get(i).get(0);
				double nullSd = nullMeanAndSd.get(i).get(1);
				double z = (exp.get(i)-nullMean)/nullSd; 
				ScoreMat.set(geneSet,sample.get(i), z);
				if(useZscore){
					vec.add(z);
				}else{
					try {
						double p = 1-ND.cumulativeProbability(Math.abs(z));
						p = p==1?0:-Math.log10(p);
						p = (Double.isInfinite(p) || Double.isNaN(p)) ? ceilingScore : p;
						p = (z > 0)?p:-p;
						vec.add(p);
					} catch (MathException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
						vec.add(0.0);
					}
				}
			}
			double Score  = MyFunc.sd(vec);
			ScoreMap.put(geneSet, Score);
			System.err.println(geneSet +  "(" + k + "/"  + n + "):"  + "\t" + Score + "\t" + vec);	
		}
	}
	
	private List <Double> getMeanExp(List <String> genes){
			List <Double> exp = new ArrayList<Double>(sample.size());
			for(int i = 0;i<sample.size(); i++){
				exp.add(0.0);
			}
			for(String gene: genes){
				for(int  i=0; i<sample.size(); i++){
					exp.set(i, exp.get(i)+ Exp.get(gene, sample.get(i)));					
				}	
			}
			for(int  i=0; i<sample.size(); i++){
				exp.set(i, exp.get(i)/genes.size());					
			}	
			return exp;
						
	}

 	private List <List <Double>> getNullMeanAndSd(int size){
 		MyMat nullDist = new MyMat(itrForZvalueCalculation,sample.size());
 		
 		for(int i=0; i<itrForZvalueCalculation; i++){
 			List <Double> tmp = getMeanExp(MyFunc.sample(allGenes,size));
 			for(int j=0; j<sample.size(); j++){
 				nullDist.set(i, j, tmp.get(j));
 			}	 			
 		}
 		
 		List <List <Double>> nullMeanAndSd = new  ArrayList<List<Double>>();
 		for(int j=0; j<sample.size(); j++){
 			List <Double> nullDistForEachSample = nullDist.getCol(j);
 			List <Double> tmp = new ArrayList<Double>();
 			tmp.add(MyFunc.mean(nullDistForEachSample));
 			tmp.add(MyFunc.sd(nullDistForEachSample));
 			nullMeanAndSd.add(tmp);
 		}
 		
 		return nullMeanAndSd;
 	}
 	
 	public Map <String, Double> getScores(){
		return ScoreMap;
	}
	
 	
 	public static void main(String [] args) throws Exception{
 		Options options = new Options();
		options.addOption("o", "outfile", true, "outFile");
		options.addOption("m", "mingeneset", true, "minGeneSetSize");
		options.addOption("M", "maxgeneset", true, "maxGeneSetSize");
		options.addOption("i", "itr", true, "# of null statistics to be calculated");
		options.addOption("s", "sample", false, "print p-value for each sample");
		options.addOption("z", "zscore", false, "use z scores instead of p values");
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
		 OneSampleOverExpression2 GSA = new OneSampleOverExpression2(new MyMat(argList.get(0)));	
		GSA.setGeneSets(MyFunc.readGeneSetFromGmtFile(argList.get(1)));

		if(commandLine.hasOption("m")){
			GSA.setMinGeneSetSize(Integer.valueOf(commandLine.getOptionValue("m")));
		}
		if(commandLine.hasOption("M")){
			GSA.setMaxGeneSetSize(Integer.valueOf(commandLine.getOptionValue("M")));
		}
		if(commandLine.hasOption("z")){
			GSA.useZscore = true;
		}
		if(commandLine.hasOption("s")){
			GSA.printScoresForEachSample = true;
		}
		
		GSA.perform();
		Map <String, Double> Zvalues = GSA.getScores();
		if(commandLine.hasOption("o")){
			PrintWriter os = new PrintWriter(new FileWriter(commandLine.getOptionValue("o")));
			if(GSA.printScoresForEachSample){
				os.print(GSA.ScoreMat);	
			}else{
			for(String s: MyFunc.sortKeysByDescendingOrderOfValues(Zvalues)){
				 os.println(s + "\t" + Zvalues.get(s));				 
			 }
			}
			 os.close();
		}else{
			if(GSA.printScoresForEachSample){
				 System.out.print(GSA.ScoreMat);	
			}else{
			 for(String s: MyFunc.sortKeysByDescendingOrderOfValues(Zvalues)){
				 System.out.println(s + "\t" + Zvalues.get(s));				 
			 }
			}
			 System.out.close();
		}
		
	}
	
	
	
	
	
}
