package bct;

import java.util.*;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

import sun.reflect.Reflection;
import utility.*;

public class BCTcompetitionTest {
	protected BCT T;
	protected Map<String, Integer> gene2index;
	protected List <String> genes;
	
	protected Map<String, List<String>> geneset;
	
	protected BCTsampler Bsamp;
	protected int NforBurnIn = 1000;
	protected int N = 10000;
	
	protected Map<String, Double> statistics;
	protected Map<String, Double> pvalue;
	
	boolean larger = false;
	
	public BCTcompetitionTest(MyMat M, Map<String, List<String>> geneset){
		setBCT(M);
		setGeneset(geneset);
	}
	
	void setBCT(MyMat M){
		T = new BCT(M);
		genes = new ArrayList<String>(M.getRowNames());
		gene2index = new HashMap<String,Integer>();
		for(int i=0;i<genes.size();i++){
			gene2index.put(genes.get(i),i);
		 }
		Bsamp = new BCTsampler(T);
	}
	
	void setGeneset(Map<String, List<String>> geneset){
		this.geneset = new HashMap <String, List<String>>();
		for(String s: geneset.keySet()){
			List <String> tmp = new ArrayList<String>();	
			for(String t: geneset.get(s)){
				if(gene2index.containsKey(t)){
					tmp.add(t);
				}
			}
			if(tmp.size() > 1){
				this.geneset.put(s, tmp);
			}
		}
	}
	
	protected double  getStatistics(List <String> genes, BCT T){
		double S = 0;
		double n = genes.size();
		for(String g1: genes){
			int i1 = gene2index.get(g1);
			for(String g2: genes){
				int i2 = gene2index.get(g2);
				if(i1 >= i2){
					continue;
				}
				double s = 0;
				for(int j =0; j < T.colSize(); j++){
					if(T.is1(i1, j) & T.is1(i2, j)){
						s++;
					}
				}
				S += Math.pow(s, 2);
			}
		}
		S /= n*(n-1);
		return S;
	}
	
	protected Map <String, Double> getStatistics(Map <String, List<String>> geneset, BCT T){
		 Map <String, Double> statistics = new HashMap <String, Double> (); 
		for(String s: geneset.keySet()){
			statistics.put(s, getStatistics(geneset.get(s), T));
		}
		return statistics;
	}
	
	protected void getStatistics(){
		statistics = getStatistics(geneset, T); 
	}
	
	protected void getPvalue(){
		pvalue = new HashMap<String, Double>();
		for(String s: geneset.keySet()){
			pvalue.put(s, 0.0);
		}
		for(int c =0; c < NforBurnIn; c++){
			Bsamp.mcmc.next();
		}
		for(int c = 0; c < N; c++){
			Bsamp.mcmc.next();
			Map <String, Double> nullStatistics = getStatistics(geneset,Bsamp.getCurrentBCT());
			if(larger){
				for(String s: geneset.keySet()){
					if(statistics.get(s) <= nullStatistics.get(s)){
						pvalue.put(s, pvalue.get(s)+1);
					}
				}
			}else{
				for(String s: geneset.keySet()){
					if(statistics.get(s) >= nullStatistics.get(s)){
						pvalue.put(s, pvalue.get(s)+1);
					}
				}
			}
		}
		for(String s: geneset.keySet()){
			if(pvalue.get(s)==0.0){
				pvalue.put(s, 1.0);
			}
			pvalue.put(s, pvalue.get(s)/N);
		}
	}
	
	protected void printResult(){
		for(String s:MyFunc.sortKeysByAscendingOrderOfValues(pvalue)){
			System.out.println(s + "\t" + pvalue.get(s) + "\t" +  statistics.get(s) + "\t" +  geneset.get(s) );
		}
	}
	
	public void perform(){
		getStatistics();
		getPvalue();
		printResult();
	}
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine = null;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] input.tab geneset.gmt", options);
			System.exit(1);
		}
		List <String> argList = commandLine.getArgList();
		if(!(argList.size() == 2)){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] input.tab geneset.gmt", options);
			System.exit(1);
		}
		BCTcompetitionTest B = new BCTcompetitionTest(new MyMat(argList.get(0)), MyFunc.readGeneSetFromGmtFile(argList.get(1)));
		B.perform();
	}
	
}

