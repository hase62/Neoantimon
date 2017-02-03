package eem;

import java.io.*;
import java.util.*;

import utility.*;

public abstract class AbstractEEMsearch extends AbstractGeneSetAnalysis  implements EEMsearch  {
	MyMat originalExp;
	int itrForPvalue2Calculation;
	Double Pvalue1Cutoff;
	StopWatch stopWatch;
	
	
	Map<String, EEM> eems;	
	List<String> seeds;
	List<String> candidates;
	String errLog; 
	String timeLog;
	Map <Integer, List <Integer>> nullDistrbutionForRecycle;
	
	
	
	public AbstractEEMsearch(){
		maxSeedGeneSize = 2000;
		minSeedGeneSize = 10;
		itrForPvalue2Calculation = 500;
		Pvalue1Cutoff = 1.0;
		stopWatch = new StopWatch();
		seedGeneSets = new HashMap<String, List<String>>();
		eems = new HashMap<String, EEM>();
		seeds = new ArrayList<String>();
		candidates  = new ArrayList<String>();
		errLog = ""; 
		timeLog = "";
		nullDistrbutionForRecycle = null;
	}
	
	
	public AbstractEEMsearch(AbstractEEMsearch A){
		originalExp = A.originalExp;
		allGenes = A.allGenes;
		maxSeedGeneSize = A.maxSeedGeneSize;
		minSeedGeneSize = A.minSeedGeneSize;
		itrForPvalue2Calculation = A.itrForPvalue2Calculation;
		Pvalue1Cutoff = A.Pvalue1Cutoff;
		stopWatch = new StopWatch();
		seedGeneSets = new HashMap<String, List<String>>(A.seedGeneSets);
		eems = new HashMap<String, EEM>();
		seeds = new ArrayList<String>();
		candidates  = new ArrayList<String>(); 
		errLog = ""; 
		timeLog = "";
		if(A.nullDistrbutionForRecycle == null){
			nullDistrbutionForRecycle = null;
		}else{
			nullDistrbutionForRecycle = new HashMap<Integer, List<Integer>>();
		}
	}
	
	
	public void setGeneSets(Map<String, List<String>> geneSets){
		for(Map.Entry<String, List<String>> e: geneSets.entrySet() ){
			List <String> tmp = MyFunc.isect(e.getValue(), allGenes);
			if(tmp.size() < minSeedGeneSize || tmp.size() > maxSeedGeneSize){
				System.err.println( e.getKey() +  ": seed geneset size ("+ tmp.size() + ") is out of range!");
				errLog += e.getKey() +  ": seed geneset size ("+ tmp.size() + ") is out of range!\n";
				continue; 
			}
			seedGeneSets.put(e.getKey(), tmp);	
		}
	}

	public void setPvalue1Cutoff(double d){
		Pvalue1Cutoff = d;
	}
	
	public void suppressPvalue1Cutoff(){
		Pvalue1Cutoff = null;
	}

	
	
	public void setItrForPvalue2Calculation(int i){
		itrForPvalue2Calculation= i;
	}
	
	public void recycleNullDistribution(){
		nullDistrbutionForRecycle = new HashMap<Integer, List<Integer>>();
	}
	
	protected  void findModuleGenes() {
		System.err.println("Findng module gene...");
		int i = 0;
		int n = eems.size();
		List <String> new_seeds = new ArrayList<String>();
		for(String s : seeds){
			try{
			i++;
			eems.get(s).findModuleGenes();
			System.err.println(s + "(" + i + "/"  + n + "): succeed! " + eems.get(s).getModuleGenes().size() + "/" + eems.get(s).getSeedGenes().size());
			new_seeds.add(s);
			}
			catch (Exception err) {
				System.err.println(s + "(" + i + "/"  + n + "): failed!");
				errLog += s +  ": unable to find module genes!\n"; 
				continue;
			}
		}
		seeds = new_seeds;
		if(seeds.isEmpty()){
			throw new MyException("findModuleGenes: Any module has passed!");
		}
		
	}
	
	protected  void calculatePvalue1(){
		if(Pvalue1Cutoff == null){
			return;
		}else{
			System.err.println("Calculating approximate P values to filter out non-significant genesets...");
			int i = 0;
			int n = seeds.size();
			List <String> new_seeds = new ArrayList<String>();
			for(String s: seeds){
				try{
					i++;
					eems.get(s).calculatePvalue1();
					System.err.println(s + "(" + i + "/"  + n + "): " + "succeed! " + (eems.get(s)).getPvalue1());
					new_seeds.add(s);
					}
					catch (Exception err) {
					System.err.println(s + "(" + i + "/"  + n + "): failed!");
					errLog +=  s +  ": unable to calculate a P value  based on hypergeometric distribution!\n";
					continue;
				}
			}
			seeds = new_seeds;
			if(seeds.isEmpty()){
				throw new MyException(" calculatePvalue1: Any module has passed!");
			}
		}
	}
	
	protected  void findCandidates(){
		if(Pvalue1Cutoff == null){
			candidates = new ArrayList<String>(seeds);
		}else{
			for(String s: seeds){
				if(eems.get(s).getPvalue1() >= Pvalue1Cutoff){
					candidates.add(s);
				}
			}
			if(candidates.isEmpty()){
				throw new MyException("findCandidates: Any module has passed!");
			}
		}
	}
	
	
	protected  void calculatePvalue2(){
		System.err.println("Calculating accurate P values...");
		int i = 0;
		int n = candidates.size();
		List <String> new_candidates = new ArrayList<String>();
		for(String s: candidates){
			try{
				i++;
				eems.get(s).calculatePvalue2();
				System.err.println(s + "(" + i + "/"  + n + "): " + "succeed! " + (eems.get(s)).getPvalue2());
				new_candidates.add(s);
			}
			catch (Exception err) {
				err.printStackTrace();
				System.err.println(s + "(" + i + "/"  + n + "): failed!");
				errLog +=  s +  ": unable to calculate a P value based on arandomization test!\n"; 
				continue;
			}
		}
		candidates = new_candidates;
		if(candidates.isEmpty()){
			throw new MyException("calculatePvalue2: Any module has passed!");
		}
		
	}
	
	
	public void printResults(String outfile) throws IOException{
		Map<String, Double> p = new HashMap<String, Double>();
		for(Map.Entry<String, EEM> e: eems.entrySet()){
			 p.put(e.getKey(), (e.getValue()).getPvalue());
		}
		List<String> tmp = MyFunc.sortKeysByDescendingOrderOfValues(p);
		PrintWriter os = new PrintWriter(new FileWriter(outfile));
		int i;
		for(i=0;i<tmp.size();i++){
			String id = tmp.get(i);
			EEM e = eems.get(id);
			os.println(id + "\t" + e.getModuleGenes().size() + "/" + e.getSeedGenes().size() + "\t" + e.getPvalue() + "\t" + MyFunc.join("\t", e.getModuleGenes()));
		}
		os.close();
		
	}
	public void printLog(String outfile) throws IOException{
		PrintWriter os = new PrintWriter(new FileWriter(outfile));
		os.println(getLog());
		os.close();
	}
	protected void setEEM(){
		throw new UnsupportedOperationException();
	}
	public ExpressionModuleSet getExpressionModuleSet(){
		ExpressionModuleSet tmp = new ExpressionModuleSet();
		for(Map.Entry<String, EEM> e: eems.entrySet()){
			if(e.getValue().getPvalue() != null){
				e.getValue().cutParent();
				tmp.add(new ExpressionModule(e.getKey(), e.getValue(), originalExp)); 
			}
		}
		return tmp;
	}
	public String getLog(){
		return timeLog + "\n" + this + "\n" + errLog;
		
	}
	protected  void writeTimeLog(){
		timeLog = "started at " +  stopWatch.getStartDate() + "\n" +
				"finished at " + stopWatch.getStopDate() + "\n"+
				"It took " + stopWatch + "\n";
	}
	
	public void perform(){
		try{
			setEEM();
			findModuleGenes();
			calculatePvalue1();
			findCandidates();
			calculatePvalue2();
			
		} catch (Exception e) {
			System.err.println(e.getMessage());
			//e.printStackTrace();
		}
		finally{
			stopWatch.stop();
			writeTimeLog();
		}
	}
	

	public Map <String, Double>getPvalues(){
		Map<String, Double> P = new HashMap<String, Double>();
		for(Map.Entry<String, EEM> e: eems.entrySet()){
			if(e.getValue().getPvalue() != null){
				P.put(e.getKey(), (e.getValue()).getPvalue());
			}
		}
		return P;
	}
	
}
