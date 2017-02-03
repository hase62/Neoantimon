package eem;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

import utility.*;

public abstract  class AbstractParameterOptimizingEEMsearch implements EEMsearch, GeneSetAnalysis{
	String timeLog;
	StopWatch stopWatch;
	List <EEMsearch> EEMsearches;
	ExpressionModuleSet expModSet;
		
	EEMsearch EEMsearchTemplate;
	
	public AbstractParameterOptimizingEEMsearch() {
		stopWatch = new StopWatch();
		EEMsearches  = new  ArrayList<EEMsearch>();
		expModSet = new ExpressionModuleSet();
		timeLog = "";
	}
	
	public void setPvalue1Cutoff(double d){
	    EEMsearchTemplate.setPvalue1Cutoff(d);
	}

	public void suppressPvalue1Cutoff(){
		EEMsearchTemplate.suppressPvalue1Cutoff();
	}
	
	public void setMaxGeneSetSize(int i){
		EEMsearchTemplate.setMaxGeneSetSize(i);
	}
	public void setMinGeneSetSize(int i){
		EEMsearchTemplate.setMinGeneSetSize(i);
	}	
	public void setItrForPvalue2Calculation(int i){
		EEMsearchTemplate.setItrForPvalue2Calculation(i);
	}

	public void setGeneSets(Map<String, List<String>> geneSets){
		EEMsearchTemplate.setGeneSets(geneSets);
	}
	public void recycleNullDistribution(){
		EEMsearchTemplate.recycleNullDistribution();
	}
	
	protected void initializeEEMsearches(){
		throw  new UnsupportedOperationException();
	}
	
	protected void unifyResults(){
		for(EEMsearch e: EEMsearches){
			expModSet.addAll(e.getExpressionModuleSet());
		}
		
	}
	public void perform(){
		initializeEEMsearches();
		for(EEMsearch e: EEMsearches){
			e.perform();
			System.err.println("");
		}
		unifyResults();
		stopWatch.stop();
		writeTimeLog();
	}
	public Map <String, Double>getPvalues(){
		Map<String, Double> P = new HashMap<String, Double>();
		for(String  s: expModSet.getIds()){
			P.put(s, expModSet.get(s).getPvalue());
		}
		return P;
	}
	public String getLog(){
		List <String> s = new ArrayList<String>();
		s.add( timeLog + "\n" + this + "\n");
		for(EEMsearch e: EEMsearches){
			s.add(e.getLog());
		}
		return MyFunc.join("***************************************\n", s);
	}
	
	private  void writeTimeLog(){
		timeLog = "started at " +  stopWatch.getStartDate() + "\n" +
				"finished at " + stopWatch.getStopDate() + "\n"+
				"It took " + stopWatch + "\n";
	}
	
	public void printLog(String outfile) throws IOException{
		PrintWriter os = new PrintWriter(new FileWriter(outfile));
		os.println("started at " +  stopWatch.getStartDate());
		os.println("finished at " + stopWatch.getStopDate());
		os.println("It took " + stopWatch);
		os.println("");
		os.println(getLog());
		os.flush();
		os.close();
	}
	
	public void printResults(String outfile) throws IOException{
		Map<String, Double> p = new HashMap<String, Double>();
		for(String s: expModSet.getIds()){
			p.put(s, expModSet.get(s).getPvalue() );
		}
		List<String> tmp = MyFunc.sortKeysByDescendingOrderOfValues(p);
		PrintWriter os = new PrintWriter(new FileWriter(outfile));
		int i;
		for(i=0;i<tmp.size();i++){
			String id = tmp.get(i);
			EEM e = expModSet.get(id).getEEM();
			os.println(id + "\t" + e.getModuleGenes().size() + "/" + e.getSeedGenes().size() + "\t" + e.getPvalue() + "\t" + MyFunc.join("\t", e.getModuleGenes()));
		}
		
	}
	
	public ExpressionModuleSet getExpressionModuleSet(){
		return expModSet;
	}
	
	
}
