package eem;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import utility.MyFunc;

public class AbstractGeneSetAnalysis implements GeneSetAnalysis  {
	List<String> allGenes;
	int maxSeedGeneSize = 2000;
	int minSeedGeneSize = 10 ;
	Map<String, List<String>> seedGeneSets = new HashMap<String, List<String>>();
	public void setGeneSets(Map<String, List<String>> geneSets){
		for(Map.Entry<String, List<String>> e: geneSets.entrySet() ){
			List <String> tmp = MyFunc.isect(e.getValue(), allGenes);
			if(tmp.size() < minSeedGeneSize || tmp.size() > maxSeedGeneSize){
				System.err.println( e.getKey() +  ": seed geneset size ("+ tmp.size() + ") is out of range!");
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
	
	public void perform(){
		throw new UnsupportedOperationException();
	}
	
	
	
	
}
