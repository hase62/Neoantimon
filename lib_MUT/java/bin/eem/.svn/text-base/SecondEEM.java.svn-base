package eem;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import utility.MyException;
import utility.MyFunc;

public class SecondEEM extends CoherenceBasedEEM{

	protected List<String> seedGenes2;
	protected List<String> moduleGenes2;
	
	
	SecondEEM(CoherenceBasedEEMsearch parent, List<String> gene) {
		super(parent, gene);
		
	}
	
	protected  void findSecondCenter(){
		int i,j;	
		moduleGenes2 = new ArrayList<String>();
		seedGenes2 = MyFunc.diff(seedGenes, moduleGenes);
		centerGene = null; 
		for(i =0; i < seedGenes2.size(); i++){
			List <String> tmp = new ArrayList<String>();	
			for(j =0; j < seedGenes2.size(); j++){
				if(i == j){
					tmp.add(seedGenes2.get(j));
					continue;
				}
				double tmp2 = parent.distfunc.get(seedGenes2.get(i), seedGenes2.get(j));
				if(tmp2 < absoluteRadius){
					tmp.add(seedGenes2.get(j));
				}
			}
			if(moduleGenes2.isEmpty() || tmp.size() > moduleGenes2.size() ){
				moduleGenes2  = tmp;
				center = parent.Exp.getRow(seedGenes2.get(i));
				centerGene = seedGenes2.get(i);
			}
		}
	}
	protected  void refineSecondCenter() throws MyException{
		int i,j;	
		if(moduleGenes2.size() < 3){
			throw  new  MyException("unnable to find an enough number of module genes!");
		}
		List <String> coreGenes ;	
		if(coreGeneSize < moduleGenes2.size()){
			Map <String, Double> dist2center  = new HashMap<String, Double>();
			for(i =0; i < moduleGenes2.size(); i++){
			dist2center.put(moduleGenes2.get(i), parent.distfunc.get(moduleGenes2.get(i), centerGene));
			}
		coreGenes = MyFunc.sortKeysByAscendingOrderOfValues(dist2center);
		coreGenes = coreGenes.subList(0, coreGeneSize-1);	
		}else{
			coreGenes = moduleGenes2;
		}
		List<List<String>> CoreTriplet = MyFunc.getAllCombination(coreGenes, 3);
		for(i=0; i<CoreTriplet.size();i++){
			List <String> tmp = new ArrayList<String>();
			List<Double> tripletCenter = parent.Exp.getRowMeans(CoreTriplet.get(i));
			for(j =0; j < seedGenes2.size(); j++){
				double tmp2 = parent.distfunc.get(tripletCenter, parent.Exp.getRow(seedGenes2.get(j)));
				if(tmp2 < absoluteRadius){
					tmp.add(seedGenes2.get(j));
				}
			}
			if(tmp.size() > moduleGenes2.size() ){
				moduleGenes2  = tmp;
				center = tripletCenter;
			}
		}
	}
	public void findModuleGenes() throws MyException {
		findCenter();
		//try {
		//	refineCenter();
		//}
		//catch (Exception err) {
		//}
		findSecondCenter();
		try{
			refineSecondCenter();
		}
		catch(Exception err) {
		}
		moduleGenes = new ArrayList<String>();
		for(int j =0; j < seedGenes.size(); j++){
			double d = parent.distfunc.get(center, parent.Exp.getRow(seedGenes.get(j)));
			if(d < absoluteRadius){
				moduleGenes.add(seedGenes.get(j));
			}
		}
	}
	

}
