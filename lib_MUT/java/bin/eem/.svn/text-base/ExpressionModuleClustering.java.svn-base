package eem;

import java.util.*;
import utility.*;

public class ExpressionModuleClustering {
	 
	private ExpressionModuleSet expModSet;
	private Double PvalueCutoff; 
	private int bgGeneNumber = 20000;
	
	public ExpressionModuleClustering(ExpressionModuleSet expModSet){
		this.expModSet = expModSet;
		PvalueCutoff = 0.001;
	}

	public void clear(){
		ExpressionModuleSet tmp = new  ExpressionModuleSet(expModSet);
		expModSet.clear();
		for(ExpressionModule e: tmp.asList()){
			if(e.expressionModuleCluster != null){
				expModSet.addAll(e.expressionModuleCluster);
			}else{
				expModSet.add(e);
			}
		}
		for(ExpressionModule e: expModSet.asList()){
			e.expressionModuleCluster = null;
		}
	}
		
	public void setPvalueCutoff(double d){
		if(d <= 0){
			return;	
		}
		if(d < 1){
			PvalueCutoff = d;
		}else{
			PvalueCutoff = Math.pow(10, -d);
		}
	}
	
	public void setBgGeneNumber(int n){
		bgGeneNumber = n;
	}
	
	private double calculatePvalue(ExpressionModule e1, ExpressionModule e2){
		List <String> module1 = e1.getModuleGenes();
		List <String> module2 = e2.getModuleGenes();
		int isect = MyFunc.isect(module1, module2).size();
		double P = (isect == 0)?1: MyFunc.calculatePvalueForSetOverlap(bgGeneNumber, module1.size(), module2.size(),isect );
		return P;
	}
	
	public void perform(){
		ExpressionModuleSet copy = new  ExpressionModuleSet(expModSet);
		expModSet.clear();
		while(!copy.isEmpty()){
			ExpressionModule topExpMod = copy.get(0);
			ExpressionModuleSet currentCluster = new ExpressionModuleSet();
			currentCluster.add(topExpMod);
			for(int i = 1, n = copy.size(); i < n; i++){
				if(calculatePvalue(topExpMod, copy.get(i)) <=  PvalueCutoff){
					currentCluster.add(copy.get(i));
				}
			}
			for(int i = 0, n = currentCluster.size(); i < n; i++){
				copy.remove(currentCluster.getIds().get(i));
			}
			for(ExpressionModule e: currentCluster.asList()){
				e.expressionModuleCluster = currentCluster;
			}
			expModSet.add(topExpMod);
		}
	}
	
	
	public String toString(){
		List <String> tmp = new ArrayList<String>();
		for(ExpressionModule e: expModSet.asList()){
			if(e.expressionModuleCluster != null){
				tmp.add(MyFunc.join("\t",e.expressionModuleCluster.getIds()));
			}
		}
		if(tmp.isEmpty()){
			return "Remained to be clustered!";
		}else{
			return MyFunc.join("\n", tmp);
		}
	}
	
	public static String toString(ExpressionModuleSet expModSet){
		List <String> tmp = new ArrayList<String>();
		for(ExpressionModule e: expModSet.asList()){
			if(e.expressionModuleCluster != null){
				tmp.add(MyFunc.join("\t",e.expressionModuleCluster.getIds()));
			}
		}
		if(tmp.isEmpty()){
			return "Remained to be clustered!";
		}else{
			return MyFunc.join("\n", tmp);
		}
	}
}
