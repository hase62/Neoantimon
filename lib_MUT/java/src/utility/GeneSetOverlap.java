package utility;
import java.util.*;
import java.io.*;


public class GeneSetOverlap {
	private Map <String, List <String>> geneSetMap1;
	private Map <String, List <String>> geneSetMap2;
	
	private List <String> bgGenes;
	
	private int bgGeneNumber = 20000;
	
	private MyMat Pvalue;
	private MyMat Qvalue;
	
	
	public GeneSetOverlap(Map <String, List <String> > geneSetMap1, Map <String, List <String> > geneSetMap2) {
		setGeneSetMap1(geneSetMap1);
		setGeneSetMap2(geneSetMap2);
		bgGenes = new ArrayList<String>();
	}	
	
	
	public void setBgGenes(List <String> Gene){
		bgGenes = new ArrayList<String>(Gene);
		bgGeneNumber = bgGenes.size();	
		if(!geneSetMap1.isEmpty()){
			Map <String, List<String>> tmp = new HashMap<String, List<String>>();	
			for(Map.Entry<String, List<String>> e: geneSetMap1.entrySet()){
				List <String> tmp2 = MyFunc.isect(bgGenes, e.getValue());
				if(!tmp2.isEmpty()){
					tmp.put(e.getKey(),tmp2);
				}			
			}
			geneSetMap1 = tmp;
		}
		if(!geneSetMap2.isEmpty()){
			Map <String, List<String>> tmp = new HashMap<String, List<String>>();	
			for(Map.Entry<String, List<String>> e: geneSetMap2.entrySet()){
				List <String> tmp2 = MyFunc.isect(bgGenes, e.getValue());
				if(!tmp2.isEmpty()){
					tmp.put(e.getKey(),tmp2);
				}			
			}
			geneSetMap2 = tmp;
		}
		
	}
	
	public void setBgGeneNumber(int n){
			bgGeneNumber = n;
	}
	
	public void setGeneSetMap1(Map <String, List <String> > geneSetMap){
		if(bgGenes.isEmpty()){
			geneSetMap1 = geneSetMap;
		}else{
			Map <String, List<String>> tmp = new HashMap<String, List<String>>();	
			for(Map.Entry<String, List<String>> e: geneSetMap.entrySet()){
				List <String> tmp2 = MyFunc.isect(bgGenes, e.getValue());
				if(!tmp2.isEmpty()){
					tmp.put(e.getKey(),tmp2);
				}			
			}
			geneSetMap1 = tmp;
		}	
	}
	
	public void setGeneSetMap2(Map <String, List <String> > geneSetMap){
		if(bgGenes.isEmpty()){
			geneSetMap2 = geneSetMap;
		}else{
			Map <String, List<String>> tmp = new HashMap<String, List<String>>();	
			for(Map.Entry<String, List<String>> e: geneSetMap.entrySet()){
				List <String> tmp2 = MyFunc.isect(bgGenes, e.getValue());
				if(!tmp2.isEmpty()){
					tmp.put(e.getKey(),tmp2);
				}			
			}
			geneSetMap2 = tmp;
		}	
	}
	
	public void calculatePvalue(){
		Pvalue = new MyMat(new ArrayList<String>(geneSetMap1.keySet()), new ArrayList<String>(geneSetMap2.keySet()));
		for(Map.Entry<String, List<String>> e1: geneSetMap1.entrySet()){
			for(Map.Entry<String, List<String>> e2: geneSetMap2.entrySet()){
				int isect = (MyFunc.isect(e1.getValue(), e2.getValue())).size();
				double P = (isect == 0) ? 1 : MyFunc.calculatePvalueForSetOverlap(bgGeneNumber, e1.getValue().size(), e2.getValue().size(),isect );
				Pvalue.set(e1.getKey(), e2.getKey(), P);
			}
		}
	}

	public void calculateQvalue(){
		Qvalue = new MyMat(new ArrayList<String>(geneSetMap1.keySet()), new ArrayList<String>(geneSetMap2.keySet()));
		Map <String, Double> PvalueMap = Pvalue.asMap();
		Map <String, Double> QvalueMap = MyFunc.calculateStoreyQvalue(PvalueMap);
		for(Map.Entry<String, Double> e: QvalueMap.entrySet()){
			List <String> tmp = MyFunc.split("\t", e.getKey());
			Qvalue.set(tmp.get(0), tmp.get(1), e.getValue());
		}
	}
	
	
	public String toString(){
		StringBuffer S = new StringBuffer("\tgene set1\tgene set2\tP value\tQ value");
		Map <String, Double> PvalueMap = Pvalue.asMap();
		Map <String, Double> QvalueMap = Qvalue.asMap();
		List<String> keys =  MyFunc.sortKeysByAscendingOrderOfValues(PvalueMap);
		for(String s: keys){
			List <String>  tmp = new ArrayList<String>();
			tmp.add(s);
			tmp.add(Double.toString(PvalueMap.get(s)));
			tmp.add(Double.toString(QvalueMap.get(s)));
			S.append(MyFunc.join("\t", tmp) + "\n");
		}
		return S.toString();
	}
	
	
	
}
