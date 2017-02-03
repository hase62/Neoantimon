package eem;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;


import utility.MyFunc;

public class SimulationDataGeneratorUsingBiclusterCompositeModule extends SimulationDataGeneratorUsingBicluster {

	private int moduleNumberInOneGeneSet = 2;

	public void setModuleNumberInOneGeneSet(int k){
		moduleNumberInOneGeneSet = k;
	}
	
	public void simulateGeneset(){
		Random rnd = new Random();
		int moduleGeneSubsetSize = (int) Math.round(geneSetSize*moduleGeneSubsetRate);
		for(int i = 0; i < positiveGeneSetNumber; i++){
			List<String> tmp = new ArrayList<String>();
			Set <Integer> seen = new HashSet<Integer>();
			for(int j = 0; j < moduleNumberInOneGeneSet; j++){
				int k = rnd.nextInt(moduleNumber);
				if(seen.contains(k)){
					j--;
					continue;
				}else{
					tmp.addAll(new ArrayList<String>(moduleGeneList.get(k).subList(0,moduleGeneSubsetSize/moduleNumberInOneGeneSet)));
					tmp = MyFunc.uniq(tmp);
					seen.add(k);
				}
			}
			tmp.addAll(MyFunc.sample(MyFunc.diff(geneName, tmp), geneSetSize-tmp.size()));
			geneSet.put("positive" + (i+1),tmp);
		}
		for(int i = 0; i < negativeGeneSetNumber; i++){
			geneSet.put("negative"+(i+1), MyFunc.sample(geneName, geneSetSize));
		}
		
		
	}
	
	
	
	
	
	
	
	
}
