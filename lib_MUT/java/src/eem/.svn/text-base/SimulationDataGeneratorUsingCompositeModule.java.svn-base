package eem;


import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

import org.apache.commons.math.random.*;

import utility.*;

public class SimulationDataGeneratorUsingCompositeModule extends SimulationDataGenerator {
	
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
			for(int j = 0; j < moduleNumberInOneGeneSet;j++){
				int k = rnd.nextInt(moduleNumber);
				if(seen.contains(k)){
					j--;
					continue;
				}else{
					tmp.addAll(geneName.subList(k*moduleSize, k*moduleSize + moduleGeneSubsetSize/moduleNumberInOneGeneSet));
					seen.add(k);
				}
			}
			tmp.addAll(MyFunc.sample(MyFunc.diff(geneName, tmp), geneSetSize - tmp.size()));
			geneSet.put("positive" + (i+1),tmp);
		}
		for(int i = 0; i < negativeGeneSetNumber; i++){
			geneSet.put("negative"+(i+1), MyFunc.sample(geneName, geneSetSize));
		}
		
	}
	
	
	
	
}
