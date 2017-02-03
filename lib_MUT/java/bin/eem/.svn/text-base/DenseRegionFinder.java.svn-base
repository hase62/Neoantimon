package eem;


import java.util.*;

import utility.*;

public class DenseRegionFinder {
	private int numberOfDenseRegions = 10;
	private List <CoherenceBasedEEM> EEMs;
	private List <String> targetGenes; 
	private CoherenceBasedEEMsearch parent;
	private List <Integer> Density;	// the number of genes in dense regions
	//private List <String> centerGenes;
	private List <String> allGenes;
	
	public DenseRegionFinder(CoherenceBasedEEMsearch parent){
		EEMs = new ArrayList<CoherenceBasedEEM>();
		targetGenes = new ArrayList<String>(parent.allGenes);	
		Density = new ArrayList<Integer>();
		this.parent = parent;
		allGenes = parent.allGenes;
	}
		
	public void find(){
		System.err.println("Finding dense regions in expression space...");
		for(int i = 0; i < numberOfDenseRegions && !targetGenes.isEmpty(); i++){
			CoherenceBasedEEM eem  = new CoherenceBasedEEM(parent, targetGenes);
			eem.findModuleGenes();
			targetGenes = MyFunc.diff(targetGenes, eem.getModuleGenes());
			EEMs.add(eem);
			List <String> tmp = new ArrayList<String>();		
			Density.add(eem.getModuleGenes().size());
		}
		numberOfDenseRegions = Density.size();
	}

	public int getNumberOfDenseRegions(){
		return numberOfDenseRegions;	
	}
	
	
	public int getNearestDenseRegion(List <String> genes){
		
		int count[] = new int[numberOfDenseRegions]; 
		for(int i = 0; i < numberOfDenseRegions; i++){
			count[i] =  MyFunc.isect(EEMs.get(i).getModuleGenes(), genes).size(); 	
		}
		
		int maxIndex = 0;
		int maxCount = count[0];
		for(int i = 1; i < numberOfDenseRegions; i++){
			if(count[i] > maxCount){
				maxCount = count[i];
				maxIndex = i;
			}
		}
		return maxIndex;	
	}
	
	
	public Map<Integer, Double> getDensityDistribution( int geneSetSize, int numberOfRandomSampling){
		
		double count[]  = new double[numberOfDenseRegions]; 
		
		for(int i = 0; i <  numberOfRandomSampling; i++){
			count[getNearestDenseRegion(MyFunc.sample(parent.allGenes,geneSetSize))]++;
		}	
		for(int i = 0; i < numberOfDenseRegions; i++){
			count[i] /= numberOfRandomSampling;
		}
		Map<Integer, Double> Dist = new HashMap<Integer, Double>();
		
		for(int i = 0; i <   numberOfDenseRegions; i++){
			
			//System.err.println(Density.get(i) + "\t"  + count[i]);
			if(count[i] == 0.0){
				continue;
			}
			
			
			if(Dist.containsKey(Density.get(i))){
				Dist.put(Density.get(i), Dist.get(Density.get(i)) + count[i]);
			}else{
				Dist.put(Density.get(i), count[i]);
			}
		}
		return Dist;
	}
	
	
	
	
	
	
	
	
	
}
