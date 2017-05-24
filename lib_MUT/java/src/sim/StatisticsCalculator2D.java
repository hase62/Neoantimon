package sim;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;
import sim.Simulator.Cell;
import utility.MyFunc;
import utility.MyMat;

public class StatisticsCalculator2D {
	
	//parameter
	int mutationProfileMatrixSize = 1000;
	int mutationProfileMatrixSizeBySliceSampling = 30;
	int sliceSize = 21;
	double emptyPositionProportionCutoff = 0.1;
	double subpopulationProportionCutoff = 0.01;
	double populationCoverageForfounderMutation = 0.95;
	
	//input
	Simulator2D S;
	
	//statistics
	int time;
	int populationSize;
	
	double mutatedGeneCount;
	double mutatedDriverGeneCount;
	double founderMutationPropotion;
	double populationGrowthRate;
	Map<String, Double> subpopulationProportion;
	Map<String, Double> majorSubpopulationProportion;
	Double populationEntropy;
	
	MyMat mutationProfileMatrix;
	MyMat mutationSimilalityMatrix;
	MyMat mutationProfileMatrixBySliceSampling;
	
	//random number generator
	Random R = new Random();
	
	public StatisticsCalculator2D (Simulator2D S){
		this.S = S;
		time = S.time;
		populationSize = S.populationSize;
	}
	 
	void getSubpopulationProportion(){
		 subpopulationProportion = new HashMap<String, Double>();
		 List <Long>  genomes = S.getGenomesAsList();
		 for(Long g:genomes){
			 String k = S.genome2mutatedGenes(g).toString();
			 if(subpopulationProportion.containsKey(k)){
				 subpopulationProportion.put(k, subpopulationProportion.get(k)+1.0);
			}else{
				subpopulationProportion.put(k, 1.0);
				}
			}
			for(String k: subpopulationProportion.keySet()){
				subpopulationProportion.put(k, subpopulationProportion.get(k)/populationSize);
			}
	}
	 
	void getMutationProfileMatrixBySliceSampling(){
		mutationProfileMatrixBySliceSampling = new MyMat(S.genomeSize,mutationProfileMatrixSizeBySliceSampling);
		int sliceSizeHalf = sliceSize/2;
		List <String> colnames = new ArrayList <String>();
		for(int k = 0; k < mutationProfileMatrixSizeBySliceSampling; k++){
			int x=0, y=0;		  
			List <Long> genomesInSlice = new ArrayList<Long>();
			int l;
			L: for(l=0; l==100; l++){
				x = R.nextInt(S.maxX-S.minX+1) +S.minX;
				y = R.nextInt(S.maxY-S.minY+1) + S.minY;
				List <Integer> I = new ArrayList<Integer>();
				for(int i = -sliceSizeHalf; i <= sliceSizeHalf; i++){
					I.add(i);
				}	
				double emptyPositionCount = 0;
				for(Integer i: I){
					for(Integer j: I){
						int x2 = x+i;
						int y2 = y+j;
						if(Math.abs(x2) >  S.limXY | Math.abs(y2) >  S.limXY){
							genomesInSlice.clear();
							continue L;
						}
						if(S.isEmpty(x2, y2)){
							emptyPositionCount++;
							continue;
						}
						genomesInSlice.add(S.get(x2,y2));
					}
				}
				double emptyPositionProportion = (double)emptyPositionCount/Math.pow(2*sliceSizeHalf+1,2);
				if(emptyPositionProportion > emptyPositionProportionCutoff){
					genomesInSlice.clear();
					continue L;
				}
				break;
			}
			if(l==100){
				System.err.println("err: cannnot get MutationProfileMatrixBySliceSampling!");
				return;
			}
			for(Long g: genomesInSlice){
				Set<Integer> tmp = S.genome2mutatedGenes(g);
				for(int i = 0; i < S.genomeSize; i++){
					if(tmp.contains(i)){
						mutationProfileMatrixBySliceSampling.set(i,k, mutationProfileMatrixBySliceSampling.get(i,k)+1);
					}
				}	 
			 }
			for(int i = 0; i < S.genomeSize; i++){
				mutationProfileMatrixBySliceSampling.set(i,k, mutationProfileMatrixBySliceSampling.get(i,k)/genomesInSlice.size());
			 }
			String X, Y;
			if(x==0){
				X = "0";
			 }else if(x<0){
				 X = "n" + -x; 
			 }else{
				 X = "p" + x;
			 }
			if(y==0){
				Y = "0";
			}else if(y<0){
				Y = "n" + -y; 
			}else{
				Y = "p" + y;
			}
			colnames.add(X + "_" + Y);			 
		}
		mutationProfileMatrixBySliceSampling.setColNames(colnames);
	}
		 			
	static MyMat getSimilarityMatrix(MyMat m){
			 MyMat M = new MyMat(m.colSize(), m.colSize());
			 for(int i = 0; i <  m.colSize(); i++){
				 for(int j = 0; j <  i; j++){
					 double tmp = 0;
					 for(int k = 0; k <  m.rowSize(); k++){
						 //tmp += Math.pow(m.get(k, i)-m.get(k, j), 2); 
						 tmp += Math.abs(m.get(k, i)-m.get(k, j)); 
					}
					 tmp = 1 - tmp/m.rowSize();
					 M.set(i, j, tmp);
					 M.set(j, i, tmp);
				 }	
			 }
			 return M;
		}
		 
	static double getEigenValueEntropy(MyMat m) throws NotConvergedException{
			 DenseMatrix  M = m.getMTJDenseMatrix();
			 SVD svd = SVD.factorize(M);
			 double[] sv =  svd.getS();
			 double ic = shannonEntropy(sv);
			 return  ic;
		}
			
	static double shannonEntropy(double[] v){
			 double e = 0;
			 double s = 0;
			 List <Double> v2 = new ArrayList<Double>();
			 for(int i = 0, n = v.length; i < n; i++){
				 s  +=  v[i];
				 v2.add(v[i]);
			 }
			 for(Double d: v2){
				 if(d > 0){
					 double p = d/s;
					 e += p*Math.log(p);
				 }
			 }
			 return - e/Math.log(v.length);   
		}
		
	void getMutationProfileMatrix(){
		mutationProfileMatrix = new MyMat(S.genomeSize,mutationProfileMatrixSize);
		for(int k = 0; k < mutationProfileMatrixSize; k++){
			int x, y;
			while(true){
				x = R.nextInt(S.maxX-S.minX+1) +S.minX;
				y = R.nextInt(S.maxY-S.minY+1) + S.minY;
				if(!S.isEmpty(x, y)){
					break;
				}
			}
			Set<Integer> tmp = S.genome2mutatedGenes(S.get(x, y));
			for(int i = 0; i < S.genomeSize; i++){
				if(tmp.contains(i)){
					mutationProfileMatrix.set(i,k, mutationProfileMatrix.get(i,k)+1);
				}
			}	 
		}
	}
					 
	void getPopulationGrothRate(){
			 populationGrowthRate = 0;
			 int n = 0;
			 for(int x  = S.minX; x <= S.maxX ; x++){
				 for(int y  = S.minY; y <= S.maxY ; y++){
					 if(S.isEmpty(x,y)){
						 continue;
					 }
					 populationGrowthRate  += 1+S.getFitness(x,y)*S.growthRate;
					 if(S.isStem(x,y)){
						 populationGrowthRate  -= S.deathRate;
					 }else{
						 populationGrowthRate  -= S.deathRateForNonStem;
					}
					 n++;
				 }
			 }
			 populationGrowthRate /=n;
		 
		 }
		 		 
	void getMutatedGeneCount(){
				mutatedDriverGeneCount = 0;
				mutatedGeneCount = 0;
				Map<String, Double> populationCoverage = new HashMap <String ,Double>();
				for(String s: subpopulationProportion.keySet()){			
					String s2 = s.substring(1,s.length()-1);
					List<String> mutatedGenes = MyFunc.split(", ", s2);
					if(s2.equals("")){
						continue;
					}
					mutatedGeneCount += mutatedGenes.size()*subpopulationProportion.get(s);
					int driver =0;
					for(String g: mutatedGenes){
						int i = Integer.valueOf(g);
						if(i < S.driverSize){
							driver++;
						}
						if(populationCoverage.containsKey(g)){
							populationCoverage.put(g, populationCoverage.get(g) + subpopulationProportion.get(s));
						}else{
							populationCoverage.put(g, subpopulationProportion.get(s));
						}
					}
					mutatedDriverGeneCount += driver*subpopulationProportion.get(s);
				}
				double founder = 0;
				for(String g:  populationCoverage.keySet()){
					if(populationCoverage.get(g) >= populationCoverageForfounderMutation){
						founder++;
					}
				}
				founderMutationPropotion = founder/mutatedGeneCount;
				
			}

	void getMajorSubpopulationProportion(){
		majorSubpopulationProportion = new LinkedHashMap<String, Double>();
		for(String k: MyFunc.sortKeysByDescendingOrderOfValues(subpopulationProportion)){
			if(subpopulationProportion.get(k) > subpopulationProportionCutoff){
				majorSubpopulationProportion.put(k, subpopulationProportion.get(k));
			}else{
				break;
			}
		}
	}
	
	void getMutationSimilarityMatrix(){
		if(mutationProfileMatrix==null){
			getMutationProfileMatrix();
		}
		mutationSimilalityMatrix = getSimilarityMatrix(mutationProfileMatrix);
	}
	
	void getPopulationEntropy() {
		try {
			if(mutationSimilalityMatrix==null){
				getMutationSimilarityMatrix();
			}
			populationEntropy = getEigenValueEntropy(mutationSimilalityMatrix);
		} catch (NotConvergedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public String toString(){
		String S = "";			
		S += "time\t" + time + "\n";
		S += "populationSize\t" + populationSize + "\n";
		S += "mutatedGeneCount\t" + mutatedGeneCount + "\n";
		S += "mutatedDriverGeneCount\t" + mutatedDriverGeneCount + "\n";
		S += "founderMutationPropotion\t" + founderMutationPropotion + "\n";
		S += "populationGrothRate\t" + populationGrowthRate + "\n";
		S += "populationEntropy\t" + populationEntropy + "\n";
		S += "majorSubpopulationProportion\t" +  majorSubpopulationProportion + "\n";	
		return S;	
	}
	
	public void calculateStatistics(){
		getSubpopulationProportion();
		getMajorSubpopulationProportion();
		getMutatedGeneCount();
		getMutationProfileMatrix();
		getMutationProfileMatrixBySliceSampling();
		getPopulationEntropy();
		getPopulationGrothRate();
		
	}
	
	
	public static void main(String [] args) throws Exception{			
		
	}

}
