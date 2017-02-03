package eem;

import java.util.*;

import utility.*;

public class generateEEMnulldist {

	
	public static void main(String[] args) throws Exception {
		MyMat E = new MyMat(args[0]);
		CoherenceBasedEEMsearch eem = new CoherenceBasedEEMsearch(E);
		//BiclusterBasedEEMsearch eem = new BiclusterBasedEEMsearch(E);
		//GeneSetExpressionCoherenceAnalysis eem = new GeneSetExpressionCoherenceAnalysis(E);
		
		
		int n = Integer.valueOf(args[1]);
		Map <String, List<String>> geneset = new HashMap<String,List<String>>();
		for(int i = 0; i < 10000; i++){
			List<String> tmp = MyFunc.sample(E.getRowNames(),n);
			geneset.put("g" + i, tmp);
		}
		eem.setGeneSets(geneset);
		eem.calculateCor();
		eem.setRelativeRadius(0.05);
		eem.setEEM();
		eem.findModuleGenes();
		
		for(String s: geneset.keySet()){
			System.out.println(eem.eems.get(s).getModuleGenes().size());
			//System.out.println(eem.getCorMean(geneset.get(s)));
		}
		
		MatrixInformationEnrichmentAnalysis  S = new MatrixInformationEnrichmentAnalysis (new MyMat("/home/niiyan/SVD/breastMiller.tab"));
		
		

		
		
		for(int i = 0; i < 10000; i++){
			List <String> tmp = MyFunc.sample(S.allGenes, 100);
			System.out.println(S.getInformationContent(tmp));
		
		}

	}
	
	
}
