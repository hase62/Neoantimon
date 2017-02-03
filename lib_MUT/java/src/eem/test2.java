package eem;

import java.applet.Applet;
import java.awt.BorderLayout;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.JFrame;
import utility.*;

public class test2   extends JFrame  {
	public test2(Applet C) throws Exception{
			setLayout(new BorderLayout());
			add(C,BorderLayout.CENTER);
			C.init();
			pack();
			setLocation(100, 100);
			setVisible(true);
	}
	public static void main(String[] args) throws Exception {
		MyMat M = new MyMat(args[0]);
		
		//M.reorderRows(MyFunc.sample(M.getRowNames(), 3000));
		
		
		Map <String,List<String>> geneset = MyFunc.readGeneSetFromGmtFile(args[1]);
		
		
		
		BiclusterBasedEEMsearch B = new BiclusterBasedEEMsearch(M);
		B.setMinGeneSetSize(10);
		B.setGeneSets(geneset);
		B.setItrForPvalue2Calculation(10);
		//B.setBiclusterType2Down();
		//B.setBiclusterType2Absolute();
		B.setTg(0.1);
		B.setTc(0.1);
	
		B.perform();
		
		
		ExpressionModuleSet ems = B.getExpressionModuleSet();
		
		
		ExpressionModule em = ems.get(0);
		
		
		/*
		em.addBiclusterInformation2SampleAnnotation();
		
		em.performGeneProfileClustering();
		
		ClusteredMyMatViewer V = em.getSeedGeneProfiles().getViewer();
		
		V.scaleColorByRow();
		
		new test2(V);
		*/
		
		
		List <String> bgene = ((BiclusterBasedEEM)em.getEEM()).getBiclusteredGenes();
		List <String> bcon = ((BiclusterBasedEEM)em.getEEM()).getBiclusteredConditions();
	 
		//M.reorderRows(MyFunc.sample(M.getRowNames(), 1000));
		M.reorderRows(em.getSeedGenes());
		ClusteredMyMatWithAnnotation C = new ClusteredMyMatWithAnnotation(M);
		
		Map <String, String> tmp = new HashMap<String, String>();
		for(String s: bgene){
			tmp.put(s, "a");
		}
		Map <String, String> tmp2 = new HashMap<String, String>();
		for(String s: bcon){
			tmp2.put(s, "a");
		}
		
		C.addRowAnnotation(new StringMat("a", tmp));
		C.addColAnnotation(new StringMat("b", tmp2));
		C.performClustering();
		
		ClusteredMyMatViewer V = new ClusteredMyMatViewer(C);
		V.scaleColorByRow();
		new test2(V);
		
		
		//List<String> tmp = MyFunc.sample(M.getRowNames(),100);
		//M = M.getSubMatByRow(tmp);
		//M.filterRowByVariance(500);
		
		
		
		
		//StringMat S = new StringMat("test_annot.tab");
		
		/*
		SampleLableCorrelation SLC = new SampleLableCorrelation(M);
		SLC.setSampleLabel(S.getColMap("ER"));
		SLC.calculateCorScore();
		SLC.calculatePvalue();
		SLC.calculateQvalue();
		
		Map <String, Double> C = SLC.getCorScore();
		M = M.getSubMatByRow(MyFunc.sortKeysByDecendingOrderOfValues(C).subList(0, 100));
		*/
		
		/*
		CoxRegression Cox = new CoxRegression(M);
		Cox.setTime(S.getColMap("time"));
		Cox.setStatus(S.getColMap("status"));
		Cox.perform();
		Map <String, Double> Z = Cox.getZscore();
		M = M.getSubMatByRow(MyFunc.sortKeysByDecendingOrderOfValues(Z).subList(0, 50));
			
		List <String> l = new ArrayList<String>();
		l.add("p53");
		l.add("ER");
		l.add("grade");
		S = S.getSubMatByCol(l);
		
		
		ClusteredMyMatWithAnnotation M2 = new ClusteredMyMatWithAnnotation(M);
		M2.normalizeRows();
		M2.addColAnnotation(S);
		M2.setRowDistFunc(HierarchicalClustering.DistFuncType.DISSIMILARITY);
		M2.setColDistFunc(HierarchicalClustering.DistFuncType.DISSIMILARITY);
		M2.setRowClusterDistFunc(HierarchicalClustering.ClusterDistFuncType.AVERAGE);
		M2.setColClusterDistFunc(HierarchicalClustering.ClusterDistFuncType.AVERAGE);
		M2.performClustering();
		ClusteredMyMatViewer V = M2.getViewer();
		V.scaleColorByRow();
		
		new test2(V);
		
		*/
		
		
		/*
		CoherenceBasedEEMsearch C = new CoherenceBasedEEMsearch(M);
		C.setSeedGeneSets("E2F_ER_17q.gmt");
		C.calculateCor();
		C.setRelativeRadius(0.03);
		C.performEEMsearch();
		ExpressionModuleSet ems = C.getExpressionModuleSet();
		ems.writeToFile("test.ob");
		ExpressionModuleSet ems2 = ExpressionModuleSet.ReadFromFile("test.ob");
		System.err.println(ems2 + "\n");
		ExpressionModule e = ems2.get("E2F");
		ClusteredMyMatWithAnnotation M2 = e.getModuleGeneProfiles();
	    M2.setRowClusterDistFunc(HierarchicalClustering.ClusterDistFuncType.COMPLETE);
		//M.setColClusterDistFunc(HierarchicalClustering.ClusterDistFuncType.AVERAGE);
		M2.performClustering();
		StringMat S = new StringMat("test_annot.tab");
		List <String> l = new ArrayList<String>();
		l.add("p53");
		l.add("ER");
		l.add("grade");
		S = S.getSubMatByCol(l);
		M2.addColAnnotation(S);
		ClusteredMyMatViewer V = M2.getViewer();
		V.scaleColorByRow();
		
		new test2(V);
		*/
		
	}

}