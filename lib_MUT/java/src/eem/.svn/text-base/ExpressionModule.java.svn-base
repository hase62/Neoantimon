package eem;

import java.io.*;
import java.sql.*;
import java.util.*;

import processing.core.PVector;

import utility.*;
import utility.HierarchicalClustering.DistFuncType;


public class ExpressionModule implements Serializable{
	private static final long serialVersionUID = -2112591160059022696L;
	private String id;
	private EEM eem;
	private Map <String, Double> activityProfile;
	private List <String> seedGenes;
	private List <String> moduleGenes;
	private ClusteredMyMatWithAnnotation moduleGeneProfiles;
	private ClusteredMyMatWithAnnotation seedGeneProfiles;
	ExpressionModuleSet expressionModuleCluster;	// for ExpressionModuleClustering
	private double minimunPvalue; // -log10 scale
	private double PvalueCorrectedForMultipleTest;  
	private int numberOfMultipleTesting;
	
	private transient MyMat Exp;
	
	public ExpressionModule(ExpressionModule e){
		id = e.id;
		eem = e.eem;
		activityProfile = e.activityProfile;
		seedGeneProfiles = e.seedGeneProfiles;
		moduleGenes = e.moduleGenes;
		seedGenes = e.seedGenes;
		expressionModuleCluster = e.expressionModuleCluster;
		minimunPvalue = e.minimunPvalue;
		PvalueCorrectedForMultipleTest = e.PvalueCorrectedForMultipleTest;
		numberOfMultipleTesting = e.numberOfMultipleTesting;	
	}
	
	private  ExpressionModule(){};
	
	static ExpressionModule getFromResultSet(ResultSet rs, MyMat Exp){
		ExpressionModule em = new ExpressionModule();
		try {
			em.id = rs.getString("id");
			em.PvalueCorrectedForMultipleTest = rs.getDouble("correctedPvalue");
			if(Double.isInfinite(em.PvalueCorrectedForMultipleTest)){
				em.PvalueCorrectedForMultipleTest = Double.MAX_VALUE;
			}
			em.numberOfMultipleTesting = rs.getInt("numberOfMultipleTesting");
			if(rs.getString("eemClass").equals("CoherenceBasedEEM") || rs.getString("eemClass").equals("eem.CoherenceBasedEEM") ){
				em.eem = CoherenceBasedEEM.getFromResultSet(rs);
			}else if (rs.getString("eemClass").equals("BiclusterBasedEEM") || rs.getString("eemClass").equals("eem.BiclusterBasedEEM") ){
				em.eem = BiclusterBasedEEM.getFromResultSet(rs);
			}else{
				throw new UnsupportedOperationException();
			}
			em.minimunPvalue = em.eem.getPvalue();
			if(Double.isInfinite(em.minimunPvalue)){
				em.minimunPvalue = Double.MAX_VALUE;
			}
			em.seedGenes = new ArrayList<String>(em.eem.getSeedGenes());
			em.moduleGenes = new ArrayList<String>(em.eem.getModuleGenes());
			em.setProfiles(Exp);
			return em;	
		} catch (SQLException e){
			e.printStackTrace();
		}
		return null;	
	}
	
	static ExpressionModule getFromString(String s, MyMat Exp){
		ExpressionModule em = new ExpressionModule();
		try{
		List <String> tmp = MyFunc.split("\t", s);
		em.id = tmp.get(0);
		em.PvalueCorrectedForMultipleTest = Double.valueOf(tmp.get(1));
		em.numberOfMultipleTesting = Integer.valueOf(tmp.get(2));
		if(tmp.get(3).startsWith("eem.CoherenceBasedEEM") || tmp.get(3).startsWith("CoherenceBasedEEM")){
			em.eem = CoherenceBasedEEM.getFromString(tmp.get(3));
		}else if (tmp.get(3).startsWith("eem.BiclusterBasedEEM") || tmp.get(3).startsWith("BiclusterBasedEEM")) {
			em.eem = BiclusterBasedEEM.getFromString(tmp.get(3));
		}else{
			throw new UnsupportedOperationException();
		}
		em.minimunPvalue = em.eem.getPvalue();
		em.seedGenes = new ArrayList<String>(em.eem.getSeedGenes());
		em.moduleGenes = new ArrayList<String>(em.eem.getModuleGenes());
		em.setProfiles(Exp);
		return em;
		}catch (Exception e){
			e.printStackTrace();
		}
		return null;
	}	
	
	
	public ExpressionModule(String id, EEM eem, MyMat Exp){
		this.id = id;
		this.eem = eem;
		minimunPvalue = eem.getPvalue();
		PvalueCorrectedForMultipleTest = eem.getPvalue();
		numberOfMultipleTesting = 1;		
		seedGenes = new ArrayList<String>(eem.getSeedGenes());
		moduleGenes = new ArrayList<String>(eem.getModuleGenes());
		setProfiles(Exp);
	}
	
	private void setProfiles(MyMat Exp){
		this.Exp = Exp;
		activityProfile = new HashMap<String, Double>();
		moduleGeneProfiles = new ClusteredMyMatWithAnnotation(Exp.getSubMatrix(moduleGenes, Exp.getColNames()));
		moduleGeneProfiles.normalizeRows();
		seedGeneProfiles = new ClusteredMyMatWithAnnotation(Exp.getSubMatrix(seedGenes, Exp.getColNames()));
		seedGeneProfiles.normalizeRows();
		seedGeneProfiles.setRowDistFunc(HierarchicalClustering.DistFuncType.CORRELATION);
		seedGeneProfiles.setRowClusterDistFunc(HierarchicalClustering.ClusterDistFuncType.COMPLETE);
		seedGeneProfiles.setColDistFunc(HierarchicalClustering.DistFuncType.DISSIMILARITY);
		seedGeneProfiles.setColClusterDistFunc(HierarchicalClustering.ClusterDistFuncType.COMPLETE);
		if(moduleGenes.size() > 0){
		if(moduleGenes.size() == 1){
			for(String s: Exp.getColNames()){
				activityProfile.put(s,Exp.get(moduleGenes.get(0), s));
			}
		}else{
			//if(!isAbsolute()){
				List <Double> tmp =  moduleGeneProfiles.getRowMeans();
				for(int i=0; i < Exp.colSize(); i++){
					activityProfile.put(Exp.getColNames().get(i),tmp.get(i));
				}
			/*}else{
				HierarchicalClustering H = new HierarchicalClustering(new Dist(getModuleGeneProfiles(),'c'));
				H.perform();
				List< List<String>> tmp = H.cutTree(2);
				List <String> subModule1 = tmp.get(0);
				List <String> subModule2 = tmp.get(1);
				List <Double> profile1 = moduleGeneProfiles.getRowMeans(subModule1);
				List <Double> profile2 = moduleGeneProfiles.getRowMeans(subModule2);
				if(subModule1.size() > subModule2.size()){
					for(int i = 0, n = Exp.colSize(); i < n; i++){
						activityProfile.put(Exp.getColNames().get(i), profile1.get(i) -  profile2.get(i));	
					}
				}else{
					for(int i = 0, n = Exp.colSize(); i < n; i++){
						activityProfile.put(Exp.getColNames().get(i), profile2.get(i) -  profile1.get(i));	
					}
				}
			}*/
		}
		}
		expressionModuleCluster = null;
		moduleGeneProfiles.supressColClustering();
		moduleGeneProfiles.reorderCols(MyFunc.sortKeysByDescendingOrderOfValues(activityProfile));
		moduleGeneProfiles.setRowDistFunc(HierarchicalClustering.DistFuncType.DISSIMILARITY);
		moduleGeneProfiles.setRowClusterDistFunc(HierarchicalClustering.ClusterDistFuncType.COMPLETE);
		Map <String, String> tmp = new HashMap<String, String>();
		for(String s: moduleGenes){
			tmp.put(s, "a");
		}
		seedGeneProfiles.setRowAnnotation(new StringMat("module genes",tmp));	
	}
	
	
	public void setId(String s){
		id = s;
	}
	
	public String getId(){
		return id;
	}
	public List <String> getModuleGenes(){
		return moduleGenes;
	}
	public List <String> getSeedGenes(){
		return seedGenes;
	}
	public Double getPvalue(){
		return PvalueCorrectedForMultipleTest;
	}
	public EEM getEEM(){
		return eem;
	}
	//public boolean isAbsolute(){
	//	return eem.isAbsolute();	
	//}
	public Map <String, Double> getActivityProfile(){
		return activityProfile;	
	}
	public ClusteredMyMatWithAnnotation getClusteredModuleGeneProfiles(){
		if(!moduleGeneProfiles.isRowClustered()){
			 performModuleGeneProfileClustering();
		}
		return moduleGeneProfiles;
	}
	public ClusteredMyMatWithAnnotation getClusteredSeedGeneProfiles(){
		if(!seedGeneProfiles.isRowClustered()){
			 performSeedGeneProfileClustering();
		}
		return seedGeneProfiles;
	}
	
	
	public ClusteredMyMatWithAnnotation getModuleGeneProfiles(){
		return moduleGeneProfiles;
	}
	public ClusteredMyMatWithAnnotation getSeedGeneProfiles(){
		return seedGeneProfiles;
	}
	public String toString(){
		return id + "\t" + PvalueCorrectedForMultipleTest + "\t" + numberOfMultipleTesting + "\t" + getEEM();
	}
	
	public void performModuleGeneProfileClustering(){
		moduleGeneProfiles.performClustering();
	}
	public void performSeedGeneProfileClustering(){
		seedGeneProfiles.performClustering();
	}
	
	public void performGeneProfileClustering(){
		performSeedGeneProfileClustering();
		performModuleGeneProfileClustering();
	}
	
	public void setSampleAnnotation(StringMat annotation){
		moduleGeneProfiles.setColAnnotation(annotation);
		seedGeneProfiles.setColAnnotation(annotation);
	}
	public void addSampleAnnotation(StringMat annotation){
		moduleGeneProfiles.addColAnnotation(annotation);
		seedGeneProfiles.addColAnnotation(annotation);
	}

	public void addBiclusterInformation2SampleAnnotation(){
		if(eem instanceof BiclusterBasedEEM){
			List <String> biclusteredSample = ((BiclusterBasedEEM)eem ).getBiclusteredConditions();
			Map <String, String> tmp = new HashMap<String, String>();
			for(String s: biclusteredSample){
				tmp.put(s, "a");
			}
			moduleGeneProfiles.addColAnnotation(new StringMat("biclustered samples", tmp));
			seedGeneProfiles.addColAnnotation(new StringMat("biclustered samples", tmp));	
		}
	}
	
	public void showBiclusteredSeedGeneProfile(){
		/*
		if(eem instanceof BiclusterBasedEEM){
			MyMat Ec = new MyMat(Exp);
			Ec.normalizeRows();
			List <String> biclusteredSample = ((BiclusterBasedEEM)eem ).getBiclusteredConditions();
			Map <String, Double>  biclusteredSampleMean = new HashMap<String, Double>();
			if(biclusteredSample.size() == 1){
				for(String s: seedGenes){
					biclusteredSampleMean.put(s,Ec.get(s, biclusteredSample.get(0)));
				}	
			}else{
				for(String s: seedGenes){
					List <Double> tmp = new ArrayList<Double>();
					for(String t:  biclusteredSample){
					tmp.add(Ec.get(s, t));
					}
					biclusteredSampleMean.put(s,MyFunc.mean(tmp));
				}
			}
			MyMat Eg = new MyMat(Exp);
			Eg.normalizeCols();
			List <String> biclusteredGene = ((BiclusterBasedEEM)eem ).getBiclusteredGenes();
			Map <String, Double>  biclusteredGeneMean = new HashMap<String, Double>();
			if(biclusteredGene.size() == 1){
				for(String s: Exp.getColNames()){
					biclusteredGeneMean.put(s,Eg.get(biclusteredGene.get(0), s));
				}	
			}else{
				for(String s: Exp.getColNames()){
					List <Double> tmp = new ArrayList<Double>();
					for(String t:  biclusteredGene ){
					tmp.add(Eg.get(t, s));
					}
					biclusteredGeneMean.put(s,MyFunc.mean(tmp));
				}
			}	
			
			
			
			seedGeneProfiles = new ClusteredMyMatWithAnnotation(Eg.getSubMatrix(seedGenes, Exp.getColNames()));	
			seedGeneProfiles.supressColClustering();
			seedGeneProfiles.reorderCols(MyFunc.sortKeysByDecendingOrderOfValues(biclusteredGeneMean));	
			seedGeneProfiles.supressRowClustering();
			seedGeneProfiles.reorderRows(MyFunc.sortKeysByDecendingOrderOfValues(biclusteredSampleMean));	
			addBiclusterInformation2SampleAnnotation();
			Map <String, String> tmp = new HashMap<String, String>();
			for(String s: moduleGenes){
				tmp.put(s, "a");
			}
			seedGeneProfiles.setRowAnnotation(new StringMat("module genes",tmp));	
		}
		*/
		addBiclusterInformation2SampleAnnotation();
	}
	
	
	public boolean hasExpressionModuleCluster(){
		return expressionModuleCluster==null?false:true;
	}
	
	public ExpressionModuleSet getExpressionModuleCluster(){
		return expressionModuleCluster;
	}
	
	public static ExpressionModule getBetterExpressionModules(ExpressionModule e1, ExpressionModule e2){
		ExpressionModule e3;
		if(e1.minimunPvalue > e2.minimunPvalue){
			e3 = new ExpressionModule(e1);
			e3.minimunPvalue = e1.minimunPvalue;
		}else{
			e3 = new ExpressionModule(e2);
			e3.minimunPvalue = e2.minimunPvalue;
		}
		e3.numberOfMultipleTesting = e1.numberOfMultipleTesting + e2.numberOfMultipleTesting;
		if(e3.minimunPvalue != Double.MAX_VALUE){
			double p = e3.minimunPvalue;	
			if(p < 5){
				p = Math.pow(10, -p);	
				p = 1- Math.pow(1-p, e3.numberOfMultipleTesting);
				p  = - Math.log10(p);
			}else{
				p = p - Math.log10(e3.numberOfMultipleTesting);
			}
			e3.PvalueCorrectedForMultipleTest = p;
		}
		return e3;
			
	}
	
	
}
