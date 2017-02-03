package eem;

import java.io.*;
import java.util.*;
import java.util.zip.DataFormatException;
import java.sql.*;

import eem.BiclusterBasedEEMsearch.BiclusterType;

import utility.*;

public class ExpressionModuleSet implements Serializable{
	private static  final long serialVersionUID = -6969609459514492054L;
	private Map <String, ExpressionModule> expressionModules;
	private List <String> ids;
	private StringMat sampleAnnotation = null;
	public ExpressionModuleSet(){
		expressionModules = new HashMap<String, ExpressionModule>();
		ids = new ArrayList<String>();
	}
	public ExpressionModuleSet(ExpressionModuleSet expModSet){
		expressionModules = new HashMap<String, ExpressionModule>(expModSet.expressionModules);
		ids = new ArrayList<String>(expModSet.ids);
		sampleAnnotation = expModSet.sampleAnnotation;
	}
	public static ExpressionModuleSet getFromMySQL(String tableName){
		try {
			ExpressionModuleSet ems = new ExpressionModuleSet();	
			List <String> tmp  = MyFunc.split("__", tableName);
			MyMat Exp = Expression.getFromMySQL(tmp.get(0));
			Connection con = MyFunc.getMySQLConnection("expression_module");
			Statement stmt = con.createStatement();
			ResultSet rs = stmt.executeQuery("select * from " + tableName );
			while(rs.next()){
				ems.add(ExpressionModule.getFromResultSet(rs, Exp));
			}
			return ems;
		} catch (SQLException e) {
			e.printStackTrace();
		}
		return null;
	}
	
	public static ExpressionModuleSet getFromTextFile(String eemFile, String expFile){
		try{
			MyMat Exp = new MyMat(expFile);
			return getFromTextFile(eemFile, Exp);
		} catch (Exception e) {
			e.printStackTrace();
		}	
		return null;		
	}
	public static ExpressionModuleSet getFromTextFile(String eemFile, MyMat Exp){
		try{
			ExpressionModuleSet ems = new ExpressionModuleSet();
			List <String> lines = MyFunc.readStringList(eemFile);
			for(String s: lines){
				ems.add(ExpressionModule.getFromString(s, Exp));
			}
			return ems;
		} catch (Exception e) {
			e.printStackTrace();
		}	
		return null;	
	}
	
	public void clear(){
		expressionModules = new HashMap<String, ExpressionModule>();
		ids = new ArrayList<String>();
		sampleAnnotation = null;
	}
	public ExpressionModuleSet getSubSet(List <String> ids){
		ExpressionModuleSet e = new ExpressionModuleSet();
		for(String s: ids){
			if(contains(s)){
				e.add(get(s));
			}
		}
		return e;
	}
	public boolean contains(String id){
		return ids.contains(id);
	}
	
	public boolean isEmpty(){
		return ids.isEmpty();
	}
	public int size(){
		return ids.size();
	}	
	public void add(ExpressionModule expMod){
		String id = expMod.getId();
		if(!expressionModules.containsKey(id)){
			expressionModules.put(id, expMod);
			if(ids.isEmpty()){
				ids.add(id);
			}else{
				for(int i = 0, n = ids.size(); i < n;i++){
					if(expMod.getPvalue() > get(i).getPvalue()){
						ids.add(i, id);
						return;
					}	
				}				
				ids.add(id);
				return;
			}
		}else{
			ExpressionModule e = ExpressionModule.getBetterExpressionModules(expressionModules.get(id),expMod);
			expressionModules.put(id, e);
			ids.remove(id);
			for(int i = 0, n = ids.size(); i < n;i++){
				if(e.getPvalue() > get(i).getPvalue()){
					ids.add(i, id);
					return;
				}	
			}				
			ids.add(id);
			return;
		}
	}
	
	public Map <String, Double> getPvalues(){
		Map <String, Double> P = new HashMap<String, Double>();
		for(String s: ids){
			P.put(s,get(s).getPvalue());
		}
		return P;
	}
	
	public List <Double> getPvalueList(){
		List <Double> P = new ArrayList<Double>();
		for(String s: ids){
			P.add(get(s).getPvalue());
		}
		return P;
	}
	
	public ExpressionModuleSet getSignificantSubset(double PvalueCutoff){
		ExpressionModuleSet ems = new ExpressionModuleSet();
		for(String s: ids){
			if(expressionModules.get(s).getPvalue() >= PvalueCutoff){
				ems.add(get(s));
			}
		}
		return ems;
	}
	
	
	public void addAll(ExpressionModuleSet expModSet){
		for(ExpressionModule e: expModSet.asList()){
			add(e);	
		}
	}
	public void remove(String id){
			expressionModules.remove(id);
			ids.remove(id);
	}
	public void remove(int i){
		remove(ids.get(i));
	}
	public ExpressionModule get(String id){
		return expressionModules.get(id);
	}
	public ExpressionModule get(int i){
		return expressionModules.get(ids.get(i));
	}
	
	public List <String> getIds(){
		return ids;
	}
	public String toString(){
		List <String> tmp  = new ArrayList<String>();
		for(String s: ids){
			ExpressionModule e = expressionModules.get(s);
			tmp.add(e.toString());
		}
		return MyFunc.join("\n", tmp);
	}
	public String toStringWithModuleGenes(){
		List <String> tmp  = new ArrayList<String>();
		for(String s: ids){
			ExpressionModule e = expressionModules.get(s);
			tmp.add(e.toString() + "\t" + MyFunc.join("\t", e.getModuleGenes()));
		}
		return MyFunc.join("\t", tmp);
	}
	
	public ClusteredMyMatWithAnnotation getActivityProfiles(){
		Set <String>  samples= new HashSet<String>();
		for(String s: ids){
			samples.addAll(get(s).getActivityProfile().keySet());
		}
		List <String> samples2 = new ArrayList<String>(samples);
		Collections.sort(samples2);
		ClusteredMyMatWithAnnotation profiles  = new ClusteredMyMatWithAnnotation(ids,  samples2);
		for(String s: ids){
			Map <String, Double> tmp = get(s).getActivityProfile();
			for(String t: profiles.getColNames()){
				profiles.set(s, t, tmp.get(t));
			}
		}
		profiles.setRowDistFunc(HierarchicalClustering.DistFuncType.CORRELATION);
		profiles.setRowClusterDistFunc(HierarchicalClustering.ClusterDistFuncType.COMPLETE);
		profiles.setColDistFunc(HierarchicalClustering.DistFuncType.DISSIMILARITY);
		profiles.setColClusterDistFunc(HierarchicalClustering.ClusterDistFuncType.COMPLETE);
		if(sampleAnnotation != null){
			profiles.setColAnnotation(sampleAnnotation);
		}
		return profiles;
	}
	
	public ClusteredMyMatWithAnnotation getBiclusterProfiles(){
		Set <String>  samples= new HashSet<String>();
		for(String s: ids){
			samples.addAll(get(s).getActivityProfile().keySet());
		}
		List <String> samples2 = new ArrayList<String>(samples);
		Collections.sort(samples2);
		ClusteredMyMatWithAnnotation profiles  = new ClusteredMyMatWithAnnotation(ids,  samples2);
		for(String s: ids){
			if(get(s).getEEM()instanceof BiclusterBasedEEM){
				if(((BiclusterBasedEEM) (get(s).getEEM())).getBiclusterType() == BiclusterType.UP){
					List <String> tmp = ((BiclusterBasedEEM) (get(s).getEEM())).getBiclusteredConditions();					
					for(String t: tmp){
						profiles.set(s, t, 1);
					}
				}
				if(((BiclusterBasedEEM) (get(s).getEEM())).getBiclusterType() == BiclusterType.DOWN){
					List <String> tmp = ((BiclusterBasedEEM) (get(s).getEEM())).getBiclusteredConditions();					
					for(String t: tmp){
						profiles.set(s, t, -1);
					}
				}
			}
		}
		profiles.setRowDistFunc(HierarchicalClustering.DistFuncType.DISSIMILARITY);
		profiles.setRowClusterDistFunc(HierarchicalClustering.ClusterDistFuncType.COMPLETE);
		profiles.setColDistFunc(HierarchicalClustering.DistFuncType.DISSIMILARITY);
		profiles.setColClusterDistFunc(HierarchicalClustering.ClusterDistFuncType.COMPLETE);
		if(sampleAnnotation != null){
			profiles.setColAnnotation(sampleAnnotation);
		}
		return profiles;
	}
	
	
	public ClusteredMyMatWithAnnotation getClusteredActivityProfiles(){
		ClusteredMyMatWithAnnotation tmp = getActivityProfiles();
		tmp.performClustering();
		return tmp;
	}
	
	public ClusteredMyMatWithAnnotation getClusteredBiclusterProfiles(){
		ClusteredMyMatWithAnnotation tmp = getBiclusterProfiles();
		tmp.performClustering();
		return tmp;
	}
	
	
	
	public static ExpressionModuleSet ReadFromFile(String file) throws Exception{
		ObjectInputStream in = new ObjectInputStream(new BufferedInputStream(new FileInputStream(file)));
		return (ExpressionModuleSet) in.readObject();
	}
	
	public void writeToFile(String file) throws Exception{
		ObjectOutputStream out = new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream(file)));
		out.writeObject(this);
		out.close();	
	}
	
	public void performModuleGeneProfileClustering(){
		for(ExpressionModule e: expressionModules.values()){
			e.performModuleGeneProfileClustering();	
		}
	}
	public  void performSeedGeneProfileClustering(){
		for(ExpressionModule e: expressionModules.values()){
			e.performSeedGeneProfileClustering();
		}
	}
	public void performGeneProfileClustering(){
		performSeedGeneProfileClustering();
		performModuleGeneProfileClustering();
	}
	public void performExpressionModuleClustering(double PvalueCutoff, int bgGeneNumber){
		ExpressionModuleClustering C = new ExpressionModuleClustering(this);
		C.setBgGeneNumber(bgGeneNumber);
		C.setPvalueCutoff(PvalueCutoff);
		C.perform();
	}
	public void setSampleAnnotation(StringMat annotation){
		for(ExpressionModule e: expressionModules.values()){
			e.setSampleAnnotation(annotation);
		}
		sampleAnnotation = annotation;
	}
	public void addSampleAnnotation(StringMat annotation){
		for(ExpressionModule e: expressionModules.values()){
			e.addSampleAnnotation(annotation);
		}
		if(sampleAnnotation != null){
			sampleAnnotation = sampleAnnotation.bind(annotation);
		}else{
			sampleAnnotation = annotation;
		}
	}
	
	public List <ExpressionModule> asList(){
		List <ExpressionModule> tmp = new ArrayList<ExpressionModule>();
		for(String s: ids){
			tmp.add(expressionModules.get(s));
		}
		return tmp;
	}
	public List <String> getModuleGeneUnion(){
		Set <String> tmp = new HashSet<String>();
		for(ExpressionModule e: asList()){
			tmp.addAll(e.getModuleGenes());	
		}
		return new ArrayList<String>(tmp);	
	}
	public List <String> getModuleGeneIntersection(){
		Set <String> tmp = new HashSet<String>(asList().get(0).getModuleGenes());
		for(int i = 1, n = size(); i < n; i++){
			tmp.retainAll(asList().get(0).getModuleGenes());	
		}
		return new ArrayList<String>(tmp);	
	}
	
	public ClusteredMyMatWithAnnotation getModuleGeneUnionProfiles(){
		Set <String>  samples= new HashSet<String>();
		for(String s: ids){
			samples.addAll(get(s).getActivityProfile().keySet());
		}
		List <String> samples2 = new ArrayList<String>(samples);
		Collections.sort(samples2);
		List <String> moduleGeneUnion  = getModuleGeneUnion();
		ClusteredMyMatWithAnnotation profiles  = new ClusteredMyMatWithAnnotation(moduleGeneUnion,  samples2);
		for(String s: moduleGeneUnion){
			for(String t: getIds()){
				if(get(t).getModuleGenes().contains(s)){
					Map <String, Double> tmp = get(t).getModuleGeneProfiles().getRowMap(s);
					for(String r: profiles.getColNames()){
						profiles.set(s, r, tmp.get(r));
					}
					break;
				}
			}
		}
		for(int i = size()-1; i >=0; i--){
			String t = getIds().get(i);
			Map <String, String> tmp = new HashMap<String, String>();
			for(String s: moduleGeneUnion){
				if(get(t).getModuleGenes().contains(s)){
					tmp.put(s, "a");
				}
			}
			profiles.addRowAnnotation(new StringMat(t ,tmp));
		}
		profiles.setRowDistFunc(HierarchicalClustering.DistFuncType.CORRELATION);
		profiles.setRowClusterDistFunc(HierarchicalClustering.ClusterDistFuncType.COMPLETE);
		profiles.setColDistFunc(HierarchicalClustering.DistFuncType.DISSIMILARITY);
		profiles.setColClusterDistFunc(HierarchicalClustering.ClusterDistFuncType.COMPLETE);
		profiles.performClustering();
		if(sampleAnnotation != null){
			profiles.setColAnnotation(sampleAnnotation);
		}
		return profiles;
	}
	
	public ClusteredMyMatWithAnnotation getClusteredModuleGeneUnionProfiles(){
		ClusteredMyMatWithAnnotation tmp = getModuleGeneUnionProfiles();
		tmp.performClustering();
		return tmp;
	}
	
	public Map <String, List<String>> getModuleGenes(){
		 Map <String, List<String>> tmp = new HashMap<String, List<String>>();
		 for(String s: getIds()){
			tmp.put(s,get(s).getModuleGenes());			 
		 }
		 return tmp;	
	}
}
