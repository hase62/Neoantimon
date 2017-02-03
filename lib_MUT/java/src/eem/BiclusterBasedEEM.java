package eem;

import java.io.*;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


import eem.BiclusterBasedEEMsearch.BiclusterType;

import utility.MyException;
import utility.MyFunc;
import utility.MyMat;

public class BiclusterBasedEEM implements EEM, Serializable{
	private static final long serialVersionUID = -1995233775200720880L;
	
	private List <String> moduleGenes;
	private List <String> biclusteredGenes;	
	private List <String> biclusteredConditions;
	private List <String> seedGenes;
	private Double Pvalue1;
	private Double Pvalue2;
	private Double Pvalue;	
	private double Tg;
	private double Tc;
	private int maxItrForISA; 
	private BiclusterType biclusterType;
	private int itrForPvalue2Calculation;
	private int itrForConvegence;
	
	private transient BiclusterBasedEEMsearch parent;
	private transient List <Integer> nullDist;
	
	private BiclusterBasedEEM(){};
	
	static BiclusterBasedEEM getFromResultSet(ResultSet rs){
		 BiclusterBasedEEM eem = new BiclusterBasedEEM();	
		 eem.parent = null;
		 eem.nullDist = new ArrayList<Integer>();
		 try {
			 eem.Tg = rs.getDouble("Tg");
			 eem.Tc = rs.getDouble("Tc");	
			
			 String tmp =  rs.getString("biclusterType");
			 if(tmp.equals("Up")){
				 eem.biclusterType = BiclusterType.UP;
			 }else if(tmp.equals("Down")){
				 eem.biclusterType = BiclusterType.DOWN;
			 }else if(tmp.equals("Absolute")){
				 eem.biclusterType = BiclusterType.ABSOLUTE;
			 }else if(tmp.equals("GeneUpButConditionDown")){
				 eem.biclusterType = BiclusterType.gUPcDOWN;
			 }else{
				eem.biclusterType = null;
			 }	 
			 eem.itrForPvalue2Calculation =  rs.getInt("itrForPvalue2Calculation");
			 eem.seedGenes = MyFunc.split(" ", rs.getString("seedGenes"));
			 eem.moduleGenes = MyFunc.split(" ", rs.getString("moduleGenes"));
			 eem.biclusteredGenes = MyFunc.split(" ", rs.getString("biclusteredGenes"));
			 eem.biclusteredConditions = MyFunc.split(" ", rs.getString("biclusteredConditions"));
			 eem.Pvalue1 = rs.getDouble("Pvalue1");
			 if(rs.wasNull()){ eem.Pvalue1 = null;}
			 else if(eem.Pvalue1.isInfinite()){ eem.Pvalue1 = Double.MAX_VALUE;}
			 eem.Pvalue2 = rs.getDouble("Pvalue2");
			 if(rs.wasNull()){ eem.Pvalue2 = null;}
			 else if(eem.Pvalue2.isInfinite()){ eem.Pvalue2 = Double.MAX_VALUE;}
			 eem.Pvalue = rs.getDouble("Pvalue");
			 if(rs.wasNull()){ eem.Pvalue = null;}
			 else if(eem.Pvalue.isInfinite()){ eem.Pvalue = Double.MAX_VALUE;}
			 eem.itrForConvegence = rs.getInt("itrForConvegence");
			 eem.maxItrForISA = rs.getInt("maxItrForISA");
			 return eem;
		 } catch (SQLException e) {
			 e.printStackTrace();
		 }
		 return null;
	 }
	
	static BiclusterBasedEEM getFromString(String s){
		BiclusterBasedEEM eem = new BiclusterBasedEEM();
		eem.parent = null;
		eem.nullDist = new ArrayList<Integer>();
		try {
		Pattern p;
		Matcher m;
		
		p = Pattern.compile("Tc=(\\S+)");
		m = p.matcher(s);
		m.find();
		eem.Tc = Double.valueOf(m.group(1));
		
		p = Pattern.compile("Tg=(\\S+)");
		m = p.matcher(s);
		m.find();
		eem.Tg = Double.valueOf(m.group(1));
		
		p = Pattern.compile("biclusterType=(\\S+)");
		m = p.matcher(s);
		m.find();
		String tmp = m.group(1);
		if(tmp.equals("Up")){
			eem.biclusterType = BiclusterType.UP;
		}else if(tmp.equals("Down")){
			eem.biclusterType = BiclusterType.DOWN;
		}else if(tmp.equals("Absolute")){
			eem.biclusterType = BiclusterType.ABSOLUTE;
		}else if(tmp.equals("GeneUpButConditionDown")){
			eem.biclusterType = BiclusterType.gUPcDOWN;
		}else{
			eem.biclusterType = null;
		}	 
		
		p = Pattern.compile("itrForPvalue2Calculation=(\\S+)");
		m = p.matcher(s);
		m.find();
		eem.itrForPvalue2Calculation = Integer.valueOf(m.group(1));
		
		
		p = Pattern.compile("itrForConvegence=(\\S+)");
		m = p.matcher(s);
		m.find();
		eem.itrForConvegence = Integer.valueOf(m.group(1));
		
		p = Pattern.compile("maxItrForISA=(\\S+)");
		m = p.matcher(s);
		m.find();
		eem.maxItrForISA = Integer.valueOf(m.group(1));
		
		p = Pattern.compile("seedGenes=\\[(.*?)\\]");
		m = p.matcher(s);
		m.find();
		eem.seedGenes = MyFunc.split(", ", m.group(1));
		
		p = Pattern.compile("moduleGenes=\\[(.*?)\\]");
		m = p.matcher(s);
		m.find();
		eem.moduleGenes = MyFunc.split(", ", m.group(1));
		
		p = Pattern.compile("biclusteredGenes=\\[(.*?)\\]");
		m = p.matcher(s);
		m.find();
		eem.biclusteredGenes = MyFunc.split(", ", m.group(1));
		
		p = Pattern.compile("biclusteredConditions=\\[(.*?)\\]");
		m = p.matcher(s);
		m.find();
		eem.biclusteredConditions = MyFunc.split(", ", m.group(1));
		
		p = Pattern.compile(".Pvalue=(\\S+)");
		m = p.matcher(s);
		m.find();
		eem.Pvalue  =m.group(1).equals("null")?null:Double.valueOf(m.group(1));
		
		p = Pattern.compile(".Pvalue1=(\\S+)");
		m = p.matcher(s);
		m.find();
		eem.Pvalue1  =m.group(1).equals("null")?null:Double.valueOf(m.group(1));
		
		p = Pattern.compile(".Pvalue2=(\\S+)");
		m = p.matcher(s);
		m.find();
		eem.Pvalue2  =m.group(1).equals("null")?null:Double.valueOf(m.group(1));
		
		
		return eem;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}	
	
	
	
	BiclusterBasedEEM(BiclusterBasedEEMsearch parent, List <String> seedGenes){
		this.parent = parent;
		Tg = parent.Tg;
		Tc = parent.Tc;
		biclusterType = new BiclusterType(parent.biclusterType);
		itrForPvalue2Calculation =  parent.itrForPvalue2Calculation;
		maxItrForISA = parent.maxItrForISA;
		this.seedGenes =  new ArrayList<String>(seedGenes);
		//biclusteredGenes = new ArrayList<String>(seedGenes);
		biclusteredConditions = new ArrayList<String>();
		moduleGenes  = new ArrayList<String>();
		Pvalue1 = null;
		Pvalue2 = null;
		Pvalue = null;
		if(parent.nullDistrbutionForRecycle !=null){
			if(parent.nullDistrbutionForRecycle.containsKey(seedGenes.size())){
				nullDist = parent.nullDistrbutionForRecycle.get(seedGenes.size());
			}else{
				nullDist = new ArrayList<Integer>();
				parent.nullDistrbutionForRecycle.put(seedGenes.size(), nullDist);
			}
		}else{	
			nullDist = new ArrayList<Integer>();
		}
		
	}
	public String toString(){
		List <String> tmp = new ArrayList<String>();
		tmp.add(getClass().getName());
		tmp.add("Tg=" + Tg);
		tmp.add("Tc=" + Tc);
		tmp.add("biclusterType=" + biclusterType);
		tmp.add("itrForPvalue2Calculation=" + itrForPvalue2Calculation);
		tmp.add("itrForConvegence=" + itrForConvegence);
		tmp.add("maxItrForISA=" + maxItrForISA);
		tmp.add("seedGenes=" + seedGenes);
		tmp.add("moduleGenes=" + moduleGenes);
		tmp.add("biclusteredGenes="  + biclusteredGenes);
		tmp.add("biclusteredConditions=" + biclusteredConditions);
		tmp.add("Pvalue=" + Pvalue);
		tmp.add("Pvalue1=" + Pvalue1);
		tmp.add("Pvalue2=" + Pvalue2);
		return MyFunc.join(" ", tmp);
	}
	
	/*
	private void getNextGenes(){
		List <Double> v = new ArrayList<Double>();
		int i,j;
		for(i=0;i<Ec.rowSize();i++){
			double tmp = 0;
			for(j=0; j < biclusteredConditions.size(); j++){
				tmp += Ec.get(i, Ec.colIndexOf(biclusteredConditions.get(j)));
			}
			v.add(tmp/biclusteredConditions.size());
		}
		List <String> nextGenes = new ArrayList<String>();
		
		for(i=0;i<v.size();i++){
			if(v.get(i) >= Tg2){
				nextGenes.add(Ec.getRowNames().get(i));
			}
		}
		if(nextGenes.isEmpty()){
			throw new MyException("nextGenes is empty!");	
		}	
		biclusteredGenes = nextGenes;
	}
	private void getNextConditions(){
		List <Double> v = new ArrayList<Double>();
		int i,j;
		for(j=0;j<Eg.colSize();j++){
			double tmp = 0;
			for(i=0; i < biclusteredGenes.size(); i++){
				tmp += Eg.get(Eg.rowIndexOf(biclusteredGenes.get(i)),j);
			}
			v.add(tmp/biclusteredGenes.size());
		}
		List <String> nextConditions = new ArrayList<String>();
		Tc2 = MyFunc.percentile(v, 1-Tc);
		for(j=0;j<v.size();j++){
			if(v.get(j) >= Tc2){
				nextConditions.add(Eg.getColNames().get(j));
			}
		}
		if(nextConditions.isEmpty()){
			throw new MyException("nextConditions is empty!");	
		}	
		biclusteredConditions = nextConditions;
	}
	*/
	
	
	private void getNextGenes(){
		List <Double> v = new ArrayList<Double>();
		int i,j;
		for(i=0;i<parent.Ec.rowSize();i++){
			double tmp = 0;
			for(j=0; j < biclusteredConditions.size(); j++){
				tmp += parent.Ec.get(i, parent.Ec.colIndexOf(biclusteredConditions.get(j)));
			}
			v.add(tmp/biclusteredConditions.size());
		}
		List <String> nextGenes = new ArrayList<String>();
		
		double threthold = MyFunc.percentile(v, 1-Tg);
		
		
		for(i=0;i<v.size();i++){
			if(v.get(i) >= threthold){
				nextGenes.add(parent.Ec.getRowNames().get(i));
			}
		}
		if(nextGenes.isEmpty()){
			throw new MyException("nextGenes is empty!");	
		}	
		biclusteredGenes = nextGenes;
	}
	private void getNextConditions(){
		List <Double> v = new ArrayList<Double>();
		int i,j;
		for(j=0;j<parent.Eg.colSize();j++){
			double tmp = 0;
			for(i=0; i < biclusteredGenes.size(); i++){
				tmp += parent.Eg.get(parent.Eg.rowIndexOf(biclusteredGenes.get(i)),j);
			}
			v.add(tmp/biclusteredGenes.size());
		}
		List <String> nextConditions = new ArrayList<String>();
		double threthold = MyFunc.percentile(v, 1-Tc);
		for(j=0;j<v.size();j++){
			if(v.get(j) >= threthold){
				nextConditions.add(parent.Eg.getColNames().get(j));
			}
		}
		if(nextConditions.isEmpty()){
			throw new MyException("nextConditions is empty!");	
		}	
		biclusteredConditions = nextConditions;
	}
	
	
	
	public void findModuleGenes(){
		List <Set<String>> beforeGeneSets = new ArrayList<Set<String>>();
		beforeGeneSets.add(new HashSet<String>());
		beforeGeneSets.add(new HashSet<String>());
		beforeGeneSets.add(new HashSet<String>());
		beforeGeneSets.add(new HashSet<String>());
		beforeGeneSets.add(new HashSet<String>());
		beforeGeneSets.add(new HashSet<String>());
		Set <String> afterGenes = new HashSet<String>();
		
		
		if(parent.coEEM != null){
		//Limit seed Genes to a coherent subset
		CoherenceBasedEEM coEEM = new CoherenceBasedEEM(parent.coEEM, seedGenes);
		coEEM.findCenter();
		biclusteredGenes = coEEM.getModuleGenes();
		}else{
			biclusteredGenes = new ArrayList<String>(seedGenes);
		}
		for(int i=0;i<maxItrForISA;i++){
			getNextConditions();
			getNextGenes();
			afterGenes = new HashSet<String>(biclusteredGenes);
			if(beforeGeneSets.get(0).equals(afterGenes)){
				itrForConvegence = i+1;
				moduleGenes = MyFunc.isect(seedGenes, biclusteredGenes);
				
				return;
			}
			for(int j = beforeGeneSets.size() -1 ; j > 0; j--){
				beforeGeneSets.set(j, beforeGeneSets.get(j-1));
			}
			beforeGeneSets.set(0,afterGenes);	
		}
		itrForConvegence = maxItrForISA;
		for(int j = 1, n = beforeGeneSets.size(); j < n; j++){
			if(beforeGeneSets.get(j).equals(afterGenes)){
				moduleGenes = MyFunc.isect(seedGenes, biclusteredGenes);
				return;
			}
		}		
		
		throw new MyException("not converged!");	
	}
	
	private void generateNullDist(int n){
		int j= 0;
		while(nullDist.size()<n){
			BiclusterBasedEEM e = new BiclusterBasedEEM(parent, MyFunc.sample(parent.allGenes, seedGenes.size()));
			try{
				e.findModuleGenes();
				nullDist.add(((Integer)e.getModuleGenes().size()));	
			}
			catch (MyException err) {
				j++;
				if(j >= 10){
					throw new MyException("Cannot generate null distribution!");			
				}
			}		
		}
	}
	public void calculatePvalue1(){
		if(moduleGenes.isEmpty()){
			Pvalue1 = 0.0;
		}else{
			double p = MyFunc.calculatePvalueForSetOverlap(parent.allGenes.size(), seedGenes.size(),biclusteredGenes.size(), moduleGenes.size() );
			Pvalue1 = -Math.log10(p);
		}
		Pvalue = Pvalue1;
	}
	
	public void calculatePvalue2(){
		generateNullDist(itrForPvalue2Calculation);
		double  p = 0;
		for(Integer d: nullDist){
			if(d >= moduleGenes.size()){
				p++;
			}
		}
		if (p == 0){
			p = Math.log10(itrForPvalue2Calculation);
		}else{
			p = -Math.log10(p/nullDist.size());
		}
		Pvalue2 = p;
		Pvalue = Pvalue2;
	}
	
	
	
	
	
	/*public void calculatePvalue3(){
		try{
			generateNullDist(itrForPvalue3Calculation);
			double p = MyFunc.calculatePvalueUsingExtremeDistribusion((double)moduleGenes.size(), nullDist);
			Pvalue3  = (p == 0)?Double.MAX_VALUE: -Math.log10(p);
			Pvalue = Pvalue3;	
		}catch (Exception e) {
			System.err.println("Cannot calculate a P value! based on a extreme value distribution!");
			double  p = 0;
			for(Double d: nullDist){
				if(d > moduleGenes.size()){
					p++;
				}
			}
			if (p == 0){
				p = 1;
			}
			Pvalue3 = -Math.log10(p/nullDist.size());
			Pvalue = Pvalue3;
		}	
	}*/
	
	
	public double getTg(){ return Tg;}
	public double getTc(){ return Tc;}
	public BiclusterType getBiclusterType(){return biclusterType;}
	public boolean isAbsolute(){
		if(biclusterType.equals(BiclusterType.ABSOLUTE)){
			return true;
		}else{
			return false;
		}
	}
	public void cutParent(){parent = null;}
	public List<String> getModuleGenes(){return moduleGenes;}
	public List<String> getSeedGenes(){return seedGenes;}
	public Double getPvalue(){return Pvalue;}
	public Double getPvalue1(){return Pvalue1;}
	public Double getPvalue2(){return Pvalue2;}
	public void setPvalue2(double p){Pvalue2 = p;}
	public void setPvalue1(double p){Pvalue1 = p;}
	public void setPvalue(double p){Pvalue = p;}
	public List<String> getBiclusteredGenes(){return biclusteredGenes;}
	public List<String> getBiclusteredConditions(){return biclusteredConditions;}
	
	
}
