package eem;

import java.io.Serializable;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.*;
import java.util.regex.*;

import eem.CoherenceBasedEEMsearch.DistFunc;
import sun.nio.cs.ext.MacHebrew;
import sun.reflect.Reflection;
import utility.Dist;
import utility.MyException;
import utility.MyFunc;
import utility.MyMat;

public class CoherenceBasedEEM implements EEM, Serializable{
	protected static final long serialVersionUID = 1652144711796321826L;
	protected double absoluteRadius;
	protected double relativeRadius;
	protected List<String> seedGenes;
	protected List<String> moduleGenes;
	
	protected Double Pvalue1;
	protected Double Pvalue2;
	protected Double Pvalue;	
	
	//protected boolean absoluteCorrelation;
	protected int itrForPvalue2Calculation;
	protected int coreGeneSize;
	protected List<Double> center;
	
	protected transient List <Integer> nullDist;
	protected transient List <Integer> densityDist;
	protected transient CoherenceBasedEEMsearch parent;	
	protected transient String centerGene;
	
	
	protected  CoherenceBasedEEM(){};
	
	static CoherenceBasedEEM getFromResultSet(ResultSet rs){
		 CoherenceBasedEEM eem = new CoherenceBasedEEM();	
		 eem.parent = null;
		 eem.center = new ArrayList<Double>();
		 eem.centerGene = null;
		 eem.nullDist = new ArrayList<Integer>();
		 eem.densityDist = new ArrayList<Integer>();
		 try {
			 eem.absoluteRadius = rs.getDouble("absoluteRadius");
			 eem.relativeRadius = rs.getDouble("relativeRadius");	
			 //eem.absoluteCorrelation = rs.getString("absoluteCorrelation").equals("true")?true:false;
			 eem.itrForPvalue2Calculation =  rs.getInt("itrForPvalue2Calculation");
			 eem.seedGenes = MyFunc.split(" ", rs.getString("seedGenes"));
			 eem.moduleGenes = MyFunc.split(" ", rs.getString("moduleGenes"));
			 eem.Pvalue1 = rs.getDouble("Pvalue1");
			 if(rs.wasNull()){ eem.Pvalue1 = null;}
			 else if(eem.Pvalue1.isInfinite()){ eem.Pvalue1 = Double.MAX_VALUE;}
			 eem.Pvalue2 = rs.getDouble("Pvalue2");
			 if(rs.wasNull()){ eem.Pvalue2 = null;}
			 else if(eem.Pvalue2.isInfinite()){ eem.Pvalue2 = Double.MAX_VALUE;}
			 eem.Pvalue = rs.getDouble("Pvalue");
			 if(rs.wasNull()){ eem.Pvalue = null;}
			 else if(eem.Pvalue.isInfinite()){ eem.Pvalue = Double.MAX_VALUE;}
			 eem.coreGeneSize = rs.getInt("coreGeneSize");
			 return eem;
		 } catch (SQLException e) {
			 e.printStackTrace();
		 }
		 return null;
	 }
	
	static CoherenceBasedEEM getFromString(String s){
		CoherenceBasedEEM eem = new CoherenceBasedEEM();
		eem.parent = null;
		eem.center = new ArrayList<Double>();
		eem.centerGene = null;
		eem.nullDist = new ArrayList<Integer>();
		 eem.densityDist = new ArrayList<Integer>();
		try {
		Pattern p;
		Matcher m;
		
		p = Pattern.compile("absoluteRadius=(\\S+)");
		m = p.matcher(s);
		m.find();
		eem.absoluteRadius = Double.valueOf(m.group(1));
		
		p = Pattern.compile("relativeRadius=(\\S+)");
		m = p.matcher(s);
		m.find();
		eem.relativeRadius = Double.valueOf(m.group(1));
		
		//p = Pattern.compile("absoluteCorrelation=(\\S+)");
		//m = p.matcher(s);
		//m.find();
		//eem.absoluteCorrelation = m.group(1).equals("true")?true:false;
		
		p = Pattern.compile("itrForPvalue2Calculation=(\\S+)");
		m = p.matcher(s);
		m.find();
		eem.itrForPvalue2Calculation = Integer.valueOf(m.group(1));
		
		p = Pattern.compile("coreGeneSize=(\\S+)");
		m = p.matcher(s);
		m.find();
		eem.coreGeneSize = Integer.valueOf(m.group(1));
		
		p = Pattern.compile("seedGenes=\\[(.*?)\\]");
		m = p.matcher(s);
		m.find();
		eem.seedGenes = MyFunc.split(", ", m.group(1));
		
		p = Pattern.compile("moduleGenes=\\[(.*?)\\]");
		m = p.matcher(s);
		m.find();
		eem.moduleGenes = MyFunc.split(", ", m.group(1));
		
		p = Pattern.compile("center=\\[(.*?)\\]");
		m = p.matcher(s);
		m.find();
		for(String S:MyFunc.split(", ", m.group(1))){
			eem.center.add(Double.valueOf(S));
		}
		
		
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
	
	CoherenceBasedEEM(CoherenceBasedEEMsearch parent, List <String> gene){
		this.parent = parent;
		absoluteRadius = parent.absoluteRadius;
		relativeRadius = parent.relativeRadius;
		//absoluteCorrelation = parent.absoluteCorrelation;
		itrForPvalue2Calculation =  parent.itrForPvalue2Calculation;
		seedGenes = new ArrayList<String>(gene);
		moduleGenes = new ArrayList<String>();
		center = new ArrayList<Double>();
		centerGene = null;
		Pvalue1 = null;
		Pvalue2 = null;
		Pvalue = null;
		coreGeneSize = parent.coreGeneSize;
		densityDist = new ArrayList<Integer>();
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
		tmp.add("absoluteRadius=" + absoluteRadius);
		tmp.add("relativeRadius=" + relativeRadius);
		//tmp.add("absoluteCorrelation=" + absoluteCorrelation);
		tmp.add("itrForPvalue2Calculation=" + itrForPvalue2Calculation);
		tmp.add("coreGeneSize=" + coreGeneSize);
		tmp.add("seedGenes=" + seedGenes);
		tmp.add("moduleGenes=" + moduleGenes);
		tmp.add("center=" + center);
		tmp.add("Pvalue=" + Pvalue);
		tmp.add("Pvalue1=" + Pvalue1);
		tmp.add("Pvalue2=" + Pvalue2);
		return MyFunc.join(" ", tmp);
	}
	protected  void findCenter(){
		int i,j;	
		moduleGenes = new ArrayList<String>();
		centerGene = null; 
		double distSum = 0;
		for(i =0; i < seedGenes.size(); i++){
			List <String> tmp = new ArrayList<String>();	
			double distSumTmp = 0; 
			for(j =0; j < seedGenes.size(); j++){
				if(i == j){
					tmp.add(seedGenes.get(j));
					continue;
				}
				double tmp2 = parent.distfunc.get(seedGenes.get(i), seedGenes.get(j));
				if(tmp2 < absoluteRadius){
					tmp.add(seedGenes.get(j));
					distSumTmp += tmp2;
				}
			}
			if(moduleGenes.isEmpty() || tmp.size() > moduleGenes.size() || (tmp.size() ==  moduleGenes.size() && distSumTmp < distSum) ){
				moduleGenes  = tmp;
				center = parent.Exp.getRow(seedGenes.get(i));
				centerGene = seedGenes.get(i);
				distSum = distSumTmp;
			}
		}
		Map <String, Double> dist2center  = new HashMap<String, Double>();
		for(i =0; i < moduleGenes.size(); i++){
		dist2center.put(moduleGenes.get(i), parent.distfunc.get(moduleGenes.get(i), centerGene));
		}
		moduleGenes = MyFunc.sortKeysByAscendingOrderOfValues(dist2center);
		
	}
	
	protected void refineCenter() throws MyException{
		int i,j;	
			
		List <String> coreGenes ;
		Map <Double, Integer> centerSum2count = new HashMap<Double, Integer>();
		
		//initialize center and core gene
		if(coreGeneSize < moduleGenes.size()){
			coreGenes = moduleGenes.subList(0, coreGeneSize-1);	
		}else{
			coreGenes = moduleGenes;
		}
		center =  parent.Exp.getRowMeans(coreGenes);
		centerSum2count.put(MyFunc.sum(center), 1);
		
		//update center and core gene;
		List<Double> bestCenter = new ArrayList<Double>();
		List<String> bestModuleGenes = new ArrayList<String>();
		for(int k = 0; k < 10; k++){
			
			Map <String, Double> dist2center  = new HashMap<String, Double>();
			for(i =0; i < seedGenes.size(); i++){
				double tmp = parent.distfunc.get(center, parent.Exp.getRow(seedGenes.get(i)));
				if(tmp < absoluteRadius){
					dist2center.put(seedGenes.get(i), tmp);
				}
			}
			moduleGenes =  MyFunc.sortKeysByAscendingOrderOfValues(dist2center);
			if(moduleGenes .size() > coreGeneSize){
				coreGenes =  moduleGenes.subList(0, coreGeneSize-1);	
			}else{
				coreGenes = moduleGenes;
			}
			center =  parent.Exp.getRowMeans(coreGenes);
			if(centerSum2count.containsKey(MyFunc.sum(center))){
				return;
			}else{
				centerSum2count.put(MyFunc.sum(center), 1);
			}
			if(bestModuleGenes.size() < moduleGenes.size()){
				bestModuleGenes = moduleGenes;
				bestCenter = center;
			}
		}
		moduleGenes = bestModuleGenes;
		center = bestCenter;
	}
		
	
	protected  void refineCenter0() throws MyException{
		int i,j;	
		if(moduleGenes.size() < 3){
			throw  new  MyException("unnable to find an enough number of module genes!");
		}
		List <String> coreGenes ;	
		if(coreGeneSize < moduleGenes.size()){
			Map <String, Double> dist2center  = new HashMap<String, Double>();
			for(i =0; i < moduleGenes.size(); i++){
			dist2center.put(moduleGenes.get(i), parent.distfunc.get(moduleGenes.get(i), centerGene));
			}
		coreGenes = MyFunc.sortKeysByAscendingOrderOfValues(dist2center);
		coreGenes = coreGenes.subList(0, coreGeneSize-1);	
		}else{
			coreGenes = moduleGenes;
		}
		List<List<String>> CoreTriplet = MyFunc.getAllCombination(coreGenes, 3);
		double distSum = 0;
		//if(!absoluteCorrelation){
			for(i=0; i<CoreTriplet.size();i++){
				double distSumTmp = 0;
				List <String> tmp = new ArrayList<String>();
				List<Double> tripletCenter = parent.Exp.getRowMeans(CoreTriplet.get(i));
				for(j =0; j < seedGenes.size(); j++){
					double tmp2 = parent.distfunc.get(tripletCenter, parent.Exp.getRow(seedGenes.get(j)));
					if(tmp2 < absoluteRadius){
						tmp.add(seedGenes.get(j));
						distSumTmp += tmp2;
					}
				}
				if(tmp.size() > moduleGenes.size() || (tmp.size() ==  moduleGenes.size() && distSumTmp < distSum)){
					moduleGenes  = tmp;
					center = tripletCenter;
					distSum = distSumTmp;
				}
			}
		/*}else{
			MyMat E = new MyMat(coreGenes,parent.Exp.getColNames());
			List <Double> centerProfile = parent.Exp.getRow(centerGene);
			for(int k = 0; k < coreGenes.size(); k++){
				List <Double> coreGeneProfile = parent.Exp.getRow(coreGenes.get(k));
				if(parent.Cor.get(centerGene, coreGenes.get(k)) > 0){
					for(int l = 0; l < centerProfile.size(); l++){
						E.set(k, l, coreGeneProfile.get(l));						
					}	
				}else{
					for(int l = 0; l < centerProfile.size(); l++){
						E.set(k, l, -coreGeneProfile.get(l));						
					}	
				}
			}
			for(i=0; i<CoreTriplet.size();i++){
				double distSumTmp = 0;
				List <String> tmp = new ArrayList<String>();
				List<Double> tripletCenter = E.getRowMeans(CoreTriplet.get(i));
				for(j =0; j < seedGenes.size(); j++){
					double tmp2 = parent.distfunc.get(tripletCenter, parent.Exp.getRow(seedGenes.get(j)));
					if(tmp2 < absoluteRadius){
						tmp.add(seedGenes.get(j));
						distSumTmp += tmp2;
					}
				}
				if(tmp.size() > moduleGenes.size() || (tmp.size() ==  moduleGenes.size() && distSumTmp < distSum)){
					moduleGenes  = tmp;
					center = tripletCenter;
					distSum = distSumTmp;
				}
			}
		}*/
	}
	
	
	
	public void findModuleGenes() throws MyException {
		findCenter();
		try {
			refineCenter();
		}
		catch (Exception err) {
		}
	}
	
	protected List <String> getGlobalCoherentGenes(){
		List <String> tmp = new ArrayList<String>();
		for(int j =0; j < parent.allGenes.size(); j++){
			double tmp2 = parent.distfunc.get(center, parent.Exp.getRow(parent.allGenes.get(j)));
			if(tmp2 < absoluteRadius){
				tmp.add(parent.allGenes.get(j));
			}
		}
		return tmp;			
	}
	
	protected void generateNullDist(int n){
		int j= 0;
		while(nullDist.size()<n){
			CoherenceBasedEEM e = new CoherenceBasedEEM(parent, MyFunc.sample(parent.allGenes, seedGenes.size()));
			try{
				e.findModuleGenes();
				//nullDist.add(e.getModuleGenes().size());
				nullDist.add(e.getGlobalCoherentGenes().size());
			}
			catch (MyException err) {
				err.printStackTrace();
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
			double p = MyFunc.calculatePvalueForSetOverlap(parent.allGenes.size(), seedGenes.size(), (int)(parent.allGenes.size()*relativeRadius), moduleGenes.size() );
			Pvalue1 = -Math.log10(p);
		}
		Pvalue = Pvalue1;
	}
	
	
	
	protected Map<Integer, Double> getDensityDist(){
		
		Map<Integer, Double> DistMap = new HashMap<Integer, Double>();
		
		for(Integer I : nullDist){
			
			if(DistMap.containsKey(I)){
				DistMap.put(I, DistMap.get(I) + 1);
			}else{
				DistMap.put(I, 1.0);
			}
		}
		return  DistMap;
		
	}
	
	
	
	
	/*// based on extreme distribution
	public void calculatePvalue2(){
		try{
			generateNullDist(itrForPvalue3Calculation);
			double p = MyFunc.calculatePvalueUsingExtremeDistribusion((double)moduleGenes.size(), nullDist);
			Pvalue2  = (p == 0)?Double.MAX_VALUE: -Math.log10(p);
			Pvalue = Pvalue2;	
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
			Pvalue2 = -Math.log10(p/nullDist.size());
			Pvalue = Pvalue2;
		}	
	}
	*/
	
	
	//based on empirical distribution
	/*public void calculatePvalue2(){
		generateNullDist(itrForPvalue2Calculation);
		double  p = 0;
		for(Integer d: nullDist){
			if(d > moduleGenes.size()){
				p++;
			}
		}
		if (p == 0){
			p = 1;
		}
		Pvalue2 = -Math.log10(p/nullDist.size());
		Pvalue = Pvalue2;
	}*/
   
    
	public void calculatePvalue2(){
		generateNullDist(itrForPvalue2Calculation);
		Map<Integer, Double> Dist =  getDensityDist();	
		double P = 0;
		for(Integer I: Dist.keySet()){
			if(moduleGenes.size() > I){
				continue;				
			}
			double p = MyFunc.calculatePvalueForSetOverlap(parent.allGenes.size(), seedGenes.size(), I, moduleGenes.size() );	
			P += p * Dist.get(I);	
		}
		Pvalue2 = -Math.log10(P/itrForPvalue2Calculation);
		Pvalue = Pvalue2;
	}
	 
	
	public List<String> getModuleGenes() {
		return moduleGenes;
	}
	public void cutParent(){parent = null;}
	public double getAbsoluteRadius(){return absoluteRadius;};
	public double getRelativeRadius(){return relativeRadius;};
	public List<Double>getCenter(){return center;};
	public String getCenterGene(){return centerGene;}; 
	//public boolean isAbsolute(){return  absoluteCorrelation;};
	void setAbsoluteRadius(double r){absoluteRadius = r;}
	public Double getPvalue(){return Pvalue;}
	public Double getPvalue1(){return Pvalue1;}
	public Double getPvalue2(){return Pvalue2;}
	public void setPvalue2(double p){Pvalue2 = p;}
	public void setPvalue1(double p){Pvalue1 = p;}
	public void setPvalue(double p){Pvalue = p;}
	public List<String> getSeedGenes(){return seedGenes;}
	
		
}
