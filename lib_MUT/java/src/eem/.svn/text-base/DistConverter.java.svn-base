package eem;

import java.util.*;
import utility.*;

public class DistConverter{
	private int maxItrForBisection ;
	private int maxItrForBisection2 ;
	private double  deltaForBisection;
	private Double upperAbsoluteDist;
	private Double lowerAbsoluteDist;
	private CoherenceBasedEEMsearch parent;
	private CoherenceBasedEEM e;
	private Map <String, Double> absolute2relative;
	public DistConverter(CoherenceBasedEEMsearch parent){
		this.parent = parent;
		maxItrForBisection = 30;
		maxItrForBisection2 = 5;
		deltaForBisection = 0.05;
		upperAbsoluteDist  = null;
		lowerAbsoluteDist  = null;
		absolute2relative = new HashMap<String, Double>();
	}
	public void setUpperAbsoluteDist(double d){
		upperAbsoluteDist = d;
	}
	public void setLowerAbsoluteDist(double d){
		lowerAbsoluteDist = d;
	}
	
	
	public double convertAbsolute2relativeDist(double absoluteDist){
		String s = Double.toString(absoluteDist);
		if(absolute2relative.containsKey(s)){
			System.err.println( absoluteDist + " in absolute distance is equal to " + absolute2relative.get(s) + " in relative distance.");
			return absolute2relative.get(s);
		}
		e = new CoherenceBasedEEM(parent, parent.allGenes);
		e.setAbsoluteRadius(absoluteDist);
		try{
			e.findModuleGenes();
		}	
		catch (Exception err) {
			throw new MyException("Unable to convert distance!");
		}
		double relativeDist = ((double) e.getModuleGenes().size())/parent.allGenes.size();	
		System.err.println( absoluteDist + " in absolute distance is equal to " + relativeDist + " in relative distance.");
		absolute2relative.put(s, relativeDist);
		return relativeDist;
	}
	public double  convertRelative2absoluteDist(double  relativeDist){
		double upperAbsoluteDist;
		double lowerAbsoluteDist; 
		if( this.upperAbsoluteDist != null){
			upperAbsoluteDist  = this.upperAbsoluteDist;
		}else{
			upperAbsoluteDist = 2;
			//upperAbsoluteDist = parent.absoluteCorrelation?1:2;
		}
		if(this.lowerAbsoluteDist != null){
			lowerAbsoluteDist  = this.lowerAbsoluteDist;
		}else{
			lowerAbsoluteDist = 0;
		}
		System.err.println("converting " + relativeDist +  " in relative distance to absolute distance...");   
		double upperRelativeDist = convertAbsolute2relativeDist(upperAbsoluteDist);
		  if(upperRelativeDist < relativeDist){
		    throw new MyException("Unable to convert distance! increase upperAbsoluteDist(upperRelativeDist=" + upperRelativeDist + " relativeDist=" + relativeDist + ").");
		  }
		  double lowerRelativeDist = convertAbsolute2relativeDist(lowerAbsoluteDist);
		  if(lowerRelativeDist > relativeDist){
			  throw new MyException("Unable to convert distance! increase lowerAbsoluteDist(lowerRelativeDist=" + lowerRelativeDist + " relativeDist=" + relativeDist + ").");
		  }
		  int i,j;
		  j = 0;
		  Double midAbsoluteDist =  null;
		  for(i = 0; i < maxItrForBisection; i++){
			  midAbsoluteDist = (upperAbsoluteDist + lowerAbsoluteDist) / 2;
			  double midRelativeDist =  convertAbsolute2relativeDist(midAbsoluteDist);
			  if(midRelativeDist > upperRelativeDist || midRelativeDist < lowerRelativeDist){
				  j++;
				  if(j == maxItrForBisection2){
					  throw new MyException("Unable to convert distance! convertAbsolute2relativeDist might not be a monolonically increasing function.");
				  }
			  }
			  if(Math.abs(midRelativeDist - relativeDist) / relativeDist < deltaForBisection){
				  return midAbsoluteDist;
			  }
			  if(midRelativeDist < relativeDist){
				  lowerAbsoluteDist = midAbsoluteDist;
				  lowerRelativeDist = midRelativeDist;
			  }else{
				  upperAbsoluteDist = midAbsoluteDist;
				  upperRelativeDist = midRelativeDist;
			  }
		  }
		  //throw new MyException("Unable to convert distance! Increase maxItrForBisection or deltaForBisection."); 
		  System.err.println("WARN: Unable to convert distance! Increase maxItrForBisection or deltaForBisection.");
		  return midAbsoluteDist;
	}
}
