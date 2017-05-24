package eem;

import java.util.*;

import processing.core.*;
import utility.*;

public class ExpressionModuleSetViewer extends MyMatViewer{
	private static final long serialVersionUID = 1946288376618794634L;
	int color = color(255,0,0);
	int lineColor = color(50);
	
	public ExpressionModuleSetViewer(ExpressionModuleSet expModSet){
		List <String> expModIds = expModSet.getIds(); 
		List <String> moduleGeneIds = expModSet.getModuleGeneUnion();
		
		M = new MyMat(expModIds, moduleGeneIds); 
		for(String s: expModIds){
			for(String t: moduleGeneIds){
				if(expModSet.get(s).getModuleGenes().contains(t)){
					M.set(s,t,1);
				}
			}
		}
		
		Map <String, Double> module2Score  = new LinkedHashMap<String,Double>();
		for(String s: moduleGeneIds){
			double tmp = 0;
			for(String t: expModIds){
				tmp += M.get(t, s)*(expModSet.get(t).getPvalue());
		
			}
			module2Score.put(s, tmp);
		}
		moduleGeneIds = MyFunc.sortKeysByDescendingOrderOfValues(module2Score);
		M.reorderCols(moduleGeneIds);
		ncol = M.colSize();
		nrow = M.rowSize();
		calculateSize();
	}	
	

	public void setup(){
		size(Math.round(Width), Math.round(Height));
		textFont(font);
		noLoop();
	}
	
	protected void drawProfile(){
		for(int i = 0; i < nrow; i++){
			for(int j = 0; j < ncol; j++){
				stroke(lineColor);
				if(M.get(i,j) == 0){
					noFill();
				}else{
					fill(color);
				}
				rect((profileX1+j*boxWidth),(profileY1+i*boxHeight),boxWidth,boxHeight);
			}
		}	
		
	}
	
	
	
	
}