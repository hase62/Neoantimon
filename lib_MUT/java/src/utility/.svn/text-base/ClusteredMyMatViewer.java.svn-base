package utility;
import java.awt.*;
import java.util.*;
import java.util.List;

import processing.core.*;


public class ClusteredMyMatViewer  extends MyMatViewer{
	private static final long serialVersionUID = 2457062611453720584L;

	protected ClusteredMyMat M;
	
	float rowDendBase = 10; 
	float rowDendHeight = 90;
	float spaceNext2RowDend = 5;
	float colDendBase = 10;
	float colDendHeight = 90;
	float spaceNext2ColDend = 5;
	
	int rowDendLineWidth = 1;
	int colDendLineWidth = 1; 
	
	
	float rowAnnotBandwidth = 20;
	float spaceNext2RowAnnot = 5;
	float colAnnotBandwidth = 20;
	float spaceNext2ColAnnot = 5;
	int rowAnnotLabelFontSize = 10;
	int colAnnotLabelFontSize = 10;
	
	
	
	int colorForAnnotBand1 = MyColor.getRGBhex("BLUE");
	int colorForAnnotBand2 = MyColor.getRGBhex("YELLOW");
	int colorForAnnotBand3 = MyColor.getRGBhex("RED");
	
	protected void setAnnotBandColor(String  lowColor, String midColor, String highColor){
		colorForAnnotBand1 = MyColor.getRGBhex(lowColor);
		colorForAnnotBand2 = MyColor.getRGBhex(midColor);
		colorForAnnotBand3 = MyColor.getRGBhex(highColor);	
	}
	
	public void setAnnotBandColor(String  lowColor, String highColor){
		colorForAnnotBand1 = MyColor.getRGBhex(lowColor);
		colorForAnnotBand2 = PApplet.lerpColor(MyColor.getRGBhex(lowColor), MyColor.getRGBhex(highColor), (float) 0.5, RGB);
		colorForAnnotBand3 = MyColor.getRGBhex(highColor);	
	}
	
	
	boolean hasRowAnnotation;
	boolean hasColAnnotation;
	
	boolean isRowClustered;
	boolean isColClustered;
	
	public ClusteredMyMatViewer(){};
	
	
	public ClusteredMyMatViewer (ClusteredMyMat M){
		this.M = M;
		ncol = M.colSize();
		nrow = M.rowSize();
		rowname = M.getRowNames();
		colname = M.getColNames();
		isRowClustered = M.isRowClustered();
		isColClustered = M.isColClustered();
		if(M instanceof ClusteredMyMatWithAnnotation){
			hasRowAnnotation = ((ClusteredMyMatWithAnnotation)M).hasRowAnnotation();
			hasColAnnotation = ((ClusteredMyMatWithAnnotation)M).hasColAnnotation();
		}
		calculateSize();
	}
	
	protected void calculateSize(){
		if(boxWidth != 0){
			profileWidth = boxWidth*ncol;	
		}
		if(boxHeight != 0){
			profileHeight = boxHeight*nrow;
		}
		
		profileX1 = marginX + (isRowClustered?rowDendHeight+rowDendBase+spaceNext2RowDend:0) 
			+ (hasRowAnnotation?rowAnnotBandwidth*((ClusteredMyMatWithAnnotation)M).getRowAnnotationTypes().size() + spaceNext2RowAnnot:0); 
		profileX2 = profileX1 + profileWidth;
		Width  =  profileX2 + spaceBetweenRowLablelAndProfile + rowLabelSpaceWidth + marginX;
		boxWidth = (profileX2-profileX1)/ncol;
		
		profileY1 = marginY + (isColClustered?colDendHeight+colDendBase+spaceNext2ColDend:0)
			+ (hasColAnnotation?colAnnotBandwidth*((ClusteredMyMatWithAnnotation)M).getColAnnotationTypes().size() + spaceNext2ColAnnot:0); 
		profileY2 = profileY1 + profileHeight;
		Height =  profileY2 + spaceBetweenColLablelAndProfile + colLabelSpaceWidth + marginY;
		boxHeight = (profileY2-profileY1)/nrow;
	}
	
	
		
	public void draw(){
		if(outFile != null){
			beginRecord(PDF, outFile);
		}
		textFont(font);
		background(bgColor);
		drawProfile();
		if(isRowClustered){
			drawRowDendrogram();
		}
		if(isColClustered){
			drawColDendrogram();
		}
		if(hasRowAnnotation){
			drawRowAnnotBand();
		}
		if(hasColAnnotation){
			drawColAnnotBand();
		}
		drawRowlabels();
		drawCollabels();
		drawIndicator();
		if(outFile != null){
			endRecord();
		}
		if(autoStop){
			this.stop();
			System.exit(0);
		}
	}
	
	
	

	protected void drawIndicator(){
		textSize(indicaterFontSize);
		textAlign(CENTER, BOTTOM);
		fill(0);
		double x = ((mouseX -profileX1)/boxWidth);
		double y = ((mouseY -profileY1)/boxHeight);
		int j = (int)x;
		int i = (int)y;
		if(x > 0 && y > 0 && i>=0 && i < nrow && j>=0 && j < ncol){
			String label = rowname.get(i) + "\n" + colname.get(j) + "\n" + M.get(i, j);
			text(label, mouseX, mouseY);
		}
		if(hasRowAnnotation){
			float rowAnnotX1 = marginX + (isRowClustered?rowDendHeight+rowDendBase+spaceNext2RowDend:0);
			float rowAnnotY1 = profileY1;
			x = ((mouseX-rowAnnotX1)/rowAnnotBandwidth);
			y = ((mouseY-rowAnnotY1)/boxHeight);
			j = (int)x;
			i = (int)y;
			StringMat A = ((ClusteredMyMatWithAnnotation)M).getRowAnnotationMatrix();
			if(x > 0 && y > 0 && i>=0 && i < nrow && j>=0 && j < A.ncol){
				String label = A.rowname.get(i) + "\n" + A.colname.get(j) + "\n" + A.get(i, j);
				text(label, mouseX, mouseY);		
			}	
		}
		if(hasColAnnotation){
			float colAnnotY1 = marginY + (isColClustered?colDendHeight+colDendBase+spaceNext2ColDend:0);
			float colAnnotX1 = profileX1;
			x = ((mouseX-colAnnotX1)/boxWidth);
			y = ((mouseY-colAnnotY1)/colAnnotBandwidth);
			i = (int)x;
			j = (int)y;
			StringMat A = ((ClusteredMyMatWithAnnotation)M).getColAnnotationMatrix();
			if(x > 0 && y > 0 && i>=0 && i < ncol && j>=0 && j < A.ncol){
				String label = A.rowname.get(i) + "\n" + A.colname.get(j) + "\n" + A.get(i, j);
				text(label, mouseX, mouseY);		
			}
		}
			
	}
	
	
	protected  void drawProfile(){
		if(colorScalingByRow){
			for(int i = 0; i < nrow; i++){
				List <Double> row = M.getRow(i);
				max = (float)MyFunc.max(row);
				min = (float)MyFunc.min(row);
				for(int j = 0; j < ncol; j++){
					float p = (float) ((M.get(i, j)-min)/(max-min));
					int color = MyColor.lerpColor(color1 ,color2,color3, p, RGB);
					fill(color);
					noStroke();
					rect((profileX1+j*boxWidth),(profileY1+i*boxHeight),boxWidth,boxHeight);
				}
			}
			return;
		}
		if(colorScalingByCol){
			for(int j = 0; j < ncol; j++){
				List <Double> col = M.getCol(j);
				max = (float)MyFunc.max(col);
				min = (float)MyFunc.min(col);
				for(int i = 0; i < nrow; i++){
					float p = (float) ((M.get(i, j)-min)/(max-min));
					int color = MyColor.lerpColor(color1 ,color2,color3, p, RGB);
					fill(color);
					noStroke();
					rect((profileX1+j*boxWidth),(profileY1+i*boxHeight),boxWidth,boxHeight);
				}
			}
			return;
		}
		max = (float)MyFunc.max(M.asList());
		min = (float)MyFunc.min(M.asList());
		for(int i = 0; i < nrow; i++){
			for(int j = 0; j < ncol; j++){
				float p = (float) ((M.get(i, j)-min)/(max-min));
				int color = MyColor.lerpColor(color1 ,color2,color3, p, RGB);
				fill(color);
				noStroke();
				rect((profileX1+j*boxWidth),(profileY1+i*boxHeight),boxWidth,boxHeight);
			}
		}
	}	
		
	private void drawRowDendrogram(){
		float rowDendX1 =  marginX;
		float rowDendX2 =  marginX + rowDendHeight+rowDendBase;
		float rowDendY1 = profileY1;
		stroke(color(0));
		strokeWeight(rowDendLineWidth);
		Map <String, Double> node2dist = M.getRowClustering().getNode2distMap();
		Map <String, String[]> node2daughter = M.getRowClustering().getNode2daughterMap();
		Map <String, String> node2parent = M.getRowClustering().getNode2parentMap();
		List <String> rownames = M.getRowNames();
		Map <String, Double> node2Y = new  HashMap<String, Double>();
		Map <String, Double> node2X1 = new  HashMap<String, Double>();
		Map <String, Double> node2X2 = new  HashMap<String, Double>();
		for(int i = 0; i < nrow; i++){
			double Y  =  (rowDendY1 + (i+0.5)*boxHeight);
			double X1 =  (rowDendX2 - rowDendBase - rowDendHeight * node2dist.get(node2parent.get(rownames.get(i))));
			double X2 = (double)rowDendX2;
			node2Y.put(rownames.get(i), Y);
			node2X1.put(rownames.get(i), X1);
			node2X2.put(rownames.get(i), X2);
			line((float)X1,(float)Y,(float)X2,(float)Y);
		}
		for(int i  = 1; i < nrow-1; i++){
			String node = "node" + i;
			double Y = ( node2Y.get(node2daughter.get(node)[0]) + node2Y.get(node2daughter.get(node)[1]) )/2;
			double X1 = (rowDendX2 - rowDendBase - rowDendHeight * node2dist.get(node2parent.get(node)));
			double X2 = node2X1.get(node2daughter.get(node)[0]);
			node2Y.put(node, Y);
			node2X1.put(node, X1);
			node2X2.put(node, X2);
			line((float)X1,(float)Y,(float)X2,(float)Y);
			line((float)X2,(float)(double)node2Y.get(node2daughter.get(node)[0]),(float)X2,(float)(double)node2Y.get(node2daughter.get(node)[1]));
		}
		line(rowDendX1, (float)(double)node2Y.get(node2daughter.get("node" + (nrow-1))[0]) , rowDendX1, (float)(double)node2Y.get(node2daughter.get("node" + (nrow-1))[1]));
	}
		
	private void drawColDendrogram(){
		float colDendX1 =  profileX1;
		float colDendY1 = marginY;
		float colDendY2 =  marginY + colDendHeight+ colDendBase;
		
		stroke(color(0));
		strokeWeight(colDendLineWidth);
		
		Map <String, Double> node2dist = M.getColClustering().getNode2distMap();
		Map <String, String[]> node2daughter = M.getColClustering().getNode2daughterMap();
		Map <String, String> node2parent = M.getColClustering().getNode2parentMap();
		List <String> colnames = M.getColNames();
		Map <String, Double> node2X = new  HashMap<String, Double>();
		Map <String, Double> node2Y1 = new  HashMap<String, Double>();
		Map <String, Double> node2Y2 = new  HashMap<String, Double>();
		
		for(int i = 0; i < ncol; i++){
			double X  =  (colDendX1 + (i+0.5)*boxWidth);
			double Y1 =  (colDendY2 - colDendBase - colDendHeight * node2dist.get(node2parent.get(colnames.get(i))));
			double Y2 = (double)colDendY2;
			node2X.put(colnames.get(i), X);
			node2Y1.put(colnames.get(i), Y1);
			node2Y2.put(colnames.get(i), Y2);
			line((float)X,(float)Y1,(float)X,(float)Y2);
		}
		for(int i  = 1; i < ncol-1; i++){
			String node = "node" + i;
			double X = ( node2X.get(node2daughter.get(node)[0]) + node2X.get(node2daughter.get(node)[1]) )/2;
			double Y1 = (colDendY2 - colDendBase - colDendHeight * node2dist.get(node2parent.get(node)));
			double Y2 = node2Y1.get(node2daughter.get(node)[0]);
			node2X.put(node, X);
			node2Y1.put(node, Y1);
			node2Y2.put(node, Y2);
			line((float)X,(float)Y1,(float)X,(float)Y2);
			line((float)(double)node2X.get(node2daughter.get(node)[0]),(float)Y2,(float)(double)node2X.get(node2daughter.get(node)[1]),(float)Y2);
		}
		line((float)(double)node2X.get(node2daughter.get("node" + (ncol-1))[0]), colDendY1, (float)(double)node2X.get(node2daughter.get("node" + (ncol-1))[1]), colDendY1);
	}

	protected void drawRowlabels(){
		textSize(rowLabelFontSize);
		textAlign(LEFT, CENTER);
		fill(0);
		for(int i = 0; i < nrow; i++){
			text(rowname.get(i), profileX2+spaceBetweenRowLablelAndProfile, (float)(profileY1 + (i+0.5)*boxHeight));
		}
	}

	protected void drawCollabels(){
		textSize(colLabelFontSize);
		textAlign(RIGHT,CENTER);
		fill(0);
		for(int i = 0; i < ncol; i++){
			float x = (float)(profileX1 + (i+0.5)*boxWidth);
			float y =  profileY2+spaceBetweenColLablelAndProfile;
			rotatedText(colname.get(i),x,y, HALF_PI*3);
		}
	}
	
	
	
	
	private void drawRowAnnotBand(){
		float rowAnnotX1 = marginX + (isRowClustered?rowDendHeight+rowDendBase+spaceNext2RowDend:0);
		float rowAnnotY1 = profileY1;
		List <String> annotTypes = ((ClusteredMyMatWithAnnotation)M).getRowAnnotationTypes();
		for(int i = 0; i < annotTypes.size(); i++){
			String annotType = annotTypes.get(i);
			List <String> annotString = ((ClusteredMyMatWithAnnotation)M).getRowAnnotation(annotType);
			List <String> annotStringNotEmpty = MyFunc.removeEmptyString(annotString);
			List <String> annotEntry = MyFunc.uniq(annotStringNotEmpty); 
			Collections.sort(annotEntry);
			if(annotEntry.size() == 1){
				for(int j = 0; j < annotString.size(); j++){
					if(!annotString.get(j).equals("")){
						fill(colorForAnnotBand3);
						noStroke();
						rect((rowAnnotX1+i*rowAnnotBandwidth),(rowAnnotY1+j*boxHeight),rowAnnotBandwidth,boxHeight);	
					}
				}
				fill(0);
				textSize(rowAnnotLabelFontSize);
				textAlign(RIGHT,CENTER);
				rotatedText(annotType, (float)(rowAnnotX1+(i+0.5)*rowAnnotBandwidth), (float)(profileY2 + spaceBetweenColLablelAndProfile),HALF_PI*3);
				continue;
			}
			if(MyFunc.canBeDouble(annotStringNotEmpty)){
				List <Double> annotDouble = MyFunc.toDouble(annotStringNotEmpty);
				float max = (float)MyFunc.max(annotDouble);
				float min = (float)MyFunc.min(annotDouble);
				for(int j = 0; j < annotString.size(); j++){
					if(!annotString.get(j).equals("")){
						float p = (float)((Double.valueOf(annotString.get(j))-min)/(max-min));
						int color = MyColor.lerpColor(colorForAnnotBand1,colorForAnnotBand2,colorForAnnotBand3, p, HSB);
						fill(color);
						noStroke();
						rect((rowAnnotX1+i*rowAnnotBandwidth),(rowAnnotY1+j*boxHeight),rowAnnotBandwidth,boxHeight);
					}else{
						fill(bgColor);
						noStroke();
						rect((rowAnnotX1+i*rowAnnotBandwidth),(rowAnnotY1+j*boxHeight),rowAnnotBandwidth,boxHeight);
					}
				}
			}else{
				float d = (float)(1.0/(annotEntry.size()-1));
				Map <String, Integer> annot2color = new HashMap<String, Integer>();
				
				for(int j = 0 ; j < annotEntry.size(); j++){
					annot2color.put(annotEntry.get(j),(Integer)MyColor.lerpColor(colorForAnnotBand1,colorForAnnotBand2,colorForAnnotBand3, d*j, HSB));
				}
				
				for(int j = 0; j < annotString.size(); j++){
					if(!annotString.get(j).equals("")){
						fill(annot2color.get(annotString.get(j)));
						noStroke();
						rect((rowAnnotX1+i*rowAnnotBandwidth),(rowAnnotY1+j*boxHeight),rowAnnotBandwidth,boxHeight);						
					}else{
						fill(bgColor);
						noStroke();
						rect((rowAnnotX1+i*rowAnnotBandwidth),(rowAnnotY1+j*boxHeight),rowAnnotBandwidth,boxHeight);
					}
				}
			}
			fill(0);
			textSize(rowAnnotLabelFontSize);
			textAlign(RIGHT,CENTER);
			rotatedText(annotType, (float)(rowAnnotX1+(i+0.5)*rowAnnotBandwidth), (float)(profileY2 + spaceBetweenColLablelAndProfile),HALF_PI*3);
		}		
	}
	
	
	
	private void drawColAnnotBand(){
		float colAnnotY1 = marginY + (isColClustered?colDendHeight+colDendBase+spaceNext2ColDend:0);
		float colAnnotX1 = profileX1;
		List <String> annotTypes = ((ClusteredMyMatWithAnnotation)M).getColAnnotationTypes();
		for(int i = 0; i < annotTypes.size(); i++){
			String annotType = annotTypes.get(i);
			List <String> annotString = ((ClusteredMyMatWithAnnotation)M).getColAnnotation(annotType);
			List <String> annotStringNotEmpty = MyFunc.removeEmptyString(annotString);
			List <String> annotEntry = MyFunc.uniq(annotStringNotEmpty); 
			Collections.sort(annotEntry);
			if(annotEntry.size()== 1){
				for(int j = 0; j < annotString.size(); j++){
					if(!annotString.get(j).equals("")){
						fill(colorForAnnotBand3);
						noStroke();
						rect((colAnnotX1+j*boxWidth),(colAnnotY1+i*colAnnotBandwidth),boxWidth,colAnnotBandwidth);		
					}
				}
				fill(0);
				textSize(colAnnotLabelFontSize);
				textAlign(LEFT, CENTER);
				text(annotType, (float)(profileX2 + spaceBetweenRowLablelAndProfile), (float)(colAnnotY1+(i+0.5)*colAnnotBandwidth));
				continue;
			}
			if(MyFunc.canBeDouble(annotStringNotEmpty)){
				List <Double> annotDouble = MyFunc.toDouble(annotStringNotEmpty);
				float max = (float)MyFunc.max(annotDouble);
				float min = (float)MyFunc.min(annotDouble);
				for(int j = 0; j < annotString.size(); j++){
					if(!annotString.get(j).equals("")){
						float p = (float)((Double.valueOf(annotString.get(j))-min)/(max-min));
						int color = MyColor.lerpColor(colorForAnnotBand1,colorForAnnotBand2,colorForAnnotBand3, p, HSB);
						fill(color);
						noStroke();
						rect((colAnnotX1+j*boxWidth),(colAnnotY1+i*colAnnotBandwidth),boxWidth,colAnnotBandwidth);
					}else{
						fill(bgColor);
						noStroke();
						rect((colAnnotX1+j*boxWidth),(colAnnotY1+i*colAnnotBandwidth),boxWidth,colAnnotBandwidth);
					}
				}
			}else{
			
				float d = (float)(1.0/(annotEntry.size()-1));
				Map <String, Integer> annot2color = new HashMap<String, Integer>();
				for(int j = 0 ; j < annotEntry.size(); j++){
					annot2color.put(annotEntry.get(j),(Integer)MyColor.lerpColor(colorForAnnotBand1,colorForAnnotBand2,colorForAnnotBand3, d*j, HSB));
				}
				for(int j = 0; j < annotString.size(); j++){
					if(!annotString.get(j).equals("")){
						fill(annot2color.get(annotString.get(j)));
						noStroke();
						rect((colAnnotX1+j*boxWidth),(colAnnotY1+i*colAnnotBandwidth),boxWidth,colAnnotBandwidth);						
					}else{
						fill(bgColor);
						noStroke();
						rect((colAnnotX1+j*boxWidth),(colAnnotY1+i*colAnnotBandwidth),boxWidth,colAnnotBandwidth);
					}
				}
			}
			fill(0);
			textSize(colAnnotLabelFontSize);
			textAlign(LEFT, CENTER);
			text(annotType, (float)(profileX2 + spaceBetweenRowLablelAndProfile), (float)(colAnnotY1+(i+0.5)*colAnnotBandwidth));
		}		
	}
	
	protected void calculateRowLabelSpaceWidth(){
		textFont(font);
		textSize(rowLabelFontSize);
		float max  =  -(float)Double.MAX_VALUE;
		for(String s: M.getRowNames()){
			float tmp = textWidth(s);
			if(tmp  > max){
				max = tmp;
			}
		}
		if(hasColAnnotation){
			textSize(colAnnotLabelFontSize);
			for(String s: ((ClusteredMyMatWithAnnotation)M).getColAnnotationTypes()){
				float tmp = textWidth(s);
				if(tmp  > max){
					max = tmp;
				}	
			}
		}
		
		rowLabelSpaceWidth = max;
	}
	protected void calculateColLabelSpaceWidth(){
		textFont(font);
		textSize(colLabelFontSize);
		float max  =  -(float)Double.MAX_VALUE;
		for(String s: M.getColNames()){
			float tmp = textWidth(s);
			if(tmp  > max){
				max = tmp;
			}
		}
		if(hasRowAnnotation){
			textSize(rowAnnotLabelFontSize);
			for(String s: ((ClusteredMyMatWithAnnotation)M).getRowAnnotationTypes()){
				float tmp = textWidth(s);
				if(tmp  > max){
					max = tmp;
				}	
			}
		}
		colLabelSpaceWidth = max;
	}
	
	                                 
	
	
}
	
	

