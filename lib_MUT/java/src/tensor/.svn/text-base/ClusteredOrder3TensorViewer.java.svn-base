package tensor;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import utility.*;

import processing.core.PApplet;

public class ClusteredOrder3TensorViewer  extends  Order3TensorViewer{

	private static final long serialVersionUID = -5287573654558352861L;
	
	
	protected int alpha = 180;
	
	protected ClusteredOrder3TensorWithAnnotation T;
	
	protected String outFile;
	protected boolean autoStop = false;
	
	protected float  profileWidth =  500;
	protected float  profileHeight = 500;
	
	protected int rowLabelFontSize = 10;
	protected int colLabelFontSize = 10;
	protected int sliceLabelFontSize = 20;
	protected float spaceBetweenRowLablelAndProfile = 10;
	protected float rowLabelSpaceWidth = 100;
	protected float spaceBetweenColLablelAndProfile = 10;
	protected float colLabelSpaceWidth = 100;
	
	protected int marginX = 20;
	protected int marginY = 20;
	
	
	protected float boxWidth = 0;
	protected float boxHeight = 0;
	
	
	protected int color1 = MyColor.getRGBhex("BLUE");
	protected int color2 = MyColor.getRGBhex("WHITE");
	protected int color3 = MyColor.getRGBhex("RED");
	
	protected int bgColor =  MyColor.getRGBhex("WHITE");
	
	protected float rowDendBase = 10; 
	protected float rowDendHeight = 90;
	protected float spaceNext2RowDend = 5;
	
	protected float colDendBase = 10;
	protected float colDendHeight = 90;
	protected float spaceNext2ColDend = 5;
	
	protected int rowDendLineWidth = 1;
	protected int colDendLineWidth = 1; 
	
	protected float rowAnnotBandwidth = 20;
	protected float spaceNext2RowAnnot = 5;
	protected float colAnnotBandwidth = 20;
	protected float spaceNext2ColAnnot = 5;
	
	protected int rowAnnotLabelFontSize = 10;
	protected int colAnnotLabelFontSize = 10;
	protected int indicaterFontSize = 10;
	
	protected int colorForAnnotBand1 = MyColor.getRGBhex("BLUE");
	protected int colorForAnnotBand2 = MyColor.getRGBhex("YELLOW");
	protected int colorForAnnotBand3 = MyColor.getRGBhex("RED");
	
	float max;
	float min;
	
	public void setAnnotBandColor(String  lowColor, String midColor, String highColor){
		colorForAnnotBand1 = MyColor.getRGBhex(lowColor);
		colorForAnnotBand2 = MyColor.getRGBhex(midColor);
		colorForAnnotBand3 = MyColor.getRGBhex(highColor);	
	}
	
	public void setAnnotBandColor(String  lowColor, String highColor){
		colorForAnnotBand1 = MyColor.getRGBhex(lowColor);
		colorForAnnotBand2 = PApplet.lerpColor(MyColor.getRGBhex(lowColor), MyColor.getRGBhex(highColor), (float) 0.5, RGB);
		colorForAnnotBand3 = MyColor.getRGBhex(highColor);	
	}
	
	protected StringMat rowAnnotation = null;
	protected StringMat colAnnotation = null;
	
	protected HierarchicalClustering rowClustering  = null;
	protected HierarchicalClustering colClustering  = null;
	
	public ClusteredOrder3TensorViewer(){
	}
	
	public ClusteredOrder3TensorViewer  (ClusteredOrder3TensorWithAnnotation T){
		this.T = T;
		setSlicedOder(3);
		calculateSize();
		max = (float) T.max();
		min = (float) T.min();
	}
	
	public void setSlicedOder(int i){
		boxWidth = 0;
		boxHeight = 0;
		if(i == 1){
			slicedOrder = 1;
			nrow = T.getDimOfOrder2();
			ncol = T.getDimOfOrder3();
			nslice = T.getDimOfOrder1();
			rowClustering = T.getOrder2Clustering();
			colClustering = T.getOrder3Clustering();
			rowAnnotation = T.getOrder2AnnotationMatrix();
			colAnnotation = T.getOrder3AnnotationMatrix();
		}
		if(i == 2){
			slicedOrder = 2;
			nrow = T.getDimOfOrder3();
			ncol = T.getDimOfOrder1();
			nslice = T.getDimOfOrder2();
			rowClustering = T.getOrder3Clustering();
			colClustering = T.getOrder1Clustering();
			rowAnnotation = T.getOrder3AnnotationMatrix();
			colAnnotation = T.getOrder1AnnotationMatrix();
		}
		if(i == 3){
			slicedOrder = 3;
			nrow = T.getDimOfOrder1();
			ncol = T.getDimOfOrder2();
			nslice = T.getDimOfOrder3();
			rowClustering = T.getOrder1Clustering();
			colClustering = T.getOrder2Clustering();
			rowAnnotation = T.getOrder1AnnotationMatrix();
			colAnnotation = T.getOrder2AnnotationMatrix();
			
		}
	}
	
	protected void calculateSize(){
		if(boxWidth != 0){
			profileWidth = boxWidth*ncol;	
		}
		if(boxHeight != 0){
			profileHeight = boxHeight*nrow;
		}
		
		profileX1 = marginX + (rowClustering!=null?rowDendHeight+rowDendBase+spaceNext2RowDend:0) 
			+ (rowAnnotation != null ?rowAnnotBandwidth*rowAnnotation.colSize() + spaceNext2RowAnnot:0); 
		profileX2 = profileX1 + profileWidth;
		Width  =  profileX2 + spaceBetweenRowLablelAndProfile + rowLabelSpaceWidth + marginX;
		boxWidth = (profileX2-profileX1)/ncol;
		
		profileY1 = marginY + (colClustering!=null?colDendHeight+colDendBase+spaceNext2ColDend:0)
			+ (colAnnotation != null ?colAnnotBandwidth*colAnnotation.colSize() + spaceNext2ColAnnot:0); 
		profileY2 = profileY1 + profileHeight;
		Height =  profileY2 + spaceBetweenColLablelAndProfile + colLabelSpaceWidth + marginY;
		boxHeight = (profileY2-profileY1)/nrow;
	}
	

	
	public void scaleColorByRow(){
		colorScalingByRow = true;
		colorScalingByCol = false;
	}
	
	public void scaleColorByColumn(){
		colorScalingByRow = false;
		colorScalingByCol = true;
	}
	

   protected void drawSliceLable(){
		textSize(sliceLabelFontSize);
		textAlign(RIGHT, BOTTOM);
		fill(0,alpha);
		text(getSliceLable(), (float)(marginX + rowDendHeight), (float)(marginY + colDendHeight));
	}
	
	protected String getSliceLable(){
		String label = "";
		if(slicedOrder == 1){
			label = T.getName1().get(sliceIndex);
		}else if(slicedOrder == 2){
			label = T.getName2().get(sliceIndex);
		}else if(slicedOrder == 3){
			label = T.getName3().get(sliceIndex);
		}
		return label;
	}
	
	
	protected MyMat getSlice(){
		if(slicedOrder == 1){
		return T.getOrder1Slice(sliceIndex);
		}else if(slicedOrder == 2){
			return T.getOrder2Slice(sliceIndex);
		}else if(slicedOrder == 3){
			return T.getOrder3Slice(sliceIndex);
		}else{
			return null;
		}
	}
	
	public void draw(){
		if(outFile != null){
			beginRecord(PDF, outFile + "." + getSliceLable() + ".pdf");
		}
		textFont(font);
		background(bgColor);
		drawProfile(getSlice());
		drawRowlabels(getSlice());
		drawCollabels(getSlice());
		drawIndicator(getSlice());
		if(rowClustering!=null){
			drawRowDendrogram();
		}
		if(colClustering!=null){
			drawColDendrogram();
		}
		if(rowAnnotation!=null){
			drawRowAnnotBand();
		}
		if(colAnnotation!=null){
			drawColAnnotBand();
		}
		drawSliceLable();
		if(outFile != null){
			endRecord();
			sliceIndex++;
		}
		if(autoStop & sliceIndex == nslice){
			this.stop();
		}
	}
	
	public void mousePressed(){
		sliceIndex++;
		if(sliceIndex == nslice){
			sliceIndex = 0;
		}	
	}
	
	protected  void drawIndicator(){
		drawIndicator(getSlice());
	}
	
	
	protected void drawIndicator(MyMat myMat){
		textSize(indicaterFontSize);
		textAlign(CENTER, BOTTOM);
		fill(0);
		double x = ((mouseX -profileX1)/boxWidth);
		double y = ((mouseY -profileY1)/boxHeight);
		int j = (int)x;
		int i = (int)y;
		if(x > 0 && y > 0 && i>=0 && i < myMat.rowSize() && j>=0 && j < myMat.colSize()){
			
			String label = myMat.getRowNames().get(i) + "\n" + myMat.getColNames().get(j) + "\n" +MyFunc.getEfficientRoundString(myMat.get(i, j),3);
			text(label, mouseX, mouseY);
		}
		if(rowAnnotation!=null){
			float rowAnnotX1 = marginX + (rowClustering!=null?rowDendHeight+rowDendBase+spaceNext2RowDend:0);
			float rowAnnotY1 = profileY1;
			x = ((mouseX-rowAnnotX1)/rowAnnotBandwidth);
			y = ((mouseY-rowAnnotY1)/boxHeight);
			j = (int)x;
			i = (int)y;
			StringMat A = rowAnnotation;
			if(x > 0 && y > 0 && i>=0 && i < myMat.rowSize() && j>=0 && j < A.colSize()){
				String label = A.getRowNames().get(i) + "\n" + A.getColNames().get(j) + "\n" + A.get(i, j);
				text(label, mouseX, mouseY);		
			}	
		}
		if(colAnnotation!=null){
			float colAnnotY1 = marginY + (colClustering!=null?colDendHeight+colDendBase+spaceNext2ColDend:0);
			float colAnnotX1 = profileX1;
			x = ((mouseX-colAnnotX1)/boxWidth);
			y = ((mouseY-colAnnotY1)/colAnnotBandwidth);
			i = (int)x;
			j = (int)y;
			StringMat A = colAnnotation;
			if(x > 0 && y > 0 && i>=0 && i < myMat.rowSize() && j>=0 && j < A.colSize()){
				String label = A.getRowNames().get(i) + "\n" + A.getColNames().get(j) + "\n" + A.get(i, j);
				text(label, mouseX, mouseY);		
			}
		}
			
	}
	
	protected  void drawProfile(){
		drawProfile(getSlice());
	}
	
	protected  void drawProfile(MyMat M){
		int ncol = M.colSize();
		int nrow = M.rowSize();
		//float max;
		//float min;
		if(colorScalingByRow){
			for(int i = 0; i < nrow; i++){
				List <Double> row = M.getRow(i);
				max = (float)MyFunc.max(row);
				min = (float)MyFunc.min(row);
				for(int j = 0; j < ncol; j++){
					float p = (float) ((M.get(i, j)-min)/(max-min));
					int color = MyColor.lerpColor(color1 ,color2,color3, p, RGB);
					fill(color,alpha);
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
					fill(color,alpha);
					noStroke();
					rect((profileX1+j*boxWidth),(profileY1+i*boxHeight),boxWidth,boxHeight);
				}
			}
			return;
		}
		//max = (float)MyFunc.max(M.asList());
		//min = (float)MyFunc.min(M.asList());
		//max = (float)MyFunc.percentile(M.asList(),0.99);
		//min = (float)MyFunc.percentile(M.asList(),0.01);
		for(int i = 0; i < nrow; i++){
			for(int j = 0; j < ncol; j++){
				float p;
				if(M.get(i, j)>=max){
					p=(float)1.0;
				}else if(M.get(i, j)<=min){
					p=(float)0.0;
				}else{
					p = (float) ((M.get(i, j)-min)/(max-min));
				}
				int color = MyColor.lerpColor(color1 ,color2,color3, p, RGB);
				fill(color,alpha);
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
		Map <String, Double> node2dist = rowClustering.getNode2distMap();
		Map <String, String[]> node2daughter = rowClustering.getNode2daughterMap();
		Map <String, String> node2parent = rowClustering.getNode2parentMap();
		List <String> rownames = rowClustering.getSortedTerminalNodes();
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
		
		Map <String, Double> node2dist = colClustering.getNode2distMap();
		Map <String, String[]> node2daughter = colClustering.getNode2daughterMap();
		Map <String, String> node2parent = colClustering.getNode2parentMap();
		List <String> colnames = colClustering.getSortedTerminalNodes();
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

	
	protected  void drawRowlabels(){
		drawRowlabels(getSlice());
	}
	
	protected void drawRowlabels(MyMat M){
		List <String> rowname = M.getRowNames();
		textSize(rowLabelFontSize);
		textAlign(LEFT, CENTER);
		fill(0,alpha);
		for(int i = 0; i < nrow; i++){
			text(rowname.get(i), profileX2+spaceBetweenRowLablelAndProfile, (float)(profileY1 + (i+0.5)*boxHeight));
		}
	}
	
	protected  void drawCollabels(){
		drawCollabels(getSlice());
	}

	protected void drawCollabels(MyMat M){
		List <String> colname = M.getColNames();
		textSize(colLabelFontSize);
		textAlign(RIGHT,CENTER);
		fill(0,alpha);
		for(int i = 0; i < ncol; i++){
			float x = (float)(profileX1 + (i+0.5)*boxWidth);
			float y =  profileY2+spaceBetweenColLablelAndProfile;
			rotatedText(colname.get(i),x,y, HALF_PI*3);
		}
	}
	
	protected void drawRowAnnotBand(){
		float rowAnnotX1 = marginX + (rowClustering!=null?rowDendHeight+rowDendBase+spaceNext2RowDend:0);
		float rowAnnotY1 = profileY1;
		List <String> annotTypes = rowAnnotation.getColNames();
		for(int i = 0; i < annotTypes.size(); i++){
			String annotType = annotTypes.get(i);
			List <String> annotString = rowAnnotation.getCol(annotType);
			List <String> annotStringNotEmpty = MyFunc.removeEmptyString(annotString);
			List <String> annotEntry = MyFunc.uniq(annotStringNotEmpty); 
			Collections.sort(annotEntry);
			if(annotEntry.size() == 1){
				for(int j = 0; j < annotString.size(); j++){
					if(!annotString.get(j).equals("")){
						fill(colorForAnnotBand3,alpha);
						noStroke();
						rect((rowAnnotX1+i*rowAnnotBandwidth),(rowAnnotY1+j*boxHeight),rowAnnotBandwidth,boxHeight);	
					}
				}
				fill(0,alpha);
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
						fill(color,alpha);
						noStroke();
						rect((rowAnnotX1+i*rowAnnotBandwidth),(rowAnnotY1+j*boxHeight),rowAnnotBandwidth,boxHeight);
					}else{
						fill(bgColor,alpha);
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
						fill(annot2color.get(annotString.get(j)), alpha);
						noStroke();
						rect((rowAnnotX1+i*rowAnnotBandwidth),(rowAnnotY1+j*boxHeight),rowAnnotBandwidth,boxHeight);						
					}else{
						fill(bgColor, alpha);
						noStroke();
						rect((rowAnnotX1+i*rowAnnotBandwidth),(rowAnnotY1+j*boxHeight),rowAnnotBandwidth,boxHeight);
					}
				}
			}
			fill(0,alpha);
			textSize(rowAnnotLabelFontSize);
			textAlign(RIGHT,CENTER);
			rotatedText(annotType, (float)(rowAnnotX1+(i+0.5)*rowAnnotBandwidth), (float)(profileY2 + spaceBetweenColLablelAndProfile),HALF_PI*3);
		}		
	}
	
	protected void drawColAnnotBand(){
		float colAnnotY1 = marginY + (colClustering!=null?colDendHeight+colDendBase+spaceNext2ColDend:0);
		float colAnnotX1 = profileX1;
		List <String> annotTypes = colAnnotation.getColNames();
		for(int i = 0; i < annotTypes.size(); i++){
			String annotType = annotTypes.get(i);
			List <String> annotString = colAnnotation.getCol(annotType);
			List <String> annotStringNotEmpty = MyFunc.removeEmptyString(annotString);
			List <String> annotEntry = MyFunc.uniq(annotStringNotEmpty); 
			Collections.sort(annotEntry);
			if(annotEntry.size()== 1){
				for(int j = 0; j < annotString.size(); j++){
					if(!annotString.get(j).equals("")){
						fill(colorForAnnotBand3,alpha);
						noStroke();
						rect((colAnnotX1+j*boxWidth),(colAnnotY1+i*colAnnotBandwidth),boxWidth,colAnnotBandwidth);		
					}
				}
				fill(0,alpha);
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
						fill(color,alpha);
						noStroke();
						rect((colAnnotX1+j*boxWidth),(colAnnotY1+i*colAnnotBandwidth),boxWidth,colAnnotBandwidth);
					}else{
						fill(bgColor,alpha);
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
						fill(annot2color.get(annotString.get(j)),alpha);
						noStroke();
						rect((colAnnotX1+j*boxWidth),(colAnnotY1+i*colAnnotBandwidth),boxWidth,colAnnotBandwidth);						
					}else{
						fill(bgColor,alpha);
						noStroke();
						rect((colAnnotX1+j*boxWidth),(colAnnotY1+i*colAnnotBandwidth),boxWidth,colAnnotBandwidth);
					}
				}
			}
			fill(0,alpha);
			textSize(colAnnotLabelFontSize);
			textAlign(LEFT, CENTER);
			text(annotType, (float)(profileX2 + spaceBetweenRowLablelAndProfile), (float)(colAnnotY1+(i+0.5)*colAnnotBandwidth));
		}		
	}
	
	protected void calculateRowLabelSpaceWidth(ClusteredMyMat M){
		textFont(font);
		textSize(rowLabelFontSize);
		float max  =  -(float)Double.MAX_VALUE;
		for(String s: M.getRowNames()){
			float tmp = textWidth(s);
			if(tmp  > max){
				max = tmp;
			}
		}
		if(colAnnotation!=null){
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
	protected void calculateColLabelSpaceWidth(ClusteredMyMat M){
		textFont(font);
		textSize(colLabelFontSize);
		float max  =  -(float)Double.MAX_VALUE;
		for(String s: M.getColNames()){
			float tmp = textWidth(s);
			if(tmp  > max){
				max = tmp;
			}
		}
		if(rowAnnotation!=null){
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
