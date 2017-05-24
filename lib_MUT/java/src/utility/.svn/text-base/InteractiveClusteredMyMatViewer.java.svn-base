package utility;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import processing.core.*;
import processing.pdf.PGraphicsPDF;
import controlP5.Button;
import controlP5.ControlEvent;
import controlP5.ControlP5;
import controlP5.Textfield;

public class InteractiveClusteredMyMatViewer  extends PApplet {
	
	protected float  profileWidth =  500;
	protected float  profileHeight = 500;
	
	protected int rowLabelFontSize = 10;
	protected int colLabelFontSize = 10;
	protected int indicaterFontSize = 10;
	protected float spaceBetweenRowLablelAndProfile = 10;
	protected float rowLabelSpaceWidth = 100;
	protected float spaceBetweenColLablelAndProfile = 10;
	protected float colLabelSpaceWidth = 100;
	
	protected int marginX = 20;
	protected int marginY = 50;
	
	protected float boxWidth = 0;
	protected float boxHeight = 0;
	
	
	protected int color1 = MyColor.getRGBhex("BLUE");
	protected int color2 = MyColor.getRGBhex("WHITE");
	protected int color3 = MyColor.getRGBhex("RED");
	
	protected int bgColor =  MyColor.getRGBhex("WHITE");
	
	
	protected ClusteredMyMatWithAnnotation originalM;  
	protected ClusteredMyMatWithAnnotation M;  
	protected ClusteredMyMatWithAnnotation preM;  
	
	protected float profileX1;
	protected float profileX2;
	protected float profileY1;
	protected float profileY2;

	protected float Width;
	protected float Height;
	
	protected PFont font = createFont("Bitstream Vera Sans",20);
	
	protected boolean colorScalingByRow = false;
	protected boolean colorScalingByCol = false;
	protected float max;
	protected float min;
	
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
	
	protected int colorForAnnotBand1 = MyColor.getRGBhex("BLUE");
	protected int colorForAnnotBand2 = MyColor.getRGBhex("YELLOW");
	protected int colorForAnnotBand3 = MyColor.getRGBhex("RED");
	
	
	public InteractiveClusteredMyMatViewer(){};
	
	public InteractiveClusteredMyMatViewer(ClusteredMyMatWithAnnotation M) {
		this.M = new ClusteredMyMatWithAnnotation(M);
		originalM = new ClusteredMyMatWithAnnotation(M);
		calculateSize();
	}

	protected int alpha = 180;
	
	protected int colorForHighlight = MyColor.getRGBhex("RED");
	
	protected String pdfFile = null;
	protected String labelFile = null;
	

	protected ControlP5 controlP5;
	
	
	
	public void setPdfFileName(String fileName){
		pdfFile = MyFunc.getPreffix(fileName);
	}
	
	public void setLabelFileName(String fileName){
		labelFile = MyFunc.getPreffix(fileName);
	}
	
	protected Set <String> highlightedRownameSet = null;
	protected Set <String> highlightedColnameSet = null;
	protected List <String> highlightedRowname = null;
	protected List <String> highlightedColname = null;
	protected String highlightedRowDendNode = null;
	protected Set <String> highlightedRowDendNodeDownstream = null;
	protected String highlightedColDendNode = null;
	protected Set <String> highlightedColDendNodeDownstream = null;
	
	
	protected boolean lockColHighlight = false;
	protected boolean lockRowHighlight = false;
	
	int colorActiveForControlP5 = color(255,182,193);
	int colorBgForControlP5 = color(200,100);
	int colorFgForControlP5 = color(255,182,193);
	
	boolean annotIndicatorDrawn = false;
	boolean refresh = false;
	
	
	protected void setButtons(){
		Button b = controlP5.addButton("Zoom",0,(int)marginX,10,30,20);
		b.setId(2);
		Button b2 = controlP5.addButton("Label",0,(int)marginX+40,10,30,20);
		b2.setId(3);
		Button b3 = controlP5.addButton("Pdf",0,(int)marginX+80,10,30,20);
		b3.setId(4);
		Button b4 = controlP5.addButton("Find",0,(int)marginX+120,10,30,20);
		b4.setId(5);
	}
	
	Textfield t;
	protected void setTextfield(){
		t = controlP5.addTextfield("search",(int)marginX+180,10,100,20);
		t.setFocus(true);
	}
	
	boolean zoomPressed = false;
	int zooming = 0;
	Integrator profileX1forZoom;
	Integrator profileX2forZoom;
	Integrator profileY1forZoom;
	Integrator profileY2forZoom;

	String highlightedLable = "";
	boolean eventHappened = false;
	
	public void controlEvent(ControlEvent theEvent) {  
		 if(theEvent.id()==2){
			 if(zoomPressed == false){
				 zoomPressed = true;
			 }
		 }else if(theEvent.id()==3){
			 if(!t.getText().equals("")){
				 setLabelFileName(t.getText());
				 t.clear();
			 }
		 }else if(theEvent.id()==4){
			 if(!t.getText().equals("")){
				 setPdfFileName(t.getText());
				 t.clear();
			 }
		 }else if(theEvent.id()==5){
			 if(!t.getText().equals("")){
				 highlightedLable = t.getText();
				 t.clear();
			 }else{
				 highlightedLable = ""; 
			 }
		 }
		 eventHappened = true;
	}	 


	private void printList(List<?> L, String outfile){
		try {
			outfile = outfile +  ".txt";
		PrintWriter os = new PrintWriter(new FileWriter(outfile));
		for(int i = 0, n = L.size(); i < n; i++){
			os.println(L.get(i).toString());
		}
		os.flush();
		os.close();
		System.err.println("print out the labels to a txt file.");
		}catch (Exception e) {
		}
	}
	
	protected void printLables(){
		if(highlightedColname != null){
			 if(!highlightedColname.isEmpty()){
				printList(highlightedColname,labelFile);
			 }
		 }
		if(highlightedRowname != null){
			 if(!highlightedRowname.isEmpty()){
				 printList(highlightedRowname,labelFile);
			 }
		 }
	}
	
	public void setupControlP5(){
		controlP5 = new ControlP5(this);
		controlP5.setColorBackground(colorBgForControlP5);
		controlP5.setColorForeground(colorFgForControlP5);
		controlP5.setColorActive(colorActiveForControlP5);
		setButtons();
		setTextfield();
	}
	
	
	public void setup(){
		size(Math.round(Width), Math.round(Height));
		setupControlP5();
		background(bgColor);
	}
	
	
	
	public void draw(){
		
		textFont(font);
		if(refresh == false & annotIndicatorDrawn == false & M.nrow * M.ncol > 2500 & eventHappened == false   & frameCount > 2){
			drawWithoutProfile();
		}else{
			annotIndicatorDrawn = false;
			refresh = false;
			background(bgColor);
			if(zoomPressed){
				zoom();	
				zoomPressed = false;
			}
			if(zooming != 0){
				if(zooming > 0){
					drawProfile(M, profileX1forZoom.value(), profileX2forZoom.value(), profileY1forZoom.value(), profileY2forZoom.value());
				}else{
					drawProfile(preM,profileX1forZoom.value(), profileX2forZoom.value(), profileY1forZoom.value(), profileY2forZoom.value());
				}
				boolean tmp1 = profileX1forZoom.update(); 
				boolean tmp2 = profileX2forZoom.update(); 
				boolean tmp3 = profileY1forZoom.update(); 
				boolean tmp4 = profileY2forZoom.update(); 
				if(!tmp1 &&  !tmp2 && !tmp3 && !tmp4){
					zooming = 0;
				}
			}else{
				if(pdfFile != null){
					beginRecord(PDF, pdfFile +  ".pdf");
					textFont(font);
				}
				drawProfile();
				
				if(M.getRowClustering()!=null){
					drawRowDendrogram();
				}
				if(M.getColClustering()!=null){
					drawColDendrogram();
				}
				drawRowlabels();
				drawCollabels();
				if(M.getRowAnnotationMatrix()!=null){
					drawRowAnnotBand();
					drawRowAnnotIndicator();
				}
				if(M.getColAnnotationMatrix()!=null){
					drawColAnnotBand();
					drawColAnnotIndicator();
				}
				if(M.nrow * M.ncol <= 2500){
					drawProfileIndicator();
				}
				if(pdfFile == null){
					controlP5.draw();
				}else{
					endRecord();
					pdfFile = null;
					System.err.println("print out the image to a pdf file.");
				}
				if(labelFile != null){
					printLables();
					labelFile = null;
				}
				eventHappened = false;
			}
		}
	}
	
	protected void drawWithoutProfile(){
		
		fill(255);
		rect(0, 0, Width,profileY1);
		rect(0, profileY1, profileX1,profileHeight);
		rect(profileX2, profileY1, Width-profileX2,profileHeight);
		rect(0, profileY2, Width,Height-profileY2);
		if(M.getRowClustering()!=null){
			drawRowDendrogram();
		}
		if(M.getColClustering()!=null){
			drawColDendrogram();
		}
		drawRowlabels();
		drawCollabels();
		if(M.getRowAnnotationMatrix()!=null){
			drawRowAnnotBand();
			drawRowAnnotIndicator();
		}
		if(M.getColAnnotationMatrix()!=null){
			drawColAnnotBand();
			drawColAnnotIndicator();
		}
	}
	
	public void mousePressed(){
		refresh = true;
		if((mouseX -  marginX)*(mouseX - (marginX + rowDendHeight+rowDendBase))<= 0 & (mouseY -profileY1)*(mouseY-profileY2)<=0){
			lockRowHighlight = (lockRowHighlight == true)?false:true;
		}
		if((mouseY -  marginY)*(mouseY - (marginY + colDendHeight+colDendBase))<= 0 & (mouseX -profileX1)*(mouseX-profileX2)<=0){
			lockColHighlight = (lockColHighlight == true)?false:true;
		}	
	}
	
	private void zoom(){
		if(highlightedRownameSet != null && !highlightedRownameSet.isEmpty()){
			zoomIn();	
		}else if(highlightedColnameSet != null && !highlightedColnameSet.isEmpty()){
			zoomIn();	
		}else{
			zoomOut();
		}
	}
	
	private void zoomIn(){	

		preM = new ClusteredMyMatWithAnnotation(M);
		
		profileX1forZoom = new  Integrator(profileX1);
		profileX2forZoom = new Integrator(profileX2);
		profileY1forZoom = new Integrator(profileY1);
		profileY2forZoom = new Integrator(profileY2);
		if(highlightedRownameSet != null){
			if(!highlightedRownameSet.isEmpty()){
				profileY1forZoom.set(profileY1+M.rowIndexOf(highlightedRowname.get(0))*boxHeight);
				profileY2forZoom.set(profileY1+(M.rowIndexOf(highlightedRowname.get(highlightedRowname.size()-1))+1)*boxHeight);
			}
		}
		if(highlightedColnameSet != null){
			if(!highlightedColnameSet.isEmpty()){
				profileX1forZoom.set(profileX1+M.colIndexOf(highlightedColname.get(0))*boxWidth);
				profileX2forZoom.set(profileX1+(M.colIndexOf(highlightedColname.get(highlightedColname.size()-1))+1)*boxWidth);
			}
		}

		profileX1forZoom.target(profileX1);
		profileX2forZoom.target(profileX2);
		profileY1forZoom.target(profileY1);
		profileY2forZoom.target(profileY2);
		
		zooming = 1;
		
		
		if(highlightedRownameSet != null){
			if(!highlightedRownameSet.isEmpty()){
				HierarchicalClustering H = M.getRowClustering().getSubHierarchicalClustering(highlightedRowDendNode);
				M.setRowClustering(H);				
				}
		}
			if(highlightedColnameSet != null){
				if(!highlightedColnameSet.isEmpty()){
					if(!highlightedColnameSet.isEmpty()){
						HierarchicalClustering H = M.getColClustering().getSubHierarchicalClustering(highlightedColDendNode);
						M.setColClustering(H);				
					}
				}
			}	
			
		
		
		
		calculateSize();
		
		highlightedRowDendNode = null;
		highlightedRowDendNodeDownstream = null;
		highlightedRownameSet =  null;
		highlightedRowname =  null;
		highlightedColDendNode =  null;
		highlightedColDendNodeDownstream =  null;
		highlightedColnameSet =  null;
		highlightedColname =  null;
		lockColHighlight = false;
		lockRowHighlight = false;
	
	}
	
	private void zoomOut(){
		
		List <String> rowname = M.getRowNames() ;
		List <String> colname = M.getColNames() ;
		preM = new ClusteredMyMatWithAnnotation(M);
		M = new ClusteredMyMatWithAnnotation(originalM);
		calculateSize();
		profileX1forZoom = new  Integrator(profileX1);
		profileX2forZoom = new Integrator(profileX2);
		profileY1forZoom = new Integrator(profileY1);
		profileY2forZoom = new Integrator(profileY2);
		profileX1forZoom.target(profileX1);
		profileX2forZoom.target(profileX2);
		profileY1forZoom.target(profileY1);
		profileY2forZoom.target(profileY2);
		
		if(rowname.size() != M.getRowNames().size()){
			profileY1forZoom.target(profileY1+M.rowIndexOf(rowname.get(0))*boxHeight);
			profileY2forZoom.target(profileY1+(M.rowIndexOf(rowname.get(rowname.size()-1))+1)*boxHeight);
		}
		if(colname.size() != M.getColNames().size()){
			profileX1forZoom.target(profileX1+M.colIndexOf(colname.get(0))*boxWidth);
			profileX2forZoom.target(profileX1+(M.colIndexOf(colname.get(colname.size()-1))+1)*boxWidth);
		}
		zooming = -1;
		
	}

	
	private void drawRowDendrogram(){
		drawRowDendrogram(M.getRowClustering());
	}	
	
	private void drawRowDendrogram(HierarchicalClustering H){
		float rowDendX1 =  marginX;
		float rowDendX2 =  marginX + rowDendHeight+rowDendBase;
		float rowDendY1 = profileY1;
		Map <String, Double> node2dist = H.getNode2distMap();
		Map <String, String[]> node2daughter = H.getNode2daughterMap();
		Map <String, String> node2parent = H.getNode2parentMap();
		List <String> rownames = H.getSortedTerminalNodes();
		int nrow = rownames.size();
		Map <String, Double> node2Y = new  HashMap<String, Double>();
		Map <String, Double> node2X1 = new  HashMap<String, Double>();
		Map <String, Double> node2X2 = new  HashMap<String, Double>();
		
		Map <String, List<Double>> node2Xline = new  HashMap<String, List<Double>>();
		Map <String, List<Double>> node2Yline = new  HashMap<String, List<Double>>();
		
		for(int i = 0; i < nrow; i++){
			double Y  =  (rowDendY1 + (i+0.5)*boxHeight);
			double X1 =  (rowDendX2 - rowDendBase - rowDendHeight * node2dist.get(node2parent.get(rownames.get(i))));
			double X2 = (double)rowDendX2;
			node2Y.put(rownames.get(i), Y);
			node2X1.put(rownames.get(i), X1);
			node2X2.put(rownames.get(i), X2);
			List <Double> tmp = new ArrayList<Double>();
			tmp.add(X1);
			tmp.add(Y);
			tmp.add(X2);
			tmp.add(Y);
			node2Xline.put(rownames.get(i), tmp);
		}
		for(int i  = 1; i < nrow-1; i++){
			String node = "node" + i;
			double Y = ( node2Y.get(node2daughter.get(node)[0]) + node2Y.get(node2daughter.get(node)[1]) )/2;
			double X1 = (rowDendX2 - rowDendBase - rowDendHeight * node2dist.get(node2parent.get(node)));
			double X2 = node2X1.get(node2daughter.get(node)[0]);
			node2Y.put(node, Y);
			node2X1.put(node, X1);
			node2X2.put(node, X2);
			List <Double> tmp = new ArrayList<Double>();
			tmp.add(X1);
			tmp.add(Y);
			tmp.add(X2);
			tmp.add(Y);
			node2Xline.put(node, tmp);
			List <Double> tmp2 = new ArrayList<Double>();
			tmp2.add(X2);
			tmp2.add((double)node2Y.get(node2daughter.get(node)[0]));
			tmp2.add(X2);
			tmp2.add((double)node2Y.get(node2daughter.get(node)[1]));
			node2Yline.put(node, tmp2);
		}
		List <Double> tmp = new ArrayList<Double>();
		tmp.add((double) rowDendX1);
		tmp.add((double)node2Y.get(node2daughter.get("node" + (nrow-1))[0]));
		tmp.add((double) rowDendX1);
		tmp.add((double)node2Y.get(node2daughter.get("node" + (nrow-1))[1]));
		node2Yline.put("node" + (nrow-1), tmp);
	  
	   if(lockRowHighlight == false){
		   highlightedRowDendNode = "";
		   highlightedRowDendNodeDownstream = new HashSet <String>();
		   highlightedRownameSet = new HashSet<String>();
		   highlightedRowname = new ArrayList<String>();
		   if(mouseX >= rowDendX1 & mouseX <= rowDendX2 & (mouseY - (double)node2Y.get(node2daughter.get("node" + (nrow-1))[1]))*( mouseY - (double)node2Y.get(node2daughter.get("node" + (nrow-1))[0])) <=0){
			   highlightedRowDendNode = "node" + (nrow-1);
		   }
		   for(int i= nrow-2; i > 0; i--){
			   String node = "node" + i;
			   if(mouseX >= node2X2.get(node) & mouseX <= rowDendX2 & (mouseY - (double)node2Y.get(node2daughter.get(node)[1]))*( mouseY - (double)node2Y.get(node2daughter.get(node)[0])) <=0){
				   highlightedRowDendNode = node;
			   }		   
		   }
		   if(!highlightedRowDendNode.equals("")){
			   highlightedRowDendNodeDownstream = new HashSet <String>(H.getDownstreamNodes(highlightedRowDendNode));
		   }
	   }	   
	 
	   strokeWeight(rowDendLineWidth);
	   for(int i = 0; i < nrow; i++){
		   List <Double> xy  = node2Xline.get(rownames.get(i));
		   if( highlightedRowDendNodeDownstream == null){
			   stroke(color(0),alpha);
		   }else if(highlightedRowDendNodeDownstream.contains(rownames.get(i))){
			   stroke(color(colorForHighlight),alpha);
			   if(lockRowHighlight == false){
				   highlightedRowname.add(rownames.get(i));
				   highlightedRownameSet.add(rownames.get(i));
			   }
		   }else{
			   stroke(color(0),alpha);
		   }
		   line((float)(double)(xy.get(0)),(float)(double)(xy.get(1)),(float)(double)(xy.get(2)),(float)(double)(xy.get(3)));   			   
	   }
	   for(int i  = 1; i < nrow-1; i++){
		   List <Double> xy  = node2Xline.get("node" + i);
		   List <Double> xy2  = node2Yline.get("node" + i);
		   if(highlightedRowDendNodeDownstream == null){
			   stroke(color(0),alpha);
		   }else if(highlightedRowDendNodeDownstream.contains("node" + i)){
			   stroke(color(colorForHighlight),alpha);
		   }else{
			   stroke(color(0),alpha);
		   }
			   
		   line((float)(double)(xy.get(0)),(float)(double)(xy.get(1)),(float)(double)(xy.get(2)),(float)(double)(xy.get(3)));  
		   
		   if(highlightedRowDendNodeDownstream == null){
			   stroke(color(0),alpha);
		   }else if(highlightedRowDendNodeDownstream.contains("node" + i) ||highlightedRowDendNode.equals("node" + i) ){
			   stroke(color(colorForHighlight),alpha);
		   }else{
			   stroke(color(0),alpha);
		   }
			   
		   line((float)(double)(xy2.get(0)),(float)(double)(xy2.get(1)),(float)(double)(xy2.get(2)),(float)(double)(xy2.get(3))); 
	   }
	   if(nrow > 2){
	   if(highlightedRowDendNode == null){
		   stroke(color(0),alpha);
	   }else if(highlightedRowDendNode.equals("node" + (nrow-1)) ){
		   stroke(color(colorForHighlight),alpha);
	   }else{
		   stroke(color(0),alpha);
	   }
	   List <Double> xy2  = node2Yline.get("node" + (nrow-1));
	   line((float)(double)(xy2.get(0)),(float)(double)(xy2.get(1)),(float)(double)(xy2.get(2)),(float)(double)(xy2.get(3)));
	   }
	}
	
	private void drawColDendrogram(){
		drawColDendrogram(M.getColClustering());
	}
	
	private void drawColDendrogram(HierarchicalClustering H){
		float colDendX1 =  profileX1;
		float colDendY1 = marginY;
		float colDendY2 =  marginY + colDendHeight+ colDendBase;
		Map <String, Double> node2dist = H.getNode2distMap();
		Map <String, String[]> node2daughter = H.getNode2daughterMap();
		Map <String, String> node2parent = H.getNode2parentMap();
		List <String> colnames = H.getSortedTerminalNodes();
		int ncol = colnames.size();
		Map <String, Double> node2X = new  HashMap<String, Double>();
		Map <String, Double> node2Y1 = new  HashMap<String, Double>();
		Map <String, Double> node2Y2 = new  HashMap<String, Double>();
		
		Map <String, List<Double>> node2Yline = new  HashMap<String, List<Double>>();
		Map <String, List<Double>> node2Xline = new  HashMap<String, List<Double>>();
		
		
		for(int i = 0; i < ncol; i++){
			double X  =  (colDendX1 + (i+0.5)*boxWidth);
			double Y1 =  (colDendY2 - colDendBase - colDendHeight * node2dist.get(node2parent.get(colnames.get(i))));
			double Y2 = (double)colDendY2;
			node2X.put(colnames.get(i), X);
			node2Y1.put(colnames.get(i), Y1);
			node2Y2.put(colnames.get(i), Y2);
			List <Double> tmp = new ArrayList<Double>();
			tmp.add(Y1);
			tmp.add(X);
			tmp.add(Y2);
			tmp.add(X);
			node2Yline.put(colnames.get(i), tmp);
		}
		for(int i  = 1; i < ncol-1; i++){
			String node = "node" + i;
			double X = ( node2X.get(node2daughter.get(node)[0]) + node2X.get(node2daughter.get(node)[1]) )/2;
			double Y1 = (colDendY2 - colDendBase - colDendHeight * node2dist.get(node2parent.get(node)));
			double Y2 = node2Y1.get(node2daughter.get(node)[0]);
			node2X.put(node, X);
			node2Y1.put(node, Y1);
			node2Y2.put(node, Y2);
			List <Double> tmp = new ArrayList<Double>();
			tmp.add(Y1);
			tmp.add(X);
			tmp.add(Y2);
			tmp.add(X);
			node2Yline.put(node, tmp);
			List <Double> tmp2 = new ArrayList<Double>();
			tmp2.add(Y2);
			tmp2.add((double)node2X.get(node2daughter.get(node)[0]));
			tmp2.add(Y2);
			tmp2.add((double)node2X.get(node2daughter.get(node)[1]));
			node2Xline.put(node, tmp2);
		}
		List <Double> tmp = new ArrayList<Double>();
		tmp.add((double) colDendY1);
		tmp.add((double)node2X.get(node2daughter.get("node" + (ncol-1))[0]));
		tmp.add((double) colDendY1);
		tmp.add((double)node2X.get(node2daughter.get("node" + (ncol-1))[1]));
		node2Xline.put("node" + (ncol-1), tmp);

		if(lockColHighlight == false){
			 highlightedColDendNode = "";
			 highlightedColDendNodeDownstream = new HashSet <String>();
			 highlightedColnameSet = new HashSet<String>();
			 highlightedColname = new ArrayList<String>();
			 if(mouseY >= colDendY1 & mouseY <= colDendY2 & (mouseX - (double)node2X.get(node2daughter.get("node" + (ncol-1))[1]))*( mouseX - (double)node2X.get(node2daughter.get("node" + (ncol-1))[0])) <=0){
				 highlightedColDendNode = "node" + (ncol-1);
			 }
			 for(int i= ncol-2; i > 0; i--){
				 String node = "node" + i;
				 if(mouseY >= node2Y2.get(node) & mouseY <= colDendY2 & (mouseX - (double)node2X.get(node2daughter.get(node)[1]))*( mouseX - (double)node2X.get(node2daughter.get(node)[0])) <=0){
					 highlightedColDendNode = node;
				 }			   
			 }
			 if(!highlightedColDendNode.equals("")){
				 highlightedColDendNodeDownstream = new HashSet <String>(H.getDownstreamNodes(highlightedColDendNode));
			 }
		}
		
		strokeWeight(colDendLineWidth);
		for(int i = 0; i < ncol; i++){
			List <Double> xy  = node2Yline.get(colnames.get(i));
			if(highlightedColDendNodeDownstream == null){
				   stroke(color(0),alpha);
			   }else if(highlightedColDendNodeDownstream.contains(colnames.get(i))){
				stroke(color(colorForHighlight),alpha);
				if(lockColHighlight == false){
					highlightedColname.add(colnames.get(i));
					highlightedColnameSet.add(colnames.get(i));
				}
			}else{
				stroke(color(0),alpha);
			}
			line((float)(double)(xy.get(1)),(float)(double)(xy.get(0)),(float)(double)(xy.get(3)),(float)(double)(xy.get(2)));		
		}
		
		for(int i  = 1; i < ncol-1; i++){
			List <Double> xy  = node2Yline.get("node" + i);
			List <Double> xy2  = node2Xline.get("node" + i);
			if(highlightedColDendNodeDownstream == null){
				   stroke(color(0),alpha);
			   }else if(highlightedColDendNodeDownstream.contains("node" + i)){
				stroke(color(colorForHighlight),alpha);
			}else{
				stroke(color(0),alpha);
			}
				   
			line((float)(double)(xy.get(1)),(float)(double)(xy.get(0)),(float)(double)(xy.get(3)),(float)(double)(xy.get(2))); 
			
			if(highlightedColDendNodeDownstream == null){
				   stroke(color(0),alpha);
			 }else if(highlightedColDendNodeDownstream.contains("node" + i) ||highlightedColDendNode.equals("node" + i) ){
				stroke(color(colorForHighlight),alpha);
			}else{
				stroke(color(0),alpha);
			}
				   
			line((float)(double)(xy2.get(1)),(float)(double)(xy2.get(0)),(float)(double)(xy2.get(3)),(float)(double)(xy2.get(2))); 
	   }
		if(ncol > 2){
		if(highlightedColDendNode==null){
			stroke(color(0),alpha);
		}else if(highlightedColDendNode.equals("node" + (ncol-1)) ){
			stroke(color(colorForHighlight),alpha);
		}else{
			stroke(color(0),alpha);
		}
		List <Double> xy2  = node2Xline.get("node" + (ncol-1));
		line((float)(double)(xy2.get(1)),(float)(double)(xy2.get(0)),(float)(double)(xy2.get(3)),(float)(double)(xy2.get(2)));
		}
	}

	protected void drawRowlabels(MyMat M){
		try {
		List <String> rowname = M.getRowNames();
		int nrow = rowname.size();
		textSize(rowLabelFontSize);
		textAlign(LEFT, CENTER);
		fill(0,alpha);
		List <Integer> tmp = new ArrayList<Integer>();
		for(int i = 0; i < nrow; i++){
			if(!highlightedLable.equals("") & rowname.get(i).matches(".*" + highlightedLable + ".*")){
				tmp.add(i);
				continue;
			}else if(highlightedRowname!=null){
				if(highlightedRowname.contains(rowname.get(i))){
					tmp.add(i);
					continue;
				}
			}
			text(rowname.get(i), profileX2+spaceBetweenRowLablelAndProfile, (float)(profileY1 + (i+0.5)*boxHeight));
		}
		
		fill(colorForHighlight,alpha);
		for(Integer i: tmp){	
			text(rowname.get(i), profileX2+spaceBetweenRowLablelAndProfile, (float)(profileY1 + (i+0.5)*boxHeight));
		}
		} catch (Exception e) {
			highlightedLable = "";
			drawRowlabels();
		}
		
	}

	protected void drawCollabels(MyMat M){
		try {
		List <String> colname = M.getColNames();
		int ncol = colname.size();
		textSize(colLabelFontSize);
		textAlign(RIGHT,CENTER);
		fill(0,alpha);
		List <Integer> tmp = new ArrayList<Integer>();
		for(int i = 0; i < ncol; i++){
			if(!highlightedLable.equals("") & colname.get(i).matches(".*" + highlightedLable + ".*")){
				tmp.add(i);
				continue;
			}else if(highlightedColname!=null){
				if(highlightedColname.contains(colname.get(i))){
					tmp.add(i);
					continue;
				}
			}
			float x = (float)(profileX1 + (i+0.5)*boxWidth);
			float y =  profileY2+spaceBetweenColLablelAndProfile;
			rotatedText(colname.get(i),x,y, HALF_PI*3);
		}
		fill(colorForHighlight,alpha);
		for(Integer i: tmp){	
			float x = (float)(profileX1 + (i+0.5)*boxWidth);
			float y =  profileY2+spaceBetweenColLablelAndProfile;
			rotatedText(colname.get(i),x,y, HALF_PI*3);
		}
		} catch (Exception e) {
			highlightedLable = "";
			drawCollabels();
		}
	}
	
	protected  void drawRowlabels(){
		drawRowlabels(M);
	}

	protected  void drawProfile(){
		drawProfile(M);
	}
	protected  void drawCollabels(){
		drawCollabels(M);
	}
	
	protected  void drawProfile(MyMat M){
		int ncol = M.colSize();
		int nrow = M.rowSize();
		float max;
		float min;
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
		max = (float)MyFunc.max(M.asList());
		min = (float)MyFunc.min(M.asList());
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
		
	
   protected  void drawProfile(MyMat M, float profileX1, float profileX2, float profileY1, float profileY2){
			int ncol = M.colSize();
			int nrow = M.rowSize();
			float boxWidth = (profileX2-profileX1)/ncol;
			float boxHeight = (profileY2-profileY1)/nrow;
			float max;
			float min;
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
			max = (float)MyFunc.percentile(M.asList(),0.99);
			min = (float)MyFunc.percentile(M.asList(),0.01);
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
	

   
   protected void drawRowAnnotBand(){
		float rowAnnotX1 = marginX + (M.getRowClustering()!=null?rowDendHeight+rowDendBase+spaceNext2RowDend:0);
		float rowAnnotY1 = profileY1;
		List <String> annotTypes = M.getRowAnnotationMatrix().getColNames();
		for(int i = 0; i < annotTypes.size(); i++){
			String annotType = annotTypes.get(i);
			List <String> annotString = M.getRowAnnotationMatrix().getCol(annotType);
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
		float colAnnotY1 = marginY + (M.getColClustering()!=null?colDendHeight+colDendBase+spaceNext2ColDend:0);
		float colAnnotX1 = profileX1;
		List <String> annotTypes = M.getColAnnotationMatrix().getColNames();
		for(int i = 0; i < annotTypes.size(); i++){
			String annotType = annotTypes.get(i);
			List <String> annotString = M.getColAnnotationMatrix().getCol(annotType);
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
	
	protected  void drawProfileIndicator(){
		drawProfileIndicator(M);
	}
	
	
	protected void drawProfileIndicator(MyMat myMat){
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
	}
	
	protected void drawRowAnnotIndicator(){
		float rowAnnotX1 = marginX + (M.getRowClustering()!=null?rowDendHeight+rowDendBase+spaceNext2RowDend:0);
		float rowAnnotY1 = profileY1;
		double x = ((mouseX-rowAnnotX1)/rowAnnotBandwidth);
		double y = ((mouseY-rowAnnotY1)/boxHeight);
		int j = (int)x;
		int i = (int)y;
		StringMat A = M.getRowAnnotationMatrix();
		if(x > 0 && y > 0 && i>=0 && i < A.rowSize() && j>=0 && j < A.colSize()){
			annotIndicatorDrawn = true;
			String label = A.getRowNames().get(i) + "\n" + A.getColNames().get(j) + "\n" + A.get(i, j);
			text(label, mouseX, mouseY);		
		}	
	}
	
	protected void drawColAnnotIndicator(){
		float colAnnotY1 = marginY + (M.getColClustering()!=null?colDendHeight+colDendBase+spaceNext2ColDend:0);
		float colAnnotX1 = profileX1;
		double x = ((mouseX-colAnnotX1)/boxWidth);
		double y = ((mouseY-colAnnotY1)/colAnnotBandwidth);
		int i = (int)x;
		int j = (int)y;
		StringMat A = M.getColAnnotationMatrix();
		if(x > 0 && y > 0 && i>=0 && i < A.rowSize() && j>=0 && j < A.colSize()){
			annotIndicatorDrawn = true;
			String label = A.getRowNames().get(i) + "\n" + A.getColNames().get(j) + "\n" + A.get(i, j);
			text(label, mouseX, mouseY);		
		}
	}
	
	protected void calculateSize(){
		/*if(boxWidth != 0){
			profileWidth = boxWidth*ncol;	
		}
		if(boxHeight != 0){
			profileHeight = boxHeight*nrow;
		}*/
		
		profileX1 = marginX + (M.getRowClustering()!=null?rowDendHeight+rowDendBase+spaceNext2RowDend:0) 
			+ (M.getRowAnnotationMatrix() != null ?rowAnnotBandwidth*M.getRowAnnotationMatrix().colSize() + spaceNext2RowAnnot:0); 
		profileX2 = profileX1 + profileWidth;
		Width  =  profileX2 + spaceBetweenRowLablelAndProfile + rowLabelSpaceWidth + marginX;
		boxWidth = (profileX2-profileX1)/M.ncol;
		
		profileY1 = marginY + (M.getColClustering()!=null?colDendHeight+colDendBase+spaceNext2ColDend:0)
			+ (M.getColAnnotationMatrix() != null ?colAnnotBandwidth*M.getColAnnotationMatrix().colSize() + spaceNext2ColAnnot:0); 
		profileY2 = profileY1 + profileHeight;
		Height =  profileY2 + spaceBetweenColLablelAndProfile + colLabelSpaceWidth + marginY;
		boxHeight = (profileY2-profileY1)/M.nrow;
	}
	
	protected void rotatedText(String s, float x, float y, float angle){
		translate(x, y);
		rotate(angle);
		text(s, 0, 0);
		rotate(-angle);
		translate(-x, -y);
	}
	
	protected void setBoxWidth(double d){
		boxWidth = (float)d;
	}
	protected void setBoxHeight(double d){
		boxHeight = (float)d;
	}
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
	public void setProfileColor(String  lowColor, String midColor, String highColor){
		color1 = MyColor.getRGBhex(lowColor);
		color2 = MyColor.getRGBhex(midColor);
		color3 = MyColor.getRGBhex(highColor);	
	}
	
	
	public void setProfileColor(String  lowColor, String highColor){
		color1 = MyColor.getRGBhex(lowColor);
		color2 = PApplet.lerpColor(MyColor.getRGBhex(lowColor), MyColor.getRGBhex(highColor), (float) 0.5, RGB);
		color3 = MyColor.getRGBhex(highColor);	
	}
	
	public void scaleColorByRow(){
		colorScalingByRow = true;
		colorScalingByCol = false;
	}
	
	public void scaleColorByColumn(){
		colorScalingByRow = false;
		colorScalingByCol = true;
	}
	

}
