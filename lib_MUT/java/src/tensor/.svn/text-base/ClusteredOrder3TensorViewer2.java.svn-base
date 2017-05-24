package tensor;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import processing.core.PFont;

import utility.HierarchicalClustering;
import utility.Integrator;
import utility.MyColor;
import utility.MyFunc;
import utility.MyMat;
import utility.StringMat;
import controlP5.Button;
import controlP5.ControlEvent;
import controlP5.ControlP5;
import controlP5.RadioButton;
import controlP5.Slider;
import controlP5.Textfield;
import controlP5.Toggle;


public class ClusteredOrder3TensorViewer2  extends  ClusteredOrder3TensorViewer{

	
	private static final long serialVersionUID = -5287573654558352861L;
	
	protected PFont font = createFont("Bitstream Vera Sans",20);
	
	protected int colorForHighlight = MyColor.getRGBhex("RED");
	
	protected String pdfFile = null;
	protected String labelFile = null;
	

	protected ControlP5 controlP5;
	protected Slider sliceIndexSlider;
	
	
	protected ClusteredOrder3TensorWithAnnotation originalT;  
	protected ClusteredOrder3TensorWithAnnotation preT;
	
	protected int marginY = 50;
	
	public ClusteredOrder3TensorViewer2  (){
	}
	
	public ClusteredOrder3TensorViewer2  (ClusteredOrder3TensorWithAnnotation T){
		super(T);
		originalT = new ClusteredOrder3TensorWithAnnotation(T);
	}
	
	
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
	
	
	protected void setSlider(){
		controlP5.remove("s");
		controlP5.addSlider("s",0,nslice-1,0,(int)marginX,10,100,20);
		sliceIndexSlider =(Slider)controlP5.controller("s");
		sliceIndexSlider.setId(0);
		sliceIndexSlider.setValue(sliceIndex);
		sliceIndexSlider.setLabelVisible(true);
		sliceIndexSlider.setNumberOfTickMarks(nslice);
		//sliceIndexSlider.setColorCaptionLabel(color(0));
		sliceIndexSlider.showTickMarks(true);
	}
	
	protected void setRadioButton(){
		RadioButton   r = controlP5.addRadioButton("radioButton",(int)marginX+130,10);
		  r.setColorLabel(color(0));
		  r.setSize(20, 20);
		  r.setItemsPerRow(3);
		  r.setLabel("sliceOrder");
		  r.setColorLabel(color(0));
		  r.setSpacingColumn(20);
		  r.setNoneSelectedAllowed(false);
		  r.setId(1);
		  addToRadioButton(r,"1",1);
		  addToRadioButton(r,"2",2);
		  addToRadioButton(r,"3",3);
		  r.activate(Integer.toString(slicedOrder));
	}
	protected void setButtons(){
		Button b = controlP5.addButton("Zoom",0,(int)marginX+260,10,30,20);
		b.setId(2);
		Button b2 = controlP5.addButton("Label",0,(int)marginX+300,10,30,20);
		b2.setId(3);
		Button b3 = controlP5.addButton("Pdf",0,(int)marginX+340,10,30,20);
		b3.setId(4);
		Button b4 = controlP5.addButton("Find",0,(int)marginX+380,10,30,20);
		b4.setId(5);
	}
	
	Textfield t;
	protected void setTextfield(){
		t = controlP5.addTextfield("search",(int)marginX+440,10,100,20);
		t.setFocus(true);
	}
	
	protected void addToRadioButton(RadioButton theRadioButton, String theName, int theValue ) {
		  Toggle t = theRadioButton.addItem(theName,theValue);
		  t.captionLabel().setColorBackground(color(255));
		  //t.captionLabel().style().
		  //t.captionLabel().style().backgroundWidth = 10;
		  //t.captionLabel().style().backgroundHeight = 17;
		  t.captionLabel().style().movePadding(5,0,0,2);
		  t.captionLabel().style().moveMargin(1,0,0,0);
		  
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
	 if(theEvent.id()==0){
		sliceIndex = (int)(theEvent.controller().value());
	 }else if(theEvent.id()==1){
		 rotateTensor((int)theEvent.group().value());
	 }else if(theEvent.id()==2){
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
	
	
	protected void rotateTensor(int sliceOrder){
		setSlicedOder(sliceOrder);
		calculateSize();
		setSlider();
		resetSliceIndex();
		lockColHighlight = false;
		lockRowHighlight = false;
	}


	public void setupControlP5(){
		controlP5 = new ControlP5(this);
		controlP5.setColorBackground(colorBgForControlP5);
		controlP5.setColorForeground(colorFgForControlP5);
		controlP5.setColorActive(colorActiveForControlP5);
		setSlider();
		setRadioButton();
		setButtons();
		setTextfield();
	}
	
	
	public void setup(){
		super.setup();
		setupControlP5();
		background(bgColor);
	}
	
	
	public void draw(){
		textFont(font);
		
		if(refresh == false & annotIndicatorDrawn == false  & nrow * ncol > 2500 & eventHappened == false   & frameCount > 2){
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
					drawProfile(getSlice(),profileX1forZoom.value(), profileX2forZoom.value(), profileY1forZoom.value(), profileY2forZoom.value());
				}else{
					drawProfile(getSlice(preT),profileX1forZoom.value(), profileX2forZoom.value(), profileY1forZoom.value(), profileY2forZoom.value());
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
				if(rowClustering!=null){
					drawRowDendrogram();
				}
				if(colClustering!=null){
					drawColDendrogram();
				}
				drawRowlabels();
				drawCollabels();
				if(rowAnnotation!=null){
					drawRowAnnotBand();
					drawRowAnnotIndicator();
				}
				if(colAnnotation!=null){
					drawColAnnotBand();
					drawColAnnotIndicator();
				}
				drawSliceLable();
				if(nrow * ncol <= 2500){
					drawProfileIndicator();
				}
				if(nrow * ncol <= 2500){
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
		if(rowClustering!=null){
			drawRowDendrogram();
		}
		if(colClustering!=null){
			drawColDendrogram();
		}
		drawRowlabels();
		drawCollabels();
		if(rowAnnotation!=null){
			drawRowAnnotBand();
			drawRowAnnotIndicator();
		}
		if(colAnnotation!=null){
			drawColAnnotBand();
			drawColAnnotIndicator();
		}
		drawSliceLable();
	}
	
	
	protected void incrementSliceIndex(){
		int tmp = sliceIndex;
		sliceIndexSlider.setValue(tmp+1);
		sliceIndex = tmp + 1;
		if(sliceIndex == nslice){
			resetSliceIndex();
		}
	}
	
	
	protected void resetSliceIndex(){
		sliceIndexSlider.setValue(0);
		sliceIndex=0;
	}
	
	public void mousePressed(){
		refresh = true;
		if(mouseX >= profileX1 & mouseX <= profileX2 & mouseY>= profileY1 & mouseY <= profileY2){
			incrementSliceIndex();
		}
		if((mouseX -  marginX)*(mouseX - (marginX + rowDendHeight+rowDendBase))<= 0 & (mouseY -profileY1)*(mouseY-profileY2)<=0){
			lockRowHighlight = (lockRowHighlight == true)?false:true;
		}
		if((mouseY -  marginY)*(mouseY - (marginY + colDendHeight+colDendBase))<= 0 & (mouseX -profileX1)*(mouseX-profileX2)<=0){
			lockColHighlight = (lockColHighlight == true)?false:true;
		}	
	}
	
	
	
	protected MyMat getSlice(ClusteredOrder3Tensor T){
		MyMat M;
		if(slicedOrder == 1){
		M =  T.getOrder1Slice(sliceIndex);
		}else if(slicedOrder == 2){
			M = T.getOrder2Slice(sliceIndex);
		}else if(slicedOrder == 3){
			M= T.getOrder3Slice(sliceIndex);
		}else{
			return null;
		}
		return M;
	}
	
	protected MyMat getSlice(){
		return getSlice(T);
	}
	
	
	
	private void zoom(){
		if((highlightedRownameSet != null && !highlightedRownameSet.isEmpty())|| (highlightedColnameSet != null  && !highlightedColnameSet.isEmpty())){
			zoomIn();
		}else{
			zoomOut();
		}
	}
	
	private void zoomIn(){	

		preT = new ClusteredOrder3TensorWithAnnotation(T);
		
		profileX1forZoom = new  Integrator(profileX1);
		profileX2forZoom = new Integrator(profileX2);
		profileY1forZoom = new Integrator(profileY1);
		profileY2forZoom = new Integrator(profileY2);
		if(highlightedRownameSet != null){
			if(!highlightedRownameSet.isEmpty()){
				profileY1forZoom.set(profileY1+getSlice().rowIndexOf(highlightedRowname.get(0))*boxHeight);
				profileY2forZoom.set(profileY1+(getSlice().rowIndexOf(highlightedRowname.get(highlightedRowname.size()-1))+1)*boxHeight);
			}
		}
		if(highlightedColnameSet != null){
			if(!highlightedColnameSet.isEmpty()){
				profileX1forZoom.set(profileX1+getSlice().colIndexOf(highlightedColname.get(0))*boxWidth);
				profileX2forZoom.set(profileX1+(getSlice().colIndexOf(highlightedColname.get(highlightedColname.size()-1))+1)*boxWidth);
			}
		}

		profileX1forZoom.target(profileX1);
		profileX2forZoom.target(profileX2);
		profileY1forZoom.target(profileY1);
		profileY2forZoom.target(profileY2);
		
		zooming = 1;
		
		
		if(slicedOrder == 1){
			if(highlightedRownameSet != null){
				if(!highlightedRownameSet.isEmpty()){
					HierarchicalClustering H = rowClustering.getSubHierarchicalClustering(highlightedRowDendNode);
					T.setOrder2Clustering(H);				
				}
			}
			if(highlightedColnameSet != null){
				if(!highlightedColnameSet.isEmpty()){
					if(!highlightedColnameSet.isEmpty()){
						HierarchicalClustering H = colClustering.getSubHierarchicalClustering(highlightedColDendNode);
						T.setOrder3Clustering(H);				
					}
				}
			}	
			setSlicedOder(1);
		}else if(slicedOrder == 2){
			if(highlightedRownameSet != null){
				if(!highlightedRownameSet.isEmpty()){
					HierarchicalClustering H = rowClustering.getSubHierarchicalClustering(highlightedRowDendNode);
					T.setOrder3Clustering(H);				
				}
			}
			if(highlightedColnameSet != null){
				if(!highlightedColnameSet.isEmpty()){
					if(!highlightedColnameSet.isEmpty()){
						HierarchicalClustering H = colClustering.getSubHierarchicalClustering(highlightedColDendNode);
						T.setOrder1Clustering(H);				
					}
				}
			}
			setSlicedOder(2);
		}else if(slicedOrder == 3){
			if(highlightedRownameSet != null){
				if(!highlightedRownameSet.isEmpty()){
					HierarchicalClustering H = rowClustering.getSubHierarchicalClustering(highlightedRowDendNode);
					T.setOrder1Clustering(H);				
				}
			}
			if(highlightedColnameSet != null){
				if(!highlightedColnameSet.isEmpty()){
					if(!highlightedColnameSet.isEmpty()){
						HierarchicalClustering H = colClustering.getSubHierarchicalClustering(highlightedColDendNode);
						T.setOrder2Clustering(H);				
					}
				}
			}
			setSlicedOder(3);
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
		
		List <String> rowname = getSlice().getRowNames() ;
		List <String> colname = getSlice().getColNames() ;
		preT = new ClusteredOrder3TensorWithAnnotation(T);
		T = new ClusteredOrder3TensorWithAnnotation(originalT);
		if(slicedOrder == 1){
			setSlicedOder(1);
		}else if(slicedOrder == 2){
			setSlicedOder(2);
		}else if(slicedOrder == 3){
			setSlicedOder(3);
		}
		calculateSize();
		
		profileX1forZoom = new  Integrator(profileX1);
		profileX2forZoom = new Integrator(profileX2);
		profileY1forZoom = new Integrator(profileY1);
		profileY2forZoom = new Integrator(profileY2);
		profileX1forZoom.target(profileX1);
		profileX2forZoom.target(profileX2);
		profileY1forZoom.target(profileY1);
		profileY2forZoom.target(profileY2);
		
		if(rowname.size() != getSlice().getRowNames().size()){
			profileY1forZoom.target(profileY1+getSlice().rowIndexOf(rowname.get(0))*boxHeight);
			profileY2forZoom.target(profileY1+(getSlice().rowIndexOf(rowname.get(rowname.size()-1))+1)*boxHeight);
		}
		if(colname.size() != getSlice().getColNames().size()){
			profileX1forZoom.target(profileX1+getSlice().colIndexOf(colname.get(0))*boxWidth);
			profileX2forZoom.target(profileX1+(getSlice().colIndexOf(colname.get(colname.size()-1))+1)*boxWidth);
		}
		zooming = -1;
		setSlider();
	}

	
	private void drawRowDendrogram(){
		drawRowDendrogram(rowClustering);
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
		drawColDendrogram(colClustering);
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
	
	protected void calculateSize(){
		/*if(boxWidth != 0){
			profileWidth = boxWidth*ncol;	
		}
		if(boxHeight != 0){
			profileHeight = boxHeight*nrow;
		}*/
		
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
	
	protected  void drawRowlabels(){
		drawRowlabels(getSlice());
	}

	protected  void drawProfile(){
		drawProfile(getSlice());
	}
	protected  void drawCollabels(){
		drawCollabels(getSlice());
	}
	
   protected  void drawProfile(MyMat M, float profileX1, float profileX2, float profileY1, float profileY2){
			int ncol = M.colSize();
			int nrow = M.rowSize();
			float boxWidth = (profileX2-profileX1)/ncol;
			float boxHeight = (profileY2-profileY1)/nrow;
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
   
   protected void drawSliceLable(){
		textSize(sliceLabelFontSize);
		textAlign(LEFT, TOP);
		//textAlign(RIGHT, BOTTOM);
		fill(0,alpha);
		//text(getSliceLable(), (float)(marginX + rowDendHeight), (float)(marginY + colDendHeight));
		text(getSliceLable(), (float)(profileX2+spaceBetweenRowLablelAndProfile), (float)(profileY2+spaceBetweenColLablelAndProfile));
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
	
	protected  void drawProfileIndicator(){
		drawProfileIndicator(getSlice());
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
		float rowAnnotX1 = marginX + (rowClustering!=null?rowDendHeight+rowDendBase+spaceNext2RowDend:0);
		float rowAnnotY1 = profileY1;
		double x = ((mouseX-rowAnnotX1)/rowAnnotBandwidth);
		double y = ((mouseY-rowAnnotY1)/boxHeight);
		int j = (int)x;
		int i = (int)y;
		StringMat A = rowAnnotation;
		if(x > 0 && y > 0 && i>=0 && i < A.rowSize() && j>=0 && j < A.colSize()){
			String label = A.getRowNames().get(i) + "\n" + A.getColNames().get(j) + "\n" + A.get(i, j);
			text(label, mouseX, mouseY);		
		}	
	}
	
	protected void drawColAnnotIndicator(){
		float colAnnotY1 = marginY + (colClustering!=null?colDendHeight+colDendBase+spaceNext2ColDend:0);
		float colAnnotX1 = profileX1;
		double x = ((mouseX-colAnnotX1)/boxWidth);
		double y = ((mouseY-colAnnotY1)/colAnnotBandwidth);
		int i = (int)x;
		int j = (int)y;
		StringMat A = colAnnotation;
		if(x > 0 && y > 0 && i>=0 && i < A.rowSize() && j>=0 && j < A.colSize()){
			String label = A.getRowNames().get(i) + "\n" + A.getColNames().get(j) + "\n" + A.get(i, j);
			text(label, mouseX, mouseY);		
		}
	}
	
}
	


