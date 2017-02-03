package utility;

import java.util.*;
import processing.core.*;


public class XYPlotViewer extends PApplet{
	private static final long serialVersionUID = 8191079932787492430L;

	protected Map <String, List<Double[]>> XY;
	
	protected float plotWidth = 1000;
	protected float plotHeight = 1000;
	
	protected int marginX = 20;
	protected int marginY = 20;
	
	protected int XVariableLabelFontSize = 20;
	protected int YVariableLabelFontSize = 20;
	protected int XAxisLabelFontsize = 20;
	protected int YAxisLabelFontsize = 20;
	
	protected int XVariableLabelSpaceWidth = 50;
	protected int YVariableLabelSpaceWidth = 100;
	protected int XAxisLabelSpaceWidth = 100;
	protected int YAxisLabelSpaceWidth = 100;
	
	
	
	protected String Xid;
	protected List <String> Yids; 
	
	protected float plotX1;
	protected float plotX2;
	protected float plotY1;
	protected float plotY2;

	
	protected float maxX;
	protected float minX;
	protected float maxY;
	protected float minY;
	
	
	protected float Width;
	protected float Height;
	
	protected int bgColor = MyColor.getRGBhex("WITHE");
	protected int defaultColor = MyColor.getRGBhex("Black");
	protected Map <String, Integer> Yid2color;
	protected int pointSize = 5;
	protected int axisWidth = 2;
	
	List <Double> XAxisScale;
	List <Double> YAxisScale;
	 
	
	
	protected PFont font = createFont("SansSerif",20);
	
	float plotMarginWitdthRate = (float)0.02;
	
	
	public XYPlotViewer(){
		XY = new HashMap<String, List<Double[]>>();
		Xid = "X";
		Yids = new ArrayList<String>();
		Yid2color = new HashMap<String, Integer>();
		caluclateSize();
	}
	
	public void addXY(String YVariableName, List <Double> X, List <Double> Y){
		int m = Math.min(X.size(), Y.size());
		List <Double[]> xyv = new ArrayList<Double[]>();
		for(int i = 0; i < m; i++){
			Double[] xy = {X.get(i), Y.get(i)};
			xyv.add(xy);	
		}
		XY.put(YVariableName, xyv);
		if(!Yids.contains(YVariableName)){
			Yids.add(YVariableName);
			Yid2color.put(YVariableName, MyColor.getRGBhex("BLACK"));
		}
	}
	
	public void addXY(String YVariableName, List <Double[]> XY){
		this.XY.put(YVariableName, XY);
		if(!Yids.contains(YVariableName)){
			Yids.add(YVariableName);
			Yid2color.put(YVariableName, MyColor.getRGBhex("BLACK"));
		}
	}
	
	public XYPlotViewer(Map <String, List<Double[]>> XY){
		this.XY = new HashMap<String, List<Double[]>>(XY);
		Xid = "X";
		Yids = new ArrayList<String>(XY.keySet());
		Collections.sort(Yids);
		Yid2color = new HashMap<String, Integer>();
		for(String s: Yids){
			Yid2color.put(s, MyColor.getRGBhex("BLACK"));
		}
		caluclateSize();
		calculateXYRange();
		calculateXAxisScale();
		calculateYAxisScale();
	}

	public void calculateXYRange(){
		List<Double> X = new ArrayList<Double>();
		for(String s: XY.keySet()){
			for(Double[] d: XY.get(s)){
				X.add(d[0]);
			}
		}
		List<Double> Y = new ArrayList<Double>();
		for(String s: XY.keySet()){
			for(Double[] d: XY.get(s)){
				Y.add(d[0]);
			}
		}
		float maxX = (float) MyFunc.max(X);
		float minX = (float) MyFunc.min(X);
		float maxY = (float) MyFunc.max(Y);
		float minY = (float) MyFunc.min(Y);
		this.maxY = maxY + (maxY - minY)*plotMarginWitdthRate;
		this.minY = minY - (maxY - minY)*plotMarginWitdthRate;
		this.maxX = maxX + (maxX - minX)*plotMarginWitdthRate;
		this.minX = minX - (maxX - minX)*plotMarginWitdthRate;
	}
	public void calculateXAxisScale(){
		double XscaleInterval = Math.pow(10, Math.floor(Math.log10((maxX - minX))));
		XAxisScale= new ArrayList <Double>();
		float minXscale = (float) (Math.ceil(minX/XscaleInterval)* XscaleInterval);
		for(float f = minXscale; f < maxX ; f += XscaleInterval){
			XAxisScale.add((double)f);
		}
		if(XAxisScale.size() < 5){
			XAxisScale.clear();
			float min  = minXscale;	
			for(float f = minXscale; f > minX; f -= 0.5*XscaleInterval){
				min = f;
			}
			for(float f = min; f < maxX ; f += 0.5*XscaleInterval){
				XAxisScale.add((double)f);
			}
		}
		if(XAxisScale.size() < 5){
			XAxisScale.clear();
			float min  = minXscale;	
			for(float f = minXscale; f > minX; f -= 0.2*XscaleInterval){
				min = f;
			}
			for(float f = min; f < maxX ; f += 0.2*XscaleInterval){
				XAxisScale.add((double)f);
			}
		}
	}
	
	public void calculateYAxisScale(){
		double YscaleInterval = Math.pow(10, Math.floor(Math.log10((maxY - minY))));
		
		YAxisScale= new ArrayList <Double>();
		float minYscale = (float) (Math.ceil(minY/YscaleInterval)* YscaleInterval);
		for(float f = minYscale; f < maxY ; f += YscaleInterval){
			YAxisScale.add((double)f);
		}
		if(YAxisScale.size() < 5){
			YAxisScale.clear();
			float min  = minYscale;	
			for(float f = minYscale; f > minY; f -= 0.5*YscaleInterval){
				min = f;
			}
			for(float f = min; f < maxY ; f += 0.5*YscaleInterval){
				YAxisScale.add((double)f);
			}
		}
		if(YAxisScale.size() < 5){
			YAxisScale.clear();
			float min  = minYscale;	
			for(float f = minYscale; f > minY; f -= 0.2*YscaleInterval){
				min = f;
			}
			for(float f = min; f < maxY ; f += 0.2*YscaleInterval){
				YAxisScale.add((double)f);
			}
		}
	}
	
	
	public void setXRange(double min, double max){
		minX = (float) min;
		maxX = (float) max;
	}
	
	public void setYRange(double min, double max){
		minY = (float) min;
		maxY = (float) max;
	}
	
	public void setXVariableName(String s){
			Xid = s;
	}
	
		
	protected  void caluclateSize(){
		plotX1 = marginX  + YVariableLabelSpaceWidth + YAxisLabelSpaceWidth;
		plotX2 = plotX1 + plotWidth;
		Width = plotX2 + marginX;
		plotY1 = marginY;
		plotY2 = plotY1 + plotHeight;
		Height = plotY2 + XAxisLabelSpaceWidth + XVariableLabelSpaceWidth + marginY;		
	}
	
	public void setColorOfY(String Yid, String color){
		if(Yid2color.containsKey(Yid) && MyColor.contains(color)){
			Yid2color.put(Yid, MyColor.getRGBhex(color));
		}
	}
	
	public void setColorOfYautomatically(){
		float d = (float)(1.0/(Yids.size()-1));
		for(int j = 0 ; j < Yids.size(); j++){
			Yid2color.put(Yids.get(j),(Integer)MyColor.lerpColor(MyColor.getRGBhex("BLUE"),MyColor.getRGBhex("YELLOW"),MyColor.getRGBhex("RED"), d*j,HSB));
		}
	}
	
	public void setup(){
		size(Math.round(Width), Math.round(Height));
		textFont(font);
		noLoop();
	}
	
	
	public void draw(){
		background(bgColor);
		drawPlot();
		drawAxes();
		drawAxisLabales();
		 drawVariableLabeles();
	}
	
	
	protected float XValue2Coordinate(float x){
		return plotX1 + plotWidth*(x-minX)/(maxX-minX);
	}
	protected float YValue2Coordinate(float y){
		return plotY2-plotHeight*(y-minY)/(maxY-minY);
	}
	
	protected void drawPlot(){
		strokeWeight(pointSize);
		for(String Yid: Yids){
			stroke(Yid2color.get(Yid));
			for(Double[] xy : XY.get(Yid)){
				float x = (float)(double) xy[0] ;
				float y = (float)(double) xy[1];
				if(x < maxX && x > minX && y < maxY && y > minY){
					point(XValue2Coordinate(x), YValue2Coordinate(y));		
				}
			}
		}
	}
		
	protected void drawAxes(){
		strokeWeight(axisWidth);
		stroke(MyColor.getRGBhex("BLACK"));
		noFill();
		rect(plotX1, plotY1, plotWidth, plotHeight);
		
		float XAxisScaleLength = XAxisLabelSpaceWidth /10;
		float YAxisScaleLength = YAxisLabelSpaceWidth /10;
		
		for(Double d: XAxisScale){
			float x = XValue2Coordinate((float)(double)d);
			line(x, plotY2, x, plotY2+XAxisScaleLength);
		}
		
		for(Double d: YAxisScale){
			float y = YValue2Coordinate((float)(double)d);
			line(plotX1-YAxisScaleLength,y, plotX1, y);
		}
	}
	
	protected void drawAxisLabales(){
		fill(MyColor.getRGBhex("BLACK"));
		textSize(XAxisLabelFontsize);
		textAlign(CENTER);
		for(Double d: XAxisScale){
			text(MyFunc.toString(d,2), XValue2Coordinate((float)(double)d), (float)(plotY2 + 0.5*XAxisLabelSpaceWidth));	
		}
		textSize(YAxisLabelFontsize);
		textAlign(CENTER,CENTER);
		for(Double d: YAxisScale){
			text(MyFunc.toString(d,2), (float)(plotX1 - 0.5*XAxisLabelSpaceWidth), YValue2Coordinate((float)(double)d));	
		}
	}
	
	protected void drawVariableLabeles(){
		fill(MyColor.getRGBhex("BLACK"));
		textSize(XVariableLabelFontSize);
		textAlign(CENTER);
		text(Xid, (float)(plotX1 + plotWidth *0.5), (float)(plotY2 +  XAxisLabelSpaceWidth + XVariableLabelSpaceWidth*0.5));
		textSize(YVariableLabelFontSize);
		textAlign(CENTER);
		if(Yids.size()%2 == 0){
			int k = 20;
			double d = (plotY2- plotY1)/k;
			double y = plotY1  +  0.5 * d  +   (((k - Yids.size())/2) -1)  *d ; 
			for(int i = 0; i < Yids.size(); i++){
				y += d;
				fill(Yid2color.get(Yids.get(i)));
				text(Yids.get(i),(float)(marginX+YVariableLabelSpaceWidth*0.5),(float)y);
			}
		}else{
			int k = 21;
			double d = (plotY2- plotY1)/k;
			double y = plotY1  + 0.5 * d +  ((k +1 - Yids.size())/2 -  1)*d ; 
			for(int i = 0; i < Yids.size(); i++){
				y += d;
				fill(Yid2color.get(Yids.get(i)));
				text(Yids.get(i),(float)(marginX+YVariableLabelSpaceWidth*0.5),(float)y);
			}
			
		}	
		
	}
	






}
