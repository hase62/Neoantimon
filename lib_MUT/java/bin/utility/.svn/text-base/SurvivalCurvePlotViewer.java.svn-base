package utility;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class SurvivalCurvePlotViewer extends XYLinePlotViewer {
	private static final long serialVersionUID = -1092057020806991908L;
	private Double Pvalue = null;
	public void setPvalue(double p){
		Pvalue = p;
	}
	public  SurvivalCurvePlotViewer() {
		super();
		setYRange(-0.02, 1.02);
		calculateYAxisScale();
		setXVariableName("time");
	}
	public  SurvivalCurvePlotViewer(Map <String, List<Double[]>> XY) {
		super(XY);
		setYRange(-0.02, 1.02);
		calculateYAxisScale();
		setXVariableName("time");
	}
	protected void drawPlot(){
		strokeWeight(pointSize); 
		for(String Yid: Yids){
			stroke(Yid2color.get(Yid));
			Map <String, Double> index2x = new LinkedHashMap<String, Double>();
			int j = 0;
			List <Double[]> xy = new ArrayList<Double[]>(XY.get(Yid));
			Double[] tmp = {0.0, 1.0};
			xy.add(tmp);
			for(Double[] d : xy){
				index2x.put(((Integer)j).toString(), d[0]);
				j++;
			}
			List <String> index = MyFunc.sortKeysByAscendingOrderOfValues(index2x);
			
			float x1 = (float)(double) xy.get(Integer.valueOf(index.get(0)))[0];
			float y1 = (float)(double) xy.get(Integer.valueOf(index.get(0)))[1];
			float x2, y2;			
			for(int i = 1, n = index.size();i < n; i++){
				x2 = (float)(double)xy.get(Integer.valueOf(index.get(i)))[0];
				y2 = (float)(double) xy.get(Integer.valueOf(index.get(i)))[1];
				line(XValue2Coordinate(x1), YValue2Coordinate(y1), XValue2Coordinate(x2), YValue2Coordinate(y1));
				line(XValue2Coordinate(x2), YValue2Coordinate(y1), XValue2Coordinate(x2), YValue2Coordinate(y2));
				x1 = x2;
				y1 = y2;
			}
		}
	}
	protected void drawVariableLabeles(){
		fill(MyColor.getRGBhex("BLACK"));
		textSize(XVariableLabelFontSize);
		textAlign(CENTER);
		text(Xid, (float)(plotX1 + plotWidth *0.5), (float)(plotY2 +  XAxisLabelSpaceWidth + XVariableLabelSpaceWidth*0.5));
		textSize(YVariableLabelFontSize);
		textAlign(CENTER);
		if(Pvalue!=null){
			String s = "P = " + MyFunc.toString(Pvalue,3);
			Yids.add(s);
			Yid2color.put(s, MyColor.getRGBhex("BLACK"));
		}
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
		if(Pvalue!=null){
		 Yids.remove(Yids.size()-1);
		}
	}	
}
