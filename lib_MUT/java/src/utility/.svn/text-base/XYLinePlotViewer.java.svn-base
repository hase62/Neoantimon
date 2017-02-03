package utility;

import java.util.*;

public class XYLinePlotViewer extends XYPlotViewer{
	private static final long serialVersionUID = 8657200209454566252L;
	public XYLinePlotViewer() {
		super();
	}
	public XYLinePlotViewer(Map <String, List<Double[]>> XY) {
		super(XY);
	}

	protected void drawPlot(){
		strokeWeight(pointSize); 
		for(String Yid: Yids){
			stroke(Yid2color.get(Yid));
			Map <String, Double> index2x = new LinkedHashMap<String, Double>();
			int j = 0;
			List <Double[]> xy = new ArrayList<Double[]>(XY.get(Yid));
			for(Double[] tmp : xy){
				index2x.put(((Integer)j).toString(), tmp[0]);
				j++;
			}
			List <String> index = MyFunc.sortKeysByAscendingOrderOfValues(index2x);
			
			float x1 = (float)(double) xy.get(Integer.valueOf(index.get(0)))[0];
			float y1 = (float)(double) xy.get(Integer.valueOf(index.get(0)))[1];
			float x2, y2;			
			for(int i = 1, n = index.size();i < n; i++){
				x2 = (float)(double)xy.get(Integer.valueOf(index.get(i)))[0];
				y2 = (float)(double) xy.get(Integer.valueOf(index.get(i)))[1];
				line(XValue2Coordinate(x1), YValue2Coordinate(y1), XValue2Coordinate(x2), YValue2Coordinate(y2));
				x1 = x2;
				y1 = y2;
			}
		}
	}
	
	
}
