package utility;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

import org.apache.commons.cli.*;

import sun.reflect.Reflection;
import utility.*;

public class ROCanalysis {
	double min;
	double max;
	List <Double> positives;
	List <Double> negatives;
	List <Double> all;
	Map <Double, Double> ROCpoints;
	List <Double[]> ROCpoints2;
	double AUC = 0;
	
	boolean usePR = false;
	
	
	public ROCanalysis(List <Double> positives, List <Double> negatives){
 		this.positives = positives;
 		this.negatives = negatives;
 		all = new ArrayList<Double>(positives);
 		all.addAll(negatives);
 		min = MyFunc.min(all);
 		max = MyFunc.max(all);
 		ROCpoints = new TreeMap<Double, Double>();
 		ROCpoints2 = new  ArrayList<Double[]>();
 	}
	
	public ROCanalysis(String inFile) throws IOException{
		BufferedReader inputStream = new BufferedReader(new FileReader(inFile));
		String line;
		positives = new ArrayList<Double>();
		positives = new ArrayList<Double>();
		while((line = inputStream.readLine())!= null){
			List<Double> tmp = MyFunc.toDouble(MyFunc.split("\t", line));
			if(tmp.get(1)==0){
				negatives.add(tmp.get(0));
			}else{
				positives.add(tmp.get(0));
			}	
		}
 		all = new ArrayList<Double>(positives);
 		all.addAll(negatives);
 		min = MyFunc.min(all);
 		max = MyFunc.max(all);
 		ROCpoints = new TreeMap<Double, Double>();
 		ROCpoints2 = new  ArrayList<Double[]>();
 	}
	
	private double countAboveCutoff(List <Double> v, double c){
		double i = 0;
		for(Double d: v){
			if(d >= c){
				i++;
			}
		}
		return i;	
	}
	 
	
	public void calculateROC(){
		Collections.sort(all);
		for(int i = 0, n = all.size(); i < n; i++){
			double x, y;
			if(!usePR){
				y =   countAboveCutoff(positives, all.get(i))/positives.size(); // sensitivity
				x =   countAboveCutoff(negatives,all.get(i))/negatives.size(); // 1 - specificity
			}else{
				y =   countAboveCutoff(positives, all.get(i))/countAboveCutoff(all, all.get(i)); // precision 
				x =   countAboveCutoff(positives, all.get(i))/positives.size(); //  recall
			}
			ROCpoints.put(x,y);
			Double[] tmp = {x, y};
			ROCpoints2.add(tmp);
		}
		
	}
	
	public void calculateAUC(){
		List <Double> x = new ArrayList<Double>(ROCpoints.keySet());
		for(int i = 1, n = x.size(); i < n; i++){
			AUC += (ROCpoints.get(x.get(i)) + ROCpoints.get(x.get(i-1)))*(x.get(i) - x.get(i-1))*0.5;
		}
	}
	
 	
	public double getAUC(){
		return AUC;
	}
 	
	public  XYLinePlotViewer  getROCplotViewer(){
		XYLinePlotViewer viewer = new XYLinePlotViewer();
		viewer.setXVariableName("1-Specificity");
		viewer.addXY("Sensitivity", ROCpoints2);
		viewer.setXRange(-0.02, 1.02);
		viewer.setYRange(-0.02, 1.02);
		viewer.calculateXAxisScale();
		viewer.calculateYAxisScale();
		return viewer;
	}
	
	public static void main(String[] args) throws Exception {
		Options options = new Options();
		options.addOption("p", "pr", false, "use pr");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine = null;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] SegFile  probeTsvFile", options);
			System.exit(1);
		}
		List <String> argList = commandLine.getArgList();
		if(!(argList.size() == 1)){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] SegFile probeTsvFile", options);
			System.exit(1);
		}
		ROCanalysis R = new ROCanalysis(argList.get(0));
		if(commandLine.hasOption("p")){
			R.usePR = true;
		}
		R.calculateROC();
		R.calculateAUC();
		System.out.println(R.getAUC());
	}
	
	
}
