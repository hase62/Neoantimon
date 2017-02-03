package mutation;

import java.applet.Applet;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.cli.*;

import sun.reflect.Reflection;
import utility.*;

public class MutHeatmap extends Heatmap{
	
	private static final long serialVersionUID = 1L;

	public MutHeatmap(Applet C) throws Exception {
		super(C);
	}

	public static void sortMatrix(MyMat M){
		List<Double> colmean = M.getColMeans();
		Map <String,Double>  tmp = new HashMap<String, Double>();
		for(int i = 0; i<M.rowSize(); i++){
			tmp.put(M.getRowNames().get(i), colmean.get(i));
		}
		List <String> sortedRowname = MyFunc.sortKeysByDescendingOrderOfValues(tmp);
		
		List <String> sortedColname = new ArrayList<String>(); 
		Set <String> seen = new HashSet<String>();
		for(int i=0; i<M.rowSize(); i++){
			String r = sortedRowname.get(i);
			List <Double> s = M.getRow(r);
			tmp.clear();
			for(int j = 0; j<M.colSize(); j++){
				if(s.get(j)>0 && !seen.contains(M.getColNames().get(j))){
					tmp.put(M.getColNames().get(j), s.get(j));
					seen.add(M.getColNames().get(j));
				}
			}
			sortedColname.addAll(MyFunc.sortKeysByDescendingOrderOfValues(tmp));
		}
		sortedColname.addAll(MyFunc.diff(M.getColNames(),sortedColname));
		M.reorderRows(sortedRowname);
		M.reorderCols(sortedColname);	
	}	
	//remove empty sub matrix
	public static void sortMatrix2(MyMat M){
		List<Double> colmean = M.getColMeans();
		Map <String,Double>  tmp = new HashMap<String, Double>();
		for(int i = 0; i<M.rowSize(); i++){
			tmp.put(M.getRowNames().get(i), colmean.get(i));
		}
		List <String> sortedRowname = MyFunc.sortKeysByDescendingOrderOfValues(tmp);
		
		List <String> sortedColname = new ArrayList<String>(); 
		Set <String> seen = new HashSet<String>();
		for(int i=0; i<M.rowSize(); i++){
			String r = sortedRowname.get(i);
			List <Double> s = M.getRow(r);
			tmp.clear();
			for(int j = 0; j<M.colSize(); j++){
				if(s.get(j)>0 && !seen.contains(M.getColNames().get(j))){
					tmp.put(M.getColNames().get(j), s.get(j));
					seen.add(M.getColNames().get(j));
				}
			}
			sortedColname.addAll(MyFunc.sortKeysByDescendingOrderOfValues(tmp));
		}
		M.reorderRows(sortedRowname);
		M.reorderCols(sortedColname);	
	}	
	
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("r", "rmenp", false, "remove empty submatrix");
		options.addOption("W", "writepdf", true, "write pdf file");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] tabfile ", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 1){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] tabfile", options);
			return;
		}
			
		ClusteredMyMatWithAnnotation M = new ClusteredMyMatWithAnnotation(argList.get(0));
		if(commandLine.hasOption("r")){
			sortMatrix2(M);
		}else{
			sortMatrix(M);
		}
		if(commandLine.hasOption("W")){
			ClusteredMyMatViewer MV = new ClusteredMyMatViewer(M);
			MV.setProfileColor("lightgray", "red");
			MV.setOutFile(commandLine.getOptionValue("W"));
			MV.useAutoStop();
			Heatmap H = new Heatmap(MV);
			H.setVisible(true);
		}else{	
			InteractiveClusteredMyMatViewer  MV = new InteractiveClusteredMyMatViewer(M);
			MV.setProfileColor("lightgray", "red");
			Heatmap H = new Heatmap(MV);
			H.setVisible(true);
		}
	}
}
