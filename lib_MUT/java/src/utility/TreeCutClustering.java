package utility;

import java.util.*;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

import sun.reflect.Reflection;
import tensor.*;
import utility.HierarchicalClustering.ClusterDistFuncType;

public class TreeCutClustering extends HierarchicalClustering{
	
	
	private List <List<String>>  clusterMembers;
	private Map <String, Integer> menber2cluster;
	
	public TreeCutClustering(Dist dist){
		super(dist);
	}
	
	public TreeCutClustering (HierarchicalClustering H){
		super(H);
	}
	
	public void  cut(int k){
		clusterMembers = cutTree(k);
		Collections.sort(clusterMembers, new MyComparator()); 
		menber2cluster = new LinkedHashMap<String, Integer>();
		for(int i =0; i < clusterMembers.size(); i++){
			for(String s: clusterMembers.get(i)){
				menber2cluster .put(s, i+1);
			}
		}
	}
	
	public static TreeCutClustering getTreeCutClusteringFromOrder3Tensor(Order3Tensor T, int sliceOrder){
		return new TreeCutClustering(Order3Tensor.getDistBetweenSlices(T, sliceOrder));
	}

	
		
	public void printGmt(){
		int i = 1;
		for(List<String> L: clusterMembers){
			System.out.println("cluster" + i + "\t" + L.size() + "\t" +MyFunc.join("\t", L));			
			i++;
		}
	}
	
	public void printAnnotMatrix(){
		System.out.println("\tcluster");
		for(String s: menber2cluster.keySet()){
			System.out.println(s + "\t" + menber2cluster.get(s));
		}
	}
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("a", "annot", false, "output annot matrix");
		options.addOption("k", "clstno", true, "the number of clusters");
		options.addOption("c", "cor", false, "use correlation");
		options.addOption("t", "tensor", true, "cluser tensor order");
		options.addOption("C", "complete", false, "perform complete clustering");
		options.addOption("w", "ward", false, "perform ward clustering");
		options.addOption("T", "colunm", false, "apply clustering to columns");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] expfile", options);
			return;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 1){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] expfile ", options);
			return;
		}
		
		MyMat M;
		int k = 5;
		if(commandLine.hasOption("k")){
			k = Integer.valueOf(commandLine.getOptionValue("k"));
		}
		TreeCutClustering T;
		if(commandLine.hasOption("t")){
			T = getTreeCutClusteringFromOrder3Tensor(new Order3Tensor(argList.get(0)), Integer.valueOf(commandLine.getOptionValue("t")));
		}else{
			M = new Expression(argList.get(0));
			if(commandLine.hasOption("w")){
				if(commandLine.hasOption("T")){
					M.transpose();
				}
				T = new TreeCutClustering(getHierarchicalClusteringBasedWardMethod(M));
			}else{
				if(commandLine.hasOption("c")){
					T = new TreeCutClustering(commandLine.hasOption("T")?MyMat.getDistBetweenCols(M, 'c'):MyMat.getDistBetweenRows(M, 'c'));
					T.setDistFunc(DistFuncType.SIMILARITY);
				}else{
					T = new TreeCutClustering(commandLine.hasOption("T")?MyMat.getDistBetweenCols(M, 'e'):MyMat.getDistBetweenRows(M, 'e'));
				}
				if(commandLine.hasOption("C")){
					T.setClusterDistFunc(ClusterDistFuncType.COMPLETE);
				}
			}
		}
		T.perform();
		T.cut(k);
		if(commandLine.hasOption("a")){
			T.printAnnotMatrix();
		}else{
			T.printGmt();
		}
	}
	
	
	
	
	
}
