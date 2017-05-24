package tensor;

import java.util.List;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

import processing.core.PApplet;
import utility.MyFunc;
import utility.StringMat;
import utility.HierarchicalClustering.ClusterDistFuncType;

public class test2 {
	public static void main(String[] args) throws Exception {
	Options options = new Options();
	options.addOption("k", "k", true, "the number of clusters");
	options.addOption("s", "sliceord", true,  "sliced order");
	HelpFormatter formatter = new HelpFormatter();
	CommandLineParser parser = new BasicParser();
	CommandLine commandLine;
	try{
		commandLine = parser.parse(options, args);
	}catch (Exception e) {
		formatter.printHelp("Expression [options] tensorBinaryFile", options);
		return;
	}
	List <String> argList = commandLine.getArgList();
	if(argList.size() != 1){
		formatter.printHelp("Expression [options] tensorBinaryFile ", options);
		return;
	}
	ClusteredOrder3TensorWithAnnotation T;
	int k = 3;
	T = new ClusteredOrder3TensorWithAnnotation(ClusteredOrder3Tensor.readFromBinary(argList.get(0)));
	

	
	if(commandLine.hasOption("k")){
		k = Integer.valueOf(commandLine.getOptionValue("k"));
	}
	
	int s = 3;
	if(commandLine.hasOption("s")){
		s = Integer.valueOf(commandLine.getOptionValue("s"));
	}
	
	if(s==1){
		System.out.print(new StringMat("order1cluster", T.getOrder1Clustering().getCutTreeMap(k)));	
	}else if(s==2){
		System.out.print(new StringMat("order1cluster", T.getOrder2Clustering().getCutTreeMap(k)));	
	}else if(s==3){
		System.out.print(new StringMat("order1cluster", T.getOrder3Clustering().getCutTreeMap(k)));	
	}
}















}