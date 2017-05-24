package tensor;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.util.List;
import javax.swing.*;
import org.apache.commons.cli.*;
import processing.core.PApplet;
import utility.*;
import utility.HierarchicalClustering.ClusterDistFuncType;



public class Order3Tensor3DHeatMap  extends JFrame{
	
	
	private static final long serialVersionUID = 1L;

	public  Order3Tensor3DHeatMap (final PApplet C) throws Exception{
		super("Order3Tensor3DHeatMap");
		JScrollPane scr = new JScrollPane(C);
		scr.setPreferredSize(new Dimension(200, 100));
		getContentPane().add(scr, BorderLayout.CENTER);
		C.init();
		setLocation(100, 100);
		setSize(500, 500);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
	
	public static void main(String[] args) throws Exception {
		Options options = new Options();
		options.addOption("s", "supclust", false, "suppress clustering");
		options.addOption("o", "outclustbin", true,  "output clustered tensor binary file");
		options.addOption("b", "inclustbin", false,  "read infile as clustered tensor binary file");
		options.addOption("a", "rannot", true, "annot matrix files (split by ':')");
		options.addOption("A", "annotGmt", true, "annot matrix files in gmt format (split by ':')");
		options.addOption("S", "sym", false, "symetrical betweem order 1 and 2");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp("Expression [options] tensorfile", options);
			return;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 1){
			formatter.printHelp("Expression [options] tensorfile ", options);
			return;
		}
		ClusteredOrder3TensorWithAnnotation T;
		if(commandLine.hasOption("b")){
			T = new ClusteredOrder3TensorWithAnnotation(ClusteredOrder3Tensor.readFromBinary(argList.get(0)));
		}else{
			T = new ClusteredOrder3TensorWithAnnotation(argList.get(0));
			if(commandLine.hasOption("S")){
				T.setClusterDistFunc(ClusterDistFuncType.COMPLETE);
				T.supressOrder2Clustering();
				T.performClustering();
				T.copyClustering(1, 2);
			}else if(!commandLine.hasOption("s")){
			T.setClusterDistFunc(ClusterDistFuncType.COMPLETE);	
			  T.performClustering();
			}
		}
		if(commandLine.hasOption("o")){
			T.printAsBinary(commandLine.getOptionValue("o"));	
		}
		if(commandLine.hasOption("a")){
			List<String> tmp = MyFunc.split(":",commandLine.getOptionValue("a"));
			for(String s: tmp){
				T.setAnnotation(new StringMat(s));
			}
		}
		if(commandLine.hasOption("A")){
			List<String> tmp = MyFunc.split(":",commandLine.getOptionValue("A"));
			for(String s: tmp){
				T.setAnnotation(StringMat.getStringMatFromTxtFile(s));
			}
		}		
		Order3TensorViewer TV = new Order3TensorViewer(T);
		Order3Tensor3DHeatMap  H = new Order3Tensor3DHeatMap (TV);
		H.setVisible(true);		
	}
	
	
	
}