package tensor;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.util.List;
import javax.swing.*;
import org.apache.commons.cli.*;

import processing.core.PApplet;
import utility.HierarchicalClustering.ClusterDistFuncType;
import utility.*;

public class Order3TensorSliceHeatMap  extends JFrame{
	
	
	private static final long serialVersionUID = 1L;

	public  Order3TensorSliceHeatMap (final PApplet C) throws Exception{
		super("Order3TensorSliceHeatMap");
		JScrollPane scr = new JScrollPane(C);
		scr.setPreferredSize(new Dimension(200, 100));
		getContentPane().add(scr, BorderLayout.CENTER);
		C.init();
		setLocation(100, 100);
		setSize(1000, 1000);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
	}

	public static void main(String[] args) throws Exception {
		Options options = new Options();
		options.addOption("s", "supclust", false, "suppress clustering");
		options.addOption("b", "writebinary", true, "write binary file");
		options.addOption("B", "readbinary", false, "read binary file");
		//options.addOption("S", "sliceord", true,  "sliced order");
		//options.addOption("k", "order1clust", true, "the number of order1 clusters");
		//options.addOption("l", "order2clust", true, "the number of order2 clusters");
		//options.addOption("m", "order3clust", true, "the number of order3 clusters");
		options.addOption("a", "annot", true, "annot matrix files (split by ':')");
		//options.addOption("A", "annotGmt", true, "annot matrix files in gmt format (split by ':')");
		options.addOption("C", "complete", false, "perform complete clustering");
		options.addOption("w", "ward", false, "perform ward clustering");
		options.addOption("r", "red", false, "use red color scale");
		//options.addOption("R", "annotred", false, "use red color scale for annot bar");
		//options.addOption("F", "fixbox", false, "fix box size");
		//options.addOption("E", "wodraw", false, "exit without drawing");
		options.addOption("y", "sym", false, "symetrical betweem order 1 and 2");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp("MCHM [options] tensorfile", options);
			return;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 1){
			formatter.printHelp("MCHM  [options] tensorfile ", options);
			return;
		}
		ClusteredOrder3TensorWithAnnotation T;
		//int k = 0;
		//int l = 0;
		//int m = 0;
		if(commandLine.hasOption("B")){
			T =  ClusteredOrder3TensorWithAnnotation.readFromBinary(argList.get(0));
		}else{
			//T = new 	ClusteredOrder3TensorWithAnnotation(argList.get(0));
			T = ClusteredOrder3TensorWithAnnotation.getBigClusteredOrder3TensorWithAnnotation(argList.get(0));
			if(!commandLine.hasOption("s")){	
				if(commandLine.hasOption("C")){
					T.setClusterDistFunc(ClusterDistFuncType.COMPLETE);
				}
				if(commandLine.hasOption("w")){
				   T.setClusterDistFunc(ClusterDistFuncType.WARD);
				}				
				if(commandLine.hasOption("y")){
					T.supressOrder2Clustering();
				}
				T.performClustering();
				if(commandLine.hasOption("y")){
					T.copyClustering(1, 2);
				}
			}
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
		
		if(commandLine.hasOption("b")){
			T.printAsBinary(commandLine.getOptionValue("b"));	
		}
		if(commandLine.hasOption("E")){
			System.exit(0);
		}
		
		/*if(commandLine.hasOption("k")){
			k = Integer.valueOf(commandLine.getOptionValue("k"));
		}
		if(commandLine.hasOption("l")){
			l = Integer.valueOf(commandLine.getOptionValue("l"));
		}
		if(commandLine.hasOption("m")){
			m = Integer.valueOf(commandLine.getOptionValue("m"));
		}
		if(k > 1){
			T.addAnnotation(new StringMat("order1cluster", T.getOrder1Clustering().getCutTreeMap(k)));	
		}
		if(l > 1){
			T.addAnnotation(new StringMat("order2cluster", T.getOrder2Clustering().getCutTreeMap(l)));
		}
		if(m > 1){
			T.addAnnotation(new StringMat("order3cluster", T.getOrder3Clustering().getCutTreeMap(m)));	
		}*/
		
		//T.compressOrder1(4);
		
		ClusteredOrder3TensorViewer2 TV = (ClusteredOrder3TensorViewer2) T.getViewer();
		
		//if(commandLine.hasOption("S")){
		//	int i = Integer.valueOf(commandLine.getOptionValue("S"));
		//	if(i==1||i==2||i==3){
		//		TV.setSlicedOder(i);
		//		TV.calculateSize();
		//	}
		//}
		if(commandLine.hasOption("r")){
			TV.setProfileColor("white", "red");
		}
		//if(commandLine.hasOption("R")){
		//	TV.setAnnotBandColor("white", "red");
		//}
		//T.print();
		//if(commandLine.hasOption("F")){
		//	TV.setBoxHeight(10);	
		//	TV.setBoxWidth(10);
		//	TV.calculateSize();
		//}
		Order3TensorSliceHeatMap H = new Order3TensorSliceHeatMap(TV);
		H.setVisible(true);
		
	}
	
	
}
