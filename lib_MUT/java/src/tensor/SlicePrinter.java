package tensor;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JScrollPane;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

import processing.core.PApplet;
import utility.MyFunc;
import utility.StringMat;
import utility.HierarchicalClustering.ClusterDistFuncType;

public class SlicePrinter extends JFrame{
	
	private static final long serialVersionUID = 1L;

	public  SlicePrinter (final PApplet C) throws Exception{
		super("SlicePrinter");
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
		options.addOption("S", "sliceord", true,  "sliced order");
		options.addOption("b", "inclustbin", false,  "read infile as clustered tensor binary file");
		options.addOption("k", "order1clust", true, "the number of order1 clusters");
		options.addOption("l", "order2clust", true, "the number of order2 clusters");
		options.addOption("m", "order3clust", true, "the number of order3 clusters");
		options.addOption("a", "annot", true, "annot matrix files (split by ':')");
		options.addOption("A", "annotGmt", true, "annot matrix files in gmt format (split by ':')");
		options.addOption("C", "complete", false, "perform complete clustering");
		options.addOption("w", "ward", false, "perform ward clustering");
		options.addOption("r", "red", false, "use red color scale");
		options.addOption("R", "annotred", false, "use red color scale for annot bar");
		options.addOption("B", "fixbox", false, "fix box size");
		options.addOption("O", "outfile", true, "output pdf file");
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
		int k = 0;
		int l = 0;
		int m = 0;
		if(commandLine.hasOption("b")){
			T = new ClusteredOrder3TensorWithAnnotation(ClusteredOrder3Tensor.readFromBinary(argList.get(0)));
		}else{
			T = new ClusteredOrder3TensorWithAnnotation(argList.get(0));
			if(!commandLine.hasOption("s")){	
				if(commandLine.hasOption("C")){
					T.setClusterDistFunc(ClusterDistFuncType.COMPLETE);
				}
				if(commandLine.hasOption("w")){
				   T.setClusterDistFunc(ClusterDistFuncType.WARD);
				}				
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
				T.setAnnotation(StringMat.getStringMatFromGmtFile(s));
			}
		}
		if(commandLine.hasOption("k")){
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
		}
		
		Order3TensorViewer TV = new Order3TensorViewer(T);
		
		if(commandLine.hasOption("S")){
			int i = Integer.valueOf(commandLine.getOptionValue("S"));
			if(i==1||i==2||i==3){
				TV.setSlicedOder(i);
				TV.calculateSize();
			}
		}
		if(commandLine.hasOption("r")){
			TV.setProfileColor("white", "red");
		}
		//if(commandLine.hasOption("R")){
		//	TV.setAnnotBandColor("white", "red");
		//}
		//T.print();
		if(commandLine.hasOption("B")){
			//TV.setBoxHeight(10);	
			//TV.setBoxWidth(10);
			//TV.calculateSize();
		}
		if(commandLine.hasOption("O")){
			TV.setOutFile(commandLine.getOptionValue("O"));
			TV.useAutoStop();
			TV.init();	
		}else{
			Order3TensorSliceHeatMap H = new Order3TensorSliceHeatMap(TV);
			H.setVisible(true);
		}
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
}
