package utility;

import java.applet.Applet;
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

import sun.reflect.Reflection;
import utility.HierarchicalClustering.ClusterDistFuncType;
import utility.HierarchicalClustering.DistFuncType;

public class Heatmap  extends JFrame  {
	

	private static final long serialVersionUID = 8212970936306216242L;

	public  Heatmap (Applet C) throws Exception{		
		super("Heatmap");
		JScrollPane scr = new JScrollPane(C);
		scr.setPreferredSize(new Dimension(200, 100));
		getContentPane().add(scr, BorderLayout.CENTER);
		C.init();
		setLocation(100, 100);
		setSize(1000, 1000);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("n", "norclst", false, "suppress row clustering");
		options.addOption("N", "nocclst", false, "suppres col clustering");
		options.addOption("s", "rscl", false, "scale color by row");
		options.addOption("S", "cscl", false, "scale color by column");
		options.addOption("s", "rscl", false, "scale color by row");
		options.addOption("S", "cscl", false, "scale color by column");
		options.addOption("t", "rsrt", true, "sort rows");
		options.addOption("T", "csrt", true, "sort columns");
		options.addOption("a", "annot", true, "annot matrix files (split by ':')");
		options.addOption("G", "annotGmt", true, "annot matrix files from gmt file (split by ':')");
		options.addOption("A", "annotTxt", true, "annot matrix files from txt file (split by ':')");
		options.addOption("c", "cor", false, "use correlation");
		options.addOption("C", "complete", false, "perform complete clustering");
		options.addOption("r", "red", false, "use red color scale");
		options.addOption("R", "annotred", false, "use red color scale for annot bar");
		options.addOption("w", "ward", false, "perform ward clustering");
		//options.addOption("k", "rclust", true, "the number of row clusters");
		//options.addOption("l", "cclust", true, "the number of col clusters");
		options.addOption("F", "fixbox", false, "fix box size");
		options.addOption("W", "writepdf", true, "write pdf file");
		options.addOption("b", "writebinary", true, "write binary file");
		options.addOption("B", "readbinary", false, "read binary file");
		options.addOption("E", "wodraw", false, "exit without drawing");
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
		
		
		ClusteredMyMatWithAnnotation M = null;
		if(commandLine.hasOption("B")){
			M = ClusteredMyMatWithAnnotation.readFromBinary(argList.get(0));
		}else{
		M = new ClusteredMyMatWithAnnotation(argList.get(0));
		
		if(commandLine.hasOption("n")){
			M.supressRowClustering();	
		}	
		if(commandLine.hasOption("N")){
			M.supressColClustering();	
		}
		if(commandLine.hasOption("c")){
			M.setRowDistFunc(DistFuncType.CORRELATION);
			M.setColDistFunc(DistFuncType.CORRELATION);
		}	
	
		if(commandLine.hasOption("C")){
			M.setRowClusterDistFunc(ClusterDistFuncType.COMPLETE);
			M.setColClusterDistFunc(ClusterDistFuncType.COMPLETE);
		}
		if(commandLine.hasOption("w")){
			M.setRowClusterDistFunc(ClusterDistFuncType.WARD);
			M.setColClusterDistFunc(ClusterDistFuncType.WARD);
		}
		if(commandLine.hasOption("t")){
			M.sortRowsByValue(commandLine.getOptionValue("t"));
			M.supressRowClustering();	
		}	
		if(commandLine.hasOption("T")){
			M.sortColsByValue(commandLine.getOptionValue("T"));
			M.supressColClustering();	
		}	
		M.performClustering();
		}
			
		/*int k = 0;
		int l = 0;
		if(commandLine.hasOption("k")){
			k = Integer.valueOf(commandLine.getOptionValue("k"));
		}
		if(commandLine.hasOption("l")){
			l = Integer.valueOf(commandLine.getOptionValue("l"));
		}*/
		
		if(commandLine.hasOption("a")){
			List<String> tmp = MyFunc.split(":",commandLine.getOptionValue("a"));
			for(String s: tmp){
				M.addAnnotation(new StringMat(s));
			}
		}
		if(commandLine.hasOption("A")){
			List<String> tmp = MyFunc.split(":",commandLine.getOptionValue("A"));
			for(String s: tmp){
				M.addAnnotation(StringMat.getStringMatFromTxtFile(s));
			}
		}
		if(commandLine.hasOption("G")){
			List<String> tmp = MyFunc.split(":",commandLine.getOptionValue("G"));
			for(String s: tmp){
				M.addAnnotation(StringMat.getStringMatFromGmtFile(s));
			}
		}
		
		/*if(k > 1){
			M.addAnnotation(new StringMat("rowcluster", M.getRowClustering().getCutTreeMap(k)));	
		}
		if(l > 1){
			M.addAnnotation(new StringMat("colcluster", M.getColClustering().getCutTreeMap(l)));
		}*/
		
		if(commandLine.hasOption("b")){
			M.printAsBinary(commandLine.getOptionValue("b"));
		}
		if(commandLine.hasOption("E")){
			System.exit(0);
		}
		
		
		if(commandLine.hasOption("W")){
			 ClusteredMyMatViewer MV = new ClusteredMyMatViewer(M);
			 if(commandLine.hasOption("s")){
				 MV.scaleColorByRow();
			 }
			 if(commandLine.hasOption("S")){
				 MV.scaleColorByColumn();
			 }
			 if(commandLine.hasOption("r")){
				 MV.setProfileColor("white", "red");
			}
			if(commandLine.hasOption("R")){
				MV.setAnnotBandColor("white", "red");
			}
			
			if(commandLine.hasOption("F")){
				MV.setBoxHeight(10);	
				MV.setBoxWidth(10);
				MV.calculateSize();
			}
			MV.setOutFile(commandLine.getOptionValue("W"));
			MV.useAutoStop();
			Heatmap H = new Heatmap(MV);
			H.setVisible(true);
		}else{	
			InteractiveClusteredMyMatViewer  MV = new InteractiveClusteredMyMatViewer(M);
		
			if(commandLine.hasOption("s")){
				MV.scaleColorByRow();
			}
			if(commandLine.hasOption("S")){
				MV.scaleColorByColumn();
			}
			if(commandLine.hasOption("r")){
				MV.setProfileColor("white", "red");
			}
			if(commandLine.hasOption("R")){
				MV.setAnnotBandColor("white", "red");
			}
			if(commandLine.hasOption("F")){
				MV.setBoxHeight(10);	
				MV.setBoxWidth(10);
				MV.calculateSize();
			}
			Heatmap H = new Heatmap(MV);
			H.setVisible(true);
		}
	}
}
