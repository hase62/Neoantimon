package snp;

import java.applet.Applet;
import java.awt.*;
import java.util.*;
import java.util.List;
import org.apache.commons.cli.*;
import javax.swing.*;
import utility.*;
import utility.HierarchicalClustering.ClusterDistFuncType;

import sun.reflect.Reflection;

public class ChromHeatmap extends JFrame{

	private static final long serialVersionUID = 1L;

	public ChromHeatmap (Applet C) throws Exception{		
		super("ChromHeatmap");
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
		options.addOption("p", "probe", true, "probe coordinate file");
		options.addOption("W", "writepdf", true, "output pdf file");
		options.addOption("g", "gene", true, "gene coordinate file");
		options.addOption("t", "thinpi", true, "thin down probe info");
		options.addOption("n", "nump", true, "# of psuedo probes");
		options.addOption("l", "locus", true, "locus (e.g., 3:120000-130000)");
		options.addOption("r", "rmxy", false, "remove XY");
		options.addOption("s", "csrt", true, "sort columns");
		options.addOption("a", "annot", true, "annot matrix files (split by ':')");
		options.addOption("A", "annotTxt", true, "annot matrix files from txt file (split by ':')");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] SegFile probeInfoFile", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 1){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] SegFile probeInfoFile", options);
			return;
		}
		SegmentContainerMap SCM = new SegmentContainerMap(argList.get(0));
		
		ProbeInfo PI = null;
		List<String> sample = new ArrayList<String>(SCM.keySet());
		if(commandLine.hasOption("p")){
			PI = 	ProbeInfo.getProbeInfoFromTsvFile(commandLine.getOptionValue("p"));
			if(commandLine.hasOption("t")){
				PI.thinDown(Integer.valueOf(commandLine.getOptionValue("t")));
			}
			PI.filter(SCM.get(sample.get(0)));
		}
		if(commandLine.hasOption("g")){
			GeneInfo GI = GeneInfo.getGeneInfoFromTsvFile(commandLine.getOptionValue("g"));
			PI = GI.toProbeInfo();
		}
		if(commandLine.hasOption("l")){
			List <String> tmp = MyFunc.split(":", commandLine.getOptionValue("l"));
			int chr = Integer.valueOf(tmp.get(0));
			if(tmp.size() == 1){
				for(SegmentContainer SC: SCM.values()){
					SC.filter(chr);
				}
			}else{
				List <String> tmp2 = MyFunc.split("-", tmp.get(1));
				if(tmp2.size()==2){
					int start = Integer.valueOf(tmp2.get(0));
					int end = Integer.valueOf(tmp2.get(1));	
					for(SegmentContainer SC: SCM.values()){
						SC.filter(chr,start,end);
					}	
				}
			}
		}
		if(commandLine.hasOption("r")){
			SCM.removeXY();
		}

		if(commandLine.hasOption("n")){
			int i = Integer.valueOf(commandLine.getOptionValue("n"));
			PI = SCM.generatePsuedoProbeInfo(i);
		}
		
		if(PI == null){
			PI = SCM.generatePsuedoProbeInfo(1000);
		}
		
		PI.filter(SCM.get(sample.get(0)));
		
		
		ClusteredMyMatWithAnnotation M = new ClusteredMyMatWithAnnotation(SCM.toMyMat(PI));
		
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
		
		if(commandLine.hasOption("s")){
			M.sortColsByValue(commandLine.getOptionValue("s"));
			M.supressColClustering();	
		}	
		
		
		if(SCM.get(sample.get(0)).chrList().size() > 1){
			Map <String, String> tmp = new HashMap<String, String>();
			for(String p: PI.getProbeSet()){
				int c = PI.chr(p);
				c = c-(c/2)*2;
				tmp.put(p,Integer.valueOf(c).toString());
			}
			M.addAnnotation(new StringMat("colcluster", tmp));
		}
		M.supressRowClustering();
		if(M.colSize() > 1000){
			M.supressColClustering();
		}
		M.setColClusterDistFunc(ClusterDistFuncType.WARD);
		M.performClustering();
		
		if(commandLine.hasOption("W")){
			ClusteredMyMatViewer MV = new ClusteredMyMatViewer(M);
			MV.setOutFile(commandLine.getOptionValue("W"));
			MV.useAutoStop();
			Heatmap H = new Heatmap(MV);
			H.setVisible(true);
		}else{
			InteractiveClusteredMyMatViewer  MV = new InteractiveClusteredMyMatViewer(M);
			Heatmap H = new Heatmap(MV);
			H.setVisible(true);
		}
	
	
	}

}
