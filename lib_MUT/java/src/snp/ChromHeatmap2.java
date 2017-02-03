package snp;

import java.applet.Applet;
import java.awt.*;
import java.util.*;
import java.util.List;
import org.apache.commons.cli.*;
import javax.swing.*;

import tensor.ClusteredOrder3TensorViewer2;
import tensor.ClusteredOrder3TensorWithAnnotation;
import utility.*;
import utility.HierarchicalClustering.ClusterDistFuncType;

import sun.reflect.Reflection;

public class ChromHeatmap2 extends JFrame{

	private static final long serialVersionUID = 1L;

	public ChromHeatmap2 (Applet C) throws Exception{		
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
		options.addOption("o", "out", true, "output pdf file");
		options.addOption("g", "gene", true, "gene coordinate file");
		options.addOption("t", "thinpi", true, "thin down probe info");
		options.addOption("n", "nump", true, "# of psuedo probes");
		options.addOption("l", "locus", true, "locus (e.g., 3:120000-130000)");
		options.addOption("r", "rmxy", false, "remove XY");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] SegFile1 SegFile2", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() <= 1){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] SegFile1 SegFile2", options);
			return;
		}
		
		List <SegmentContainerMap> SCMList = new ArrayList<SegmentContainerMap>();
		Map <String, Integer> seen = new HashMap <String, Integer>();
		
		for(int i = 0; i< argList.size(); i++){
			SCMList.add(new SegmentContainerMap(argList.get(i)));
			for(String s: SCMList.get(i).keySet()){
				if(!seen.containsKey(s)){
					seen.put(s, 1);
				}else{
					seen.put(s, seen.get(s)+1);
				}				
			}
		}
		if(seen.isEmpty()){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] SegFile1 SegFile2", options);
			return;
		}
		List <String> sample  = new ArrayList<String>();
		for(String s: seen.keySet()){
			if(seen.get(s) == SCMList.size()){
				sample.add(s);
			}
		}
		for(SegmentContainerMap S: SCMList){
			S = S.getSubMap(sample);
		}
		
		SegmentContainerMap SCM0 = SCMList.get(0);
	
		ProbeInfo PI = null;
		
		if(commandLine.hasOption("p")){
			PI = 	ProbeInfo.getProbeInfoFromTsvFile(commandLine.getOptionValue("p"));
			if(commandLine.hasOption("t")){
				PI.thinDown(Integer.valueOf(commandLine.getOptionValue("t")));
			}
			PI.filter(SCM0.get(sample.get(0)));
		}
		if(commandLine.hasOption("g")){
			GeneInfo GI = GeneInfo.getGeneInfoFromTsvFile(commandLine.getOptionValue("g"));
			PI = GI.toProbeInfo();
		}
		if(commandLine.hasOption("l")){
			List <String> tmp = MyFunc.split(":", commandLine.getOptionValue("l"));
			int chr = Integer.valueOf(tmp.get(0));
			if(tmp.size() == 1){
				for(SegmentContainerMap SCM: SCMList){
					for(SegmentContainer SC: SCM.values()){
						SC.filter(chr);
					}
				}
			}else{
				List <String> tmp2 = MyFunc.split("-", tmp.get(1));
				if(tmp2.size()==2){
					int start = Integer.valueOf(tmp2.get(0));
					int end = Integer.valueOf(tmp2.get(1));
					for(SegmentContainerMap SCM: SCMList){
						for(SegmentContainer SC: SCM.values()){
							SC.filter(chr,start,end);
						}	
					}
				}
			}
		}
		if(commandLine.hasOption("r")){
			for(SegmentContainerMap SCM: SCMList){
				SCM.removeXY();
			}
		}

		if(commandLine.hasOption("n")){
			int interval = Integer.valueOf(commandLine.getOptionValue("n"));
			PI = ProbeInfo.generatePsuedoProbeInfo(interval);
			for(SegmentContainerMap SCM: SCMList){
				PI.filter(SCM.get(sample.get(0)));
			}
		}
		
		if(PI == null){
			PI = ProbeInfo.generatePsuedoProbeInfo(1000);
			PI.filter(SCM0.get(sample.get(0)));
		}
		
		PI.filter(SCM0.get(sample.get(0)));
		
		
		List <MyMat> M = new ArrayList<MyMat>();
		for(SegmentContainerMap S: SCMList){
			M.add(S.toMyMat(PI));
		}
		
		
		ClusteredOrder3TensorWithAnnotation T = new ClusteredOrder3TensorWithAnnotation(M);
		
		if(SCM0.get(sample.get(0)).chrList().size() > 1){
			Map <String, String> tmp = new HashMap<String, String>();
			for(String p: PI.getProbeSet()){
				int c = PI.chr(p);
				c = c-(c/2)*2;
				tmp.put(p,Integer.valueOf(c).toString());
			}
			T.addAnnotation(new StringMat("colcluster", tmp));
		}
		T.reorderOrder1(M.get(0).getRowNames());
		T.setName3(argList);
		T.supressOrder1Clustering();
		T.supressOrder3Clustering();
		if(T.getDimOfOrder2()> 1000){
			T.supressOrder2Clustering();
		}
		T.setOrder2ClusterDistFunc(ClusterDistFuncType.WARD);
		T.performClustering();
		ClusteredOrder3TensorViewer2 TV = (ClusteredOrder3TensorViewer2) T.getViewer();
		if(commandLine.hasOption("o")){
			TV.setPdfFileName(commandLine.getOptionValue("o"));
		}
		Heatmap H = new Heatmap(TV);
		H.setVisible(true);
	}
}
