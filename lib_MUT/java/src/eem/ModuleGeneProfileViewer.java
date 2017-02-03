package eem;


import sun.reflect.Reflection;
import utility.*;

import java.applet.Applet;
import java.awt.BorderLayout;
import java.util.List;

import javax.swing.JFrame;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;


public class  ModuleGeneProfileViewer   {
	
	private static final long serialVersionUID = -1320748718945288860L;


	public static void main(String[] args) throws Exception {
		Options options = new Options();
		options.addOption("o", "outfile", true, "output pdf file");
		options.addOption("s", "seed", false, "show seed genes");
		options.addOption("a", "annot", true, "sample annotation file");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] eem_file exp_file module_id    ", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 2  && argList.size() != 3 ){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] eem_file exp_file module_id ", options);
			return;
		}
		
		ExpressionModuleSet ems;
		ExpressionModule em;
		if(argList.size() == 2 ){
			ems = ExpressionModuleSet.getFromMySQL(argList.get(0));
			em = ems.get(argList.get(1));
		}else{
			ems = ExpressionModuleSet. getFromTextFile(argList.get(0), new MyMat(argList.get(1)));
			em = ems.get(argList.get(2));
		}
		em.showBiclusteredSeedGeneProfile();
		//em.addBiclusterInformation2SampleAnnotation();
		ClusteredMyMatViewer C;
		if(commandLine.hasOption("a")){
			em.addSampleAnnotation(new StringMat(commandLine.getOptionValue("a")));
		}
		if(commandLine.hasOption("s")){
			C = new ClusteredMyMatViewer(em.getClusteredSeedGeneProfiles());
		}else{
			C = new ClusteredMyMatViewer(em.getClusteredModuleGeneProfiles());
		}
		C.scaleColorByRow();
		if(commandLine.hasOption("o")){
			C.setOutFile(commandLine.getOptionValue("o"));
			C.useAutoStop();
			C.init();
		}else{
			Heatmap V =  new Heatmap(C);
			V.setVisible(true);
		}
	}
}


