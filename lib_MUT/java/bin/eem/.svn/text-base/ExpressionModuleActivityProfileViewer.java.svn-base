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


public class ExpressionModuleActivityProfileViewer   extends JFrame   {
	
	private static final long serialVersionUID = -1320748718945288860L;



	public ExpressionModuleActivityProfileViewer(Applet C) throws Exception{
		setLayout(new BorderLayout());
		add(C,BorderLayout.CENTER);
		C.init();
		pack();
		setLocation(100, 100);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}

	public static void main(String[] args) throws Exception {
		Options options = new Options();
		options.addOption("p", "pcutoff", true, "p-value cutoff");
		options.addOption("o", "outfile", true, "output pdf file");
		options.addOption("a", "annot", true, "sample annotation file");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options]  eem_file exp_file", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 1  && argList.size() != 2){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options]  eem_file exp_file", options);
			return;
		}
		ExpressionModuleSet ems;
		if(argList.size() == 1){
			ems = ExpressionModuleSet.getFromMySQL(argList.get(0));
		}else{
			ems = ExpressionModuleSet. getFromTextFile(argList.get(0), new MyMat(argList.get(1)));	
		}
		if(commandLine.hasOption("p")){
			ems = ems.getSignificantSubset(Double.valueOf(commandLine.getOptionValue("p")));
		}	
		if(commandLine.hasOption("a")){
			ems.addSampleAnnotation(new StringMat(commandLine.getOptionValue("a")));
		}
		
		
		System.out.print(ems.getClusteredActivityProfiles());
		
		
		ClusteredMyMatViewer C = new ClusteredMyMatViewer(ems.getClusteredActivityProfiles());
		//ClusteredMyMatViewer C = new ClusteredMyMatViewer(ems.getClusteredBiclusterProfiles());
		C.scaleColorByRow();
		if(commandLine.hasOption("o")){
			C.setOutFile(commandLine.getOptionValue("o"));
		}
		ExpressionModuleActivityProfileViewer V =  new ExpressionModuleActivityProfileViewer(C);
		V.setVisible(true);
		
		
	}
	
	
}
