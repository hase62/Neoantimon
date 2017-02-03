package eem;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

import sun.reflect.Reflection;
import utility.*;

public class ExpressionModuleOverlapTest {
	public static void main(String[] args) throws Exception {
		Options options = new Options();
		options.addOption("p", "pcutoff", true, "expression module p-value cutoff (-log scale)");
		options.addOption("b", "bg", true, "# of background genes");
		options.addOption("P", "pcutoff", true, "overlap p-value cutoff  (-log scale)");
		options.addOption("Q", "pcutoff", true, "overlap q-value cutoff  (-log scale)");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
	
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] table_name gmt_file", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 2 && argList.size() != 1){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] table_name gmt_file", options);
			return;
		}
		ExpressionModuleSet ems = ExpressionModuleSet.getFromMySQL(argList.get(0));
		if(commandLine.hasOption("p")){
			ems = ems.getSignificantSubset(Double.valueOf(commandLine.getOptionValue("p")));
		}
		
		GeneSetOverlapTest GSO;
		if(argList.size() != 1){
			GSO = new GeneSetOverlapTest(ems.getModuleGenes(), MyFunc.readGeneSetFromGmtFile(argList.get(1)));
		}else{
			GSO = new GeneSetOverlapTest(ems.getModuleGenes());
		}
		GSO.outputInMinusLogScale();
		if(commandLine.hasOption("b")){
			GSO.setBgGeneNumber(Integer.valueOf(commandLine.getOptionValue("b")));
		}
		if(commandLine.hasOption("P")){
			GSO.setPvalueCutoff(Double.valueOf(commandLine.getOptionValue("P")));
		}
		if(commandLine.hasOption("Q")){
			GSO.setQvalueCutoff(Double.valueOf(commandLine.getOptionValue("Q")));
		}
		GSO.calculatePvalue();
		GSO.calculateQvalue();
		System.out.print(GSO);
	}		
	
}	
