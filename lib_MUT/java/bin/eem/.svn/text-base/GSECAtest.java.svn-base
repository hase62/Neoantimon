package eem;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

import no.uib.cipr.matrix.NotConvergedException;
import sun.reflect.Reflection;
import utility.MyFunc;
import utility.MyMat;

public class GSECAtest extends GeneSetExpressionCoherenceAnalysis {
	
	
	
	
	int geneSetSize = 100;
	int n = 100;
	
	
	public GSECAtest(MyMat E) {
		super(E);
	}

	
	public void perform(){
		
		if(geneSetSize >= Exp.rowSize()){
			n = 1;
			geneSetSize = Exp.rowSize();
		}
	
		
		calculateCor();
		List <Double> corMean = new ArrayList<Double>();
		for(int i = 0; i < n ; i++){
			List <String> gs = MyFunc.sample(Exp.getRowNames(), geneSetSize);
			
				corMean.add(getCorMean(gs));
			
		}
		
	
		System.out.println(MyFunc.join("\n", corMean));	
		
	
	}
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("n", "n", true,  "sample number");
		options.addOption("g", "gs", true, "gene set size");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] expfile", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 1){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] expfile", options);
			return;
		}
		
	
		GSECAtest S = new GSECAtest(new MyMat(argList.get(0)));
		
		
		if(commandLine.hasOption("n")){
			S.n = Integer.valueOf(commandLine.getOptionValue("n"));
		}
		if(commandLine.hasOption("g")){
			S.geneSetSize = Integer.valueOf(commandLine.getOptionValue("g"));
		}
		
		S.perform();
	
	}
	

}
