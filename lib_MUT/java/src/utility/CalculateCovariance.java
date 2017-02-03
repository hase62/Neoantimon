package utility;

import java.util.List;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

import sun.reflect.Reflection;

public class CalculateCovariance {
	
	public static void main(String[] args) throws Exception {
		Options options = new Options();
		options.addOption("c", "col", false, "get cov for column vectors");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] mat", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 1){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] mat", options);
			return;
		}
		MyMat M = new MyMat(argList.get(0));
		if(commandLine.hasOption("c")){
			M = M.getCovMatForCol();
		}else{
			M = M.getCovMatForRow();
		}
		System.out.print(M);		
	}
}
