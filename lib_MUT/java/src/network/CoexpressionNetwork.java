package network;

import java.io.IOException;
import java.util.*;

import org.apache.commons.cli.*;

import sun.reflect.Reflection;
import utility.*;

public class CoexpressionNetwork {
	
	private MyMat E;
	private Dist C;
	private Link L;
	private double cutoff = 0;

	public CoexpressionNetwork(MyMat E){
		this.E = new MyMat(E);
		this.E.normalizeRows();
	}
	public void getCor(){
		C = new Dist(E, 'C');
	}
	
	
	public void convert2FisherZ(){
		int i, j;
		int n = C.size();
		double max = 0.999999999999;
		for(i=0; i<n; i++){
			for(j=0;j<i;j++){
				double  c = C.get(i, j);	
				if(c  > max){
					c = max;
				}
				if(c  < -max){
					c = -max;
				}
				C.set(i, j,0.5 * Math.log((1+c)/(1-c))); 
			}
		}
	}
	
	public void normalizeCor(){
		int i, j;
		int n = C.size();
		List<Double> tmp = new ArrayList<Double>();
		for(i=0; i<n; i++){
			for(j=0;j<i;j++){
				tmp.add(C.get(i, j));
			}
		}
		double mean = MyFunc.mean(tmp);
		double sd = MyFunc.sd(tmp);	
		for(i=0; i<n; i++){
			for(j=0;j<i;j++){
				C.set(i, j, (C.get(i, j)-mean)/sd); 
			}
		}	

	}
	
	public void print() throws IOException{
		C.print();
	}
	
	public void getLink(){
		List <String> gene = E.getRowNames();
		L = new Link(gene);
		for(int i = 0; i < gene.size(); i++){
			for(int j = 0; j < i; j++){
				if(C.get(gene.get(i), gene.get(j))>= cutoff){
					L.set(gene.get(i), gene.get(j), true);
				}
			}
		}
	}
	
	public void printResult() throws IOException{
		L.printDataInSifFormat("coexp");
	}
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("c", "cutoff", true,  "correlation cutoff");
		options.addOption("n", "normalize", false,  "normalize");
		options.addOption("f", "fisherz", false,  "convert fisherz");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options]  tabFile", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 1){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] tabFile", options);
			return;
		}
		CoexpressionNetwork C = new CoexpressionNetwork(new MyMat(argList.get(0)));
		if(commandLine.hasOption("c")){
			C.cutoff = Double.valueOf(commandLine.getOptionValue("c"));
		}
		
		C.getCor();
		if(commandLine.hasOption("f")){
			C.convert2FisherZ();
		}
		if(commandLine.hasOption("n")){
			C.normalizeCor();
		}
		if(!commandLine.hasOption("c")){
			C.print();
		}else{
			C.getLink();
			C.printResult();
		}
	}
}
