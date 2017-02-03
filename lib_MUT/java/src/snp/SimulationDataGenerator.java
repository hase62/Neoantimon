package snp;

import java.util.*;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.math.distribution.PoissonDistributionImpl;

import sun.reflect.Reflection;
import utility.MyFunc;
import utility.MyMat;

public class SimulationDataGenerator {
	protected int sampleNumber = 100;
	protected int probeNumber = 1000;
	protected List <Integer> concordantProbe; 
	protected MyMat D;
	protected SegmentContainerMap SCM;
	
	protected double concordantSampleRatio = 0.3;
	protected int concordantSegmentLength = 5;
	protected int nonconcordantSegmentLengthMean = 5;
	protected int nonconcordantSegmentNumberMean = 30;

	
	public void setConcordantProbe(String s){
		concordantProbe = new ArrayList<Integer>();
		for(String s2: MyFunc.split(",",s)){
			concordantProbe.add(Integer.valueOf(s2));
		}
	}
	
	
	public void setSampleNumber(int i){
		sampleNumber = i;
	}
	public void setProbeNumber(int i){
		probeNumber = i;
	}
	public void setNonconcordantSegmentLengthMean(double d){
		if(d < 1){
			nonconcordantSegmentLengthMean = (int)(d*probeNumber);
		}else{
			nonconcordantSegmentLengthMean = (int)d;
		}
	}
	
	public void setConcordantSampleRatio(double i){
		concordantSampleRatio = i;
	}
	public void setConcordantSegmentLength(int i){
		concordantSegmentLength = i;
	}
	public void setNonconcordantSegmentNumberMean(int i){
		nonconcordantSegmentNumberMean = i;
	}
	
	protected static int getPoisson(double lambda) {
		  double L = Math.exp(-lambda);
		  double p = 1.0;
		  int k = 0;

		  do {
		    k++;
		    p *= Math.random();
		  } while (p > L);

		  return k - 1;
	}
	
	protected static int getGeometric(double p){
		double u = Math.random();
		return (int)Math.floor(Math.log(1-u)/Math.log(1-p));
		
	}
	
	public SimulationDataGenerator(){
		concordantProbe = new ArrayList<Integer>();
		concordantProbe.add(200);
		concordantProbe.add(400);
		concordantProbe.add(600);
		concordantProbe.add(800);
		
		//concordantProbe.add(250);
		//concordantProbe.add(500);
		//concordantProbe.add(750);	
	}
	
	protected void calculateMatrix(){
		D = new MyMat(probeNumber, sampleNumber);
		
		for(int j = 0; j < sampleNumber; j++){
			
			for(int center: concordantProbe){
				center -= 1;
				if(Math.random() < concordantSampleRatio){
					int l = concordantSegmentLength;
					if(l>probeNumber){
						l=probeNumber;
					}
					if(l/2 == l/2.0){
						l += 1;
					}				
					int h = (l - 1)/2; 
					int start = center - h;
					int end = center + h;
					
					if(start < 0){
						start = 0;
					}
					if(end >= probeNumber){
						end = probeNumber-1;
					}
					for(int i = start; i <=end; i++){
						D.set(i,j,1);
					}	
				}
			}
			int NonConcordantSegmentNumber = getPoisson(nonconcordantSegmentNumberMean);
			int k = 0;
			L:while(k < NonConcordantSegmentNumber){
				int l = getGeometric(1.0/nonconcordantSegmentLengthMean);
				if(l>probeNumber){
					l=probeNumber;
				}
				if(l/2 != l/2.0){
					l += 1;
				}				
				int h = (l - 1)/2; 
				int center = (int)(Math.random()*probeNumber);
				int start = center - h;
				int end = center + h;
				if(start < 0){
					start = 0;
				}
				if(end >= probeNumber){
					end = probeNumber-1;
				}
				//for(int i = start; i <= end; i++){
					//if( D.get(i, j) == 1){
						//continue L;
					//}
				//}
				for(int i = start; i <=end; i++){
					D.set(i,j,1);
				}
				k++;
			}
		}
	}
	
	protected void calculateSegmentContainerMap(){
		SCM = calculateSegmentContainerMap(D);	
	}
	
	protected static SegmentContainerMap calculateSegmentContainerMap(MyMat D){
		SegmentContainerMap SCM = new SegmentContainerMap();
		for(int j = 0; j < D.colSize(); j++){
			double preValue = D.get(0,j);
			int start = 0;
			int	 end = 0;
			SegmentContainer SC = new SegmentContainer();
			for(int i = 0; i < D.rowSize(); i++){
				if(preValue != D.get(i,j)){
					end=i-1;
					SC.add(new Segment(1,start+1,end+1,end-start+1,preValue));
					preValue=D.get(i,j);
					start=i;
				}
			}
			end=D.rowSize()-1;
			SC.add(new Segment(1,start+1,end+1,end-start+1,preValue));
			SCM.put("sample" + (j+1), SC);
		}
		return SCM;
	}

	protected void generateData(){
		calculateMatrix();
		calculateSegmentContainerMap();
	}
	
	public static void main(String [] args) throws Exception{Options options = new Options();
	
		SimulationDataGenerator S = new SimulationDataGenerator();
		options.addOption("r", "conrario", true, "concordant sample ratio (" + S.concordantSampleRatio + ")");
		options.addOption("c", "conpos", true, "concordant  position (" + MyFunc.join(",",S.concordantProbe)+ ")");
		options.addOption("s", "lnonc", true, "nonconcordant segment length average (" + S.nonconcordantSegmentLengthMean + ")");
		options.addOption("n", "nnoc", true, "nonconcordant segment number avarage(" + S.nonconcordantSegmentNumberMean + ")");
		options.addOption("S", "sampsize", true, "sample size (" + S.sampleNumber + ")");
		options.addOption("P", "prbsize", true, "probe size (" + S.probeNumber + ")");
		
	
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName(), options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 0){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName(), options);
			return;
		}
	
		if(commandLine.hasOption("r")){
			S.setConcordantSampleRatio(Double.valueOf(commandLine.getOptionValue("r")));
		}
		if(commandLine.hasOption("s")){
			S.setNonconcordantSegmentLengthMean(Integer.valueOf(commandLine.getOptionValue("s")));
		}
		if(commandLine.hasOption("n")){
			S.setNonconcordantSegmentNumberMean(Integer.valueOf(commandLine.getOptionValue("n")));
		}
		if(commandLine.hasOption("S")){
			S.setSampleNumber(Integer.valueOf(commandLine.getOptionValue("S")));
		}
		if(commandLine.hasOption("P")){
			S.setProbeNumber(Integer.valueOf(commandLine.getOptionValue("P")));
		}
		if(commandLine.hasOption("c")){
				S.setConcordantProbe(commandLine.getOptionValue("c"));
		}
		S.generateData();
		S.SCM.print();
		//S.D.print();
	}
}
