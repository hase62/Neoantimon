package snp;

import java.io.IOException;
import java.util.*;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

import sun.reflect.Reflection;
import utility.MyFunc;
import utility.MyMat;


public class SimulationDataGeneratorDecomposed extends SimulationDataGenerator{
	
	protected MyMat R;
	protected MyMat Rc;
	protected MyMat Rn;
	protected MyMat Rf;
	protected SegmentContainerMap SCM;
	protected SegmentContainerMap SCMc;
	protected SegmentContainerMap SCMn;
	protected SegmentContainerMap SCMf;
	
	
	List <Integer> concordantProbe2; 
	double concordantSampleRatio2 = 0.3;
	int concordantSegmentLength2 = 5;
	int concordantSegmentNumber2 = 10;
	
	
	public void setConcordantProbe2(String s){
		concordantProbe2 = new ArrayList<Integer>();
		for(String s2: MyFunc.split(",",s)){
			concordantProbe2.add(Integer.valueOf(s2));
		}
	}
	
	
	public void setConcordantSampleRatio2(double i){
		concordantSampleRatio2 = i;
	}
	public void setConcordantSegmentLength2(int i){
		concordantSegmentLength2 = i;
	}
	
	public void setConcordantSegmentNumber2(int k){
		concordantSegmentNumber2 = k;
		concordantProbe2 = new ArrayList<Integer>();
		for(int i=0; i< concordantSegmentNumber2; i++){
			concordantProbe2.add((int)(Math.random()*probeNumber));
		}
	}
	
	
	public SimulationDataGeneratorDecomposed(){
		super();
		setConcordantSegmentNumber2(concordantSegmentNumber2);
		//concordantProbe2.add(100);
		//concordantProbe2.add(300);
		//concordantProbe2.add(500);
		//concordantProbe2.add(700);	
		//concordantProbe2.add(900);	
	}
	
	protected void calculateMatrixRc(){
		Rc = new MyMat(probeNumber, sampleNumber);
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
						Rc.set(i,j,1);
					}	
				}
			}
		}
	}
	protected void calculateMatrixRf(){
		Rf = new MyMat(probeNumber, sampleNumber);
			for(int j = 0; j < sampleNumber; j++){
			for(int center: concordantProbe2){
				center -= 1;
				if(Math.random() < concordantSampleRatio2){
					int l = concordantSegmentLength2;
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
						Rf.set(i,j,1);
					}	
				}
			}
			}
	}
	protected void calculateMatrixRn(){
		Rn = new MyMat(probeNumber, sampleNumber);
			for(int j = 0; j < sampleNumber; j++){
			int aberrationNonConcordantSegmentNumber = SimulationDataGenerator.getPoisson(nonconcordantSegmentNumberMean);
			int k = 0;
			L:while(k < aberrationNonConcordantSegmentNumber){
				int l = SimulationDataGenerator.getGeometric(1.0/nonconcordantSegmentLengthMean);
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
					Rn.set(i,j,1);
				}
				k++;
			}
		}
	}
	protected void calculateMatrixR(){
		R = new MyMat(probeNumber, sampleNumber);
		for(int j = 0; j < sampleNumber; j++){
			for(int i = 0; i < probeNumber; i++){
				if(Rc.get(i,j)+Rf.get(i, j)+Rn.get(i,j) > 0){
					R.set(i, j, 1);
				}
			}
		}
	}
	

	protected void calculateMatrix(){
		calculateMatrixRc();
		calculateMatrixRn();
		calculateMatrixRf();
		calculateMatrixR();
	}
	
	protected void calculateSegmentContainerMap(){
		SCM = calculateSegmentContainerMap(R);
		SCMc = calculateSegmentContainerMap(Rc);	
		SCMn = calculateSegmentContainerMap(Rn);	
		SCMf = calculateSegmentContainerMap(Rf);	
	}
	protected void generateData(){
		calculateMatrix();
		calculateSegmentContainerMap();
	}
	
	
	protected void print() throws IOException{
		//SCM.print("R.seg");
		//SCMc.print("Rc.seg");
		//SCMn.print("Rn.seg");
		//SCMf.print("Rf.seg");
		R.print("R.tab");
		Rc.print("Rc.tab");
		Rn.print("Rn.tab");
		Rf.print("Rf.tab");
	}
	
	public static void main(String [] args) throws Exception{Options options = new Options();
	SimulationDataGeneratorDecomposed  S = new SimulationDataGeneratorDecomposed();
		options.addOption("r", "conrario", true, "concordant sample ratio (" + S.concordantSampleRatio + ")");
		options.addOption("c", "conpos", true, "concordant  position (" + MyFunc.join(",",S.concordantProbe)+ ")");
		options.addOption("R", "conrario2", true, "concordant sample ratio for N (" + S.concordantSampleRatio2 + ")");
		options.addOption("C", "conpos2", true, "concordant  position for N(" + MyFunc.join(",",S.concordantProbe2)+ ")");
		options.addOption("l", "lcon", true, "concordant segment length(" + S.concordantSegmentLength + ")");
		options.addOption("L", "lcoc2", true, "concordant segment length for N(" + S.concordantSegmentLength2 + ")");
		options.addOption("s", "lnonc", true, "nonconcordant segment length average (" + S.nonconcordantSegmentLengthMean + ")");
		options.addOption("n", "nnoc", true, "nonconcordant segment number avarage(" + S.nonconcordantSegmentNumberMean + ")");
		options.addOption("N", "ncon2", true, "concordant segment number for N (" + S.concordantSegmentNumber2 + ")");
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
		if(commandLine.hasOption("l")){
			S.setConcordantSegmentLength(Integer.valueOf(commandLine.getOptionValue("l")));
		}
		if(commandLine.hasOption("L")){
			S.setConcordantSegmentLength2(Integer.valueOf(commandLine.getOptionValue("L")));
		}
		if(commandLine.hasOption("r")){
			S.setConcordantSampleRatio(Double.valueOf(commandLine.getOptionValue("r")));
		}
		if(commandLine.hasOption("R")){
			S.setConcordantSampleRatio2(Double.valueOf(commandLine.getOptionValue("R")));
		}
		if(commandLine.hasOption("s")){
			S.setNonconcordantSegmentLengthMean(Integer.valueOf(commandLine.getOptionValue("s")));
		}
		if(commandLine.hasOption("n")){
			S.setNonconcordantSegmentNumberMean(Integer.valueOf(commandLine.getOptionValue("n")));
		}
		if(commandLine.hasOption("N")){
			S.setConcordantSegmentNumber2(Integer.valueOf(commandLine.getOptionValue("N")));
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
		if(commandLine.hasOption("C")){
			S.setConcordantProbe2(commandLine.getOptionValue("C"));
		}
		
		S.generateData();
		S.print();
	}
	
}
