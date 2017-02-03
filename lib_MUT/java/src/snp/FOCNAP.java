package snp;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.DataFormatException;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;


import sun.reflect.Reflection;
import utility.MyFunc;
import utility.MyMat;

public class FOCNAP {

	SegmentContainerMap CN;
	GeneInfo GI;
	List <String> sample;
	List <String> gene;
	int d =  1000000;  //1M
	int n = 100;
	MyMat Score;
	double median;
	
	public FOCNAP(String segFile, String geneBedFile) throws IOException, DataFormatException{
		CN = new SegmentContainerMap(segFile);
		GI =GeneInfo.getGeneInfoFromBedFile(geneBedFile);
		sample = new ArrayList<String>(new ArrayList<String>(CN.keySet()));
		gene = new ArrayList<String>(new ArrayList<String>(GI.getGeneSet()));
		Score = new MyMat(gene, sample);
		median = CN.percentile(0.5);
	}
	void resetScore(){
		Score = new MyMat(gene, sample);
	}
	
	void perform(){
		int j = 1;
		List<String> gene2 = new ArrayList<String>();
		for(String g: gene){
			System.err.println(g + " " + j + "/" + gene.size());
			j++;
			int c = GI.chr(g); 
			int s = GI.start(g); 
			int e = GI.end(g);
			int M = (s + e)/2;
			int S = M-d*5;
			int E = M+d*5;	
			int tmp = E - S;
			boolean hasValue = false; 
			tmp /= n;
			List <Integer> p = new ArrayList<Integer>();
			List <Double> P = new ArrayList<Double>();
			for(int i = 0; i < n; i++){
				int k = S+i*tmp;
				p.add(k);
				P.add(kernel(k,M));
			}
			double m = MyFunc.mean(P); 
			for(String samp: sample){
				SegmentContainer cn = CN.get(samp);
				if(!cn.contains(c)){
					continue;
				}
				double score = 0;
				List <Double> v = new ArrayList<Double>();
				for(int i = 0; i < n; i++){
					Segment seg = cn.get(c,p.get(i));
					if(seg != null){
						v.add(seg.value());
					}else{
						v.add(median);
					}
				}
				double m2 = MyFunc.mean(v);
				for(int i = 0; i < v.size(); i++){
					score += (v.get(i)-m2)*(P.get(i)-m) ;
					hasValue = true;
				}
				score /= v.size();
				Score.set(g, samp,score);
			}
			if(hasValue){
				gene2.add(g);
			}
		}
		if(gene.size() != gene2.size()){
			Score = Score.getSubMatByRow(gene2);
			gene=gene2;
		}
	}
	
	private double  kernel(int x, int M){
		if(Math.abs(M-x)<=d/2){
			return 1;
		}else{
			return 0;
		}
	}
	
	
	void calculateZscores(){
		List <Double> tmp = Score.asList();
		double m = MyFunc.mean(tmp);
		double sd = MyFunc.sd(tmp);
		for(int i = 0; i < Score.rowSize(); i++){
			for(int j = 0; j < Score.colSize(); j++){
				Score.set(i, j, (Score.get(i, j)-m)/sd);
			}
		}
	}
	
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("k", "kernel", true, "kernel wideth (defoult: 1M)");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] segFile geneBedFile", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 2){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] segFile geneBedFile", options);
			return;
		}	
		FOCNAP F = new FOCNAP(argList.get(0),argList.get(1));
		if(commandLine.hasOption("k")){
			F.d = Integer.valueOf(commandLine.getOptionValue("k"));
		}
		F.perform();
		F.calculateZscores();
		F.Score.print();
	}
	
	
}
