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

public class FOCNAP2 {
	List <Integer> D;
	FOCNAP F;
	MyMat Score;
	MyMat CNA;
	
	double zcutoff = 5;
	double ampcutoff = 0.1;
	double delcutoff = 0.1;
	
	public FOCNAP2(String segFile, String geneBedFile) throws IOException, DataFormatException {
		F = new FOCNAP(segFile, geneBedFile);
		D = new ArrayList<Integer>();
		D.add(100000);
		D.add(500000);
		D.add(1000000);
	}
	
	
	@SuppressWarnings("null")
	void detectCNA(){
		delcutoff = F.CN.percentile(delcutoff);
		ampcutoff = F.CN.percentile(1-ampcutoff);
		//System.err.println(ampcutoff + " " + delcutoff);
		CNA = new MyMat(Score.getRowNames(),Score.getColNames());
		for(String g: CNA.getRowNames()){
			int c = F.GI.chr(g); 
			int s = F.GI.start(g); 
			int e = F.GI.end(g);
			int M = (s + e)/2;
			for(String samp: CNA.getColNames()){
				SegmentContainer cn = F.CN.get(samp);
				if(!cn.contains(c)){
					continue;
				}
				Segment seg = cn.get(c,M);
				if(seg != null){
					if(seg.value() >= ampcutoff){
						CNA.set(g, samp, 1.0);
					}else if(seg.value() <= delcutoff){
						CNA.set(g, samp, -1.0);
					}
				}
			}
		}
	}
	
	void perform(){
		List <MyMat> Scores = new  ArrayList<MyMat>();
		for(Integer d: D){
			F.d = d;
			F.perform();
			Scores.add(new MyMat(F.Score));
			F.resetScore();
		}
		Score = Scores.get(0);
		for(int k=1; k<D.size();k++){
			MyMat tmp = Scores.get(k);
			for(int i=0; i<Score.rowSize();i++){
				for(int j=0; j<Score.colSize();j++){
					if(Math.abs(Score.get(i, j)) < Math.abs(tmp.get(i, j))){
						Score.set(i,j, tmp.get(i, j));
					}
				}
			}
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
	
	void discretizeScores(){
		calculateZscores();
		detectCNA();
		for(int i = 0; i < Score.rowSize(); i++){
			for(int j = 0; j < Score.colSize(); j++){
				if(Score.get(i, j) >= zcutoff & CNA.get(i, j) >= 0.5){
				//if(Score.get(i, j) >= zcutoff ){
					Score.set(i, j, 1);
				}else if(Score.get(i, j) <= -zcutoff & CNA.get(i, j) <= -0.5){
				//}else if(Score.get(i, j) <= -zcutoff){
					Score.set(i, j, -1);
				}else{
					Score.set(i, j, 0);
				}
			}
		}
	}
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("k", "kernel", true, "kernel width (defoult: 100000;500000;1000000)");
		options.addOption("a", "ampcutoff", true, "amplification cutoff (0.1)");
		options.addOption("d", "delcutoff", true, "deletion cutoff (0.1)");
		options.addOption("z", "zcutoff", true, "zscore cutoff (3)");
		options.addOption("o", "outfile", true, "output file");
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
		FOCNAP2 F2 = new FOCNAP2(argList.get(0),argList.get(1));
		if(commandLine.hasOption("k")){
			List <String> tmp =  MyFunc.split(";", commandLine.getOptionValue("k"));
			F2.D.clear();
			for(String S:tmp){
				F2.D.add(Integer.valueOf(S));
			}
		}
		if(commandLine.hasOption("a")){
			F2.ampcutoff = Double.valueOf(commandLine.getOptionValue("k"));
		}
		if(commandLine.hasOption("d")){
			F2.delcutoff = Double.valueOf(commandLine.getOptionValue("d"));
		}
		if(commandLine.hasOption("z")){
			F2.zcutoff = Double.valueOf(commandLine.getOptionValue("z"));
		}
		F2.perform();
		//F2.calculateZscores();
		F2.discretizeScores();
		if(commandLine.hasOption("o")){
			F2.Score.print(commandLine.getOptionValue("o"));
		}else{
			F2.Score.print();
		}
	}

}
