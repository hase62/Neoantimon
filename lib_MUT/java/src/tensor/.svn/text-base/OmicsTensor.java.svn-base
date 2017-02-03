package tensor;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.DataFormatException;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

import utility.*;

public class OmicsTensor extends Order3Tensor {
	private static final long serialVersionUID = 1L;
	//order1: gene
	//order2: sample
	//order3: omics
	
	public OmicsTensor(int n1, int n2, int n3) {
		super(n1, n2, n3);
	}

	public OmicsTensor(List<String> name1,List<String> name3, List<String> name2){
		super(name1, name2, name3);
	}
	
	public OmicsTensor(List<MyMat> Mlist){
		super(Mlist);
	}
	
	
	public OmicsTensor(String infile) throws IOException, DataFormatException{
		super(infile);
	}
	
	public static OmicsTensor readOmicsTensorFromTxet(List<String> infiles) throws IOException, DataFormatException{
		List <MyMat> Mlist = new ArrayList<MyMat>();
		List <String> name3 = new ArrayList<String>();
		Pattern p = Pattern.compile("^.*/(.*)");
		Pattern p2 = Pattern.compile("([a-zA-Z0-9]+?)[.].*");
		for(int i = 0;i < infiles.size();i++){
			Mlist.add(MyMat.readMyMatFromText(infiles.get(i)));
			 Matcher m = p.matcher(infiles.get(i));
			 if(m.matches()){
				 m = p2.matcher(m.group(1));
			 }else{
				 m = p2.matcher(infiles.get(i));
			 }
			 if(m.matches()){
				 name3.add(m.group(1));
			 }else{
				 name3.add(infiles.get(i));
			 }
		}
		OmicsTensor T = new OmicsTensor(Mlist);
		T.setName3(name3);
		return T;
	}
	
	public void filterGenes(double d){
		MyMat VarMat = new MyMat(name1, name3);
		int i,j,k;
		for(k=0; k<n3; k++){
			for(i=0; i<n1; i++){
				List <Double> tmp = new ArrayList<Double>();
				for(j=0; j<n2; j++){
					tmp.add(M[i][j][k]);
				}
				double s = MyFunc.sd(tmp);
				VarMat.set(i,k,s);
			}
		}
		VarMat.normalizeCols();
		SortedMap < Double, List<Integer>> sm = new TreeMap<Double, List<Integer>>();
		for(i=0; i<n1; i++){
			double s = MyFunc.sum(VarMat.getRow(i));
			if(sm.containsKey(s)){
				 (sm.get(s)).add(i);
			 }else{
				 List<Integer> tmp = new ArrayList<Integer>();
				 tmp.add(i);
				 sm.put(s, tmp);
			 }
		}
		
		 List< Double > tmp = new ArrayList<Double>(sm.keySet());
		  Collections.reverse(tmp);
		 List<Integer> tmp2 = new ArrayList<Integer>(); 
		  if(d < 1){
		    d = Math.round(n3*d);
		  } 
		 for(i=0;;i++){
			 List<Integer> tmp3 = sm.get(tmp.get(i));
			 if(tmp3.size()+tmp2.size() <= d){
			  tmp2.addAll(sm.get(tmp.get(i)));
			 }else{
				break;
			 }
		 }
		 reorderOrder1(tmp2);
	}
	
	public void extractGenes(List <String> S){
		reorderOrder1(S);
	}	
	
	public void normalizeAcrossSamples(){
		int i,j,k;
		for(k=0; k<n3; k++){
			for(i=0; i<n1; i++){
				List <Double> tmp = new ArrayList<Double>();
				for(j=0; j<n2; j++){
					tmp.add(M[i][j][k]);
				}
				tmp = MyFunc.normalize(tmp);
				for(j=0; j<n2; j++){
					M[i][j][k] = tmp.get(j);
				}
			}
		}
	}
	public void normalizeAcrossGenes(){
		int i,j,k;
		for(k=0; k<n3; k++){
			for(j=0; j<n2; j++){
				List <Double> tmp = new ArrayList<Double>();
				for(i=0; i<n1; i++){
					tmp.add(M[i][j][k]);
				}
				tmp = MyFunc.normalize(tmp);
				for(i=0; i<n1; i++){
					M[i][j][k] = tmp.get(i);
				}
			}
		}
	}

	public void normalizeWithinOmicsSlice(){
		int i,j,k;
		for(k=0; k<n3; k++){
			MyMat slice = getOrder3Slice(k);
			double m = MyFunc.mean(slice.asList());
			double sd = MyFunc.sd(slice.asList());
			for(i=0; i<n1; i++){
				for(j=0; j<n2; j++){
					M[i][j][k] = (M[i][j][k] - m)/sd;
				}
			}
		}
	}
	
	public void printDimension(){
		System.out.println("gene:\t" + n1);
		System.out.println("sample:\t" + n2);
		System.out.println("omics:\t" + n3);
	}
	public void printStatistics(){
		List <Double> l  = asList();
 		System.out.println("mean:\t" + MyFunc.mean(l));
		System.out.println("sd:\t" + MyFunc.sd(l));
		System.out.println("max:\t" + MyFunc.max(l));
		System.out.println("min:\t" + MyFunc.min(l));
	}
	
	public void print(){
		System.out.print(this);
	}
	
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("v", "varfilt", true, "variation filter for genes");
		options.addOption("s", "normsamp", false,  "normalize across samples");
		options.addOption("g", "normgen", false,  "normalize across genes");
		options.addOption("o", "normomics", false,  "normalize within omics slice");
		options.addOption("p", "prop", false,  "show property");
		options.addOption("m", "multifile", false, "combine multiple matrix files");
		options.addOption("G", "exgen", true, "extract genes");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp("Expression [options] tensorfile", options);
			return;
		}
		List <String> argList = commandLine.getArgList();
		if(!commandLine.hasOption("m")  && argList.size() != 1){
			formatter.printHelp("Expression [options] tensorfile ", options);
			return;
		}
		OmicsTensor T;
		if(commandLine.hasOption("m")){
			T = readOmicsTensorFromTxet(argList);
		}else{
			T = new  OmicsTensor(argList.get(0));
		}	
		if(commandLine.hasOption("p")){
			T.printDimension();
			T.printStatistics();
			return;
		}	
		if(commandLine.hasOption("v")){
			T.filterGenes(Double.valueOf(commandLine.getOptionValue("v")));
		}
		if(commandLine.hasOption("s")){
			T.normalizeAcrossSamples();
		}
		if(commandLine.hasOption("g")){
			T.normalizeAcrossGenes();
		}
		if(commandLine.hasOption("G")){
			T.extractGenes(MyFunc.readStringList(commandLine.getOptionValue("G")));
		}
		if(commandLine.hasOption("o")){
			T.normalizeWithinOmicsSlice();
		}
		T.print();
	}
	
		
		
		
	
	

}
	
	
	
	
	

