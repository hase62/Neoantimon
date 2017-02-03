package snp;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.util.zip.DataFormatException;

import utility.MyFunc;

public class GeneInfo {

	Map <String, Integer[]> gene;
	
	public void sort(){
		List <String> tmp = new ArrayList<String>(gene.keySet());
		java.util.Collections.sort(tmp,new GeneComp());
		Map <String, Integer[]> probe2 = new LinkedHashMap<String, Integer[]>();
		for(String s: tmp){
			probe2.put(s, gene.get(s));
		}
		gene = probe2;
	}
	
	public class GeneComp implements Comparator<String> {
        public int compare(String g1, String g2) {
        		if(chr(g1) == chr(g2)){
        			if(start(g1) < start(g2)){
        				return -1;
        			}
        			if(start(g1) > start(g2)){
        				return 1;
        			}
        				
        		}else{
        			if(chr(g1) < chr(g2)){
        				return -1;
        			}
        			if(chr(g1) > chr(g2)){
        				return 1;
        			}
        		}
				return 0;
        }
	}
	
	public GeneInfo(){
		gene = new LinkedHashMap<String, Integer[]>();
	}
	public void clear(){
		gene = new LinkedHashMap<String, Integer[]>();
	}
	
	public GeneInfo(GeneInfo GI){
		gene = new LinkedHashMap<String, Integer[]>(GI.gene);
	}
	
	public void put(String id, int chr, int start, int end, int strand){
		Integer[] tmp = {chr, start, end, strand};
		gene.put(id, tmp);
	}
	
	public int chr(String id){
		return gene.get(id)[0];
	}
	
	public int start(String id){
		return gene.get(id)[1];
	}
	
	public int end(String id){
		return gene.get(id)[2];
	}
	

	public int strand(String id){
		return gene.get(id)[3];
	}
	
	public Set<String> getGeneSet(){
		return gene.keySet();
		
	}
	public static GeneInfo getGeneInfoFromTsvFile(String inFile) throws IOException, DataFormatException{
		BufferedReader inputStream = new BufferedReader(new FileReader(inFile));
		String line = inputStream.readLine();
		GeneInfo GI = new GeneInfo();
		while((line = inputStream.readLine()) != null){
			List<String> tmp = Arrays.asList((line.split("\t")));
			String id = tmp.get(0);
			int chr = Integer.valueOf(tmp.get(1));
			int start = Integer.valueOf(tmp.get(2));
			int end = Integer.valueOf(tmp.get(3));
			Integer strand = 0;
			if(tmp.get(4).equals("1")|tmp.get(4).equals("+")){
				strand = 1;
			}else if (tmp.get(4).equals("-1")|tmp.get(4).equals("-")){
				strand = -1;
			}
			GI.put(id, chr, start, end, strand);
		}
		GI.sort();
		return GI;
	}

	public static int chrInt(String s){
		String regex = "chr";
		Pattern p = Pattern.compile(regex);
		Matcher m = p.matcher(s);
		String s2 = m.replaceFirst("");
		if(s2.equals("X")){
			s2 = "23";
		}
		if(s2.equals("Y")){
			s2 = "24";
		}
		return Integer.valueOf(s2);
	}
	
	
	public static GeneInfo getGeneInfoFromBedFile(String inFile) throws IOException, DataFormatException{
		BufferedReader inputStream = new BufferedReader(new FileReader(inFile));
		//String line = inputStream.readLine();
		String line = null;
		GeneInfo GI = new GeneInfo();
		while((line = inputStream.readLine()) != null){
			List<String> tmp = Arrays.asList((line.split("\t")));
			if(tmp.size()<3){
				continue;	
			}
			//chrom start end  name score strand 
			int chr = chrInt(tmp.get(0));
			int start = Integer.valueOf(tmp.get(1));
			int end = Integer.valueOf(tmp.get(2));
			String id;
			if(tmp.size()==3){
				id = MyFunc.join("_", tmp);
			}else{
				id = tmp.get(3);
			}
			Integer strand = 0;
			if(tmp.size()>4){
				if(tmp.get(4).equals("1")|tmp.get(4).equals("+")){
					strand = 1;
				}else if (tmp.get(4).equals("-1")|tmp.get(4).equals("-")){
					strand = -1;
				}
			}
			GI.put(id, chr, start, end, strand);
		}
		GI.sort();
		return GI;
	}
	
	
	public void filter(SegmentContainer S){
		GeneInfo tmp = new GeneInfo(this);
		clear();
		for(String s: tmp.getGeneSet()){
			int chr = tmp.chr(s);
			int start  = tmp.start(s);
			int end = tmp.end(s);
			int strand = tmp.strand(s);
			if(S.contains(chr,start) | S.contains(chr,end)){
				put(s,chr,start, end, strand);
			}
		}
	}
	
	public void filter(int chr){
		GeneInfo tmp = new GeneInfo(this);
		clear();
		for(String s: tmp.getGeneSet()){
			int c = tmp.chr(s);
			int start  = tmp.start(s);
			int end = tmp.end(s);
			int strand = tmp.strand(s);
			if(chr == c){
				put(s,chr,start, end, strand);
			}
		}
	}
	

	public void filter(int chr, int start, int end){
		GeneInfo tmp = new GeneInfo(this);
		clear();
		for(String s: tmp.getGeneSet()){
			int c = tmp.chr(s);
			int st  = tmp.start(s);
			int en = tmp.end(s);
			int str = tmp.strand(s);
			if(chr == c & st >= start & en <= end){
				put(s,c,st,en,str);
			}
		}
	}
		
	public String  toString(){
		List <String> out = new ArrayList<String>();
		out.add("Gene" + "\t" +  "Chrom" + "\t" +  "Start" + "\t" +  "End" + "\t" +  "Strand");
		for(String p: getGeneSet()){
			out.add(p + "\t" +  chr(p) + "\t" + start(p) + "\t" + end(p) + "\t" + strand(p));
		}
		return MyFunc.join("\n", out);
	}
	
	public ProbeInfo toProbeInfo(){
		ProbeInfo tmp = new ProbeInfo();
		for(String s: getGeneSet()){
			int c = chr(s);
			int st  = start(s);
			int en = end(s);
			tmp.put(s,c,(st+en)/2);
		}
		tmp.sort();
		return tmp;
	}
	
	

}
