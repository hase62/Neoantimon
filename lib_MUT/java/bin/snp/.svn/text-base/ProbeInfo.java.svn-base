package snp;

import java.io.*;
import java.util.*;
import java.util.zip.DataFormatException;

import utility.MyFunc;
import utility.MyMat;

public class ProbeInfo {

	Map <String, Integer[]> probe;
	
	public void sort(){
		List <String> tmp = new ArrayList<String>(probe.keySet());
		java.util.Collections.sort(tmp,new ProbeComp());
		Map <String, Integer[]> probe2 = new LinkedHashMap<String, Integer[]>();
		for(String s: tmp){
			probe2.put(s, probe.get(s));
		}
		probe = probe2;
	}
	
	public List <Integer> getChrList(){
		Set<Integer> chr = new LinkedHashSet<Integer>();
		for(String p: getProbeSet()){
			if(!chr.contains(chr(p))){
				chr.add(chr(p));
			}
		}
		return new ArrayList<Integer>(chr);
	}
	
	public List <String> getProbeSetBychr(int i){
		List<String> tmp = new ArrayList<String>();
		for(String p: getProbeSet()){
			if(chr(p) == i){
				tmp.add(p);
			}
		}
		return tmp;
	}
	
	
	public class ProbeComp implements Comparator<String> {
        public int compare(String p1, String p2) {
        		if(chr(p1) == chr(p2)){
        			if(pos(p1) < pos(p2)){
        				return -1;
        			}
        			if(pos(p1) > pos(p2)){
        				return 1;
        			}
        				
        		}else{
        			if(chr(p1) < chr(p2)){
        				return -1;
        			}
        			if(chr(p1) > chr(p2)){
        				return 1;
        			}
        		}
				return 0;
        }
	}
	
	
	
	
	public ProbeInfo(){
		probe = new LinkedHashMap<String, Integer[]>();
	}
	
	public void clear(){
		probe = new LinkedHashMap<String, Integer[]>();
	}
	
	public ProbeInfo(ProbeInfo PI){
		probe = new LinkedHashMap<String, Integer[]>(PI.probe);
	}
	
	public ProbeInfo(Map<String, Integer[]> probe){
		this.probe = probe;
	}
	
	public void put(String id, int chr, int pos){
		Integer[] tmp = {chr, pos};
		probe.put(id, tmp);
	}
	

	
	public String getPreviousNearestProbe(int chr, int pos){
		int min = Integer.MAX_VALUE;
		String nearest = "";
		for(String s: getProbeSet()){
			if(chr(s) == chr){
				int tmp = Math.abs(pos(s) -pos);
				if(pos(s) <= pos & tmp < min){
					min = tmp;
					nearest = s;
				}
			}
		}
		return nearest;
	}
	
	public String getFollowingNearestProbe(int chr, int pos){
		int min = Integer.MAX_VALUE;
		String nearest = "";
		for(String s: getProbeSet()){
			if(chr(s) == chr){
				int tmp = Math.abs(pos(s) -pos);
				if(pos(s) >= pos & tmp < min){
					min = tmp;
					nearest = s;
				}
			}
		}
		return nearest;
	}
	
	public int getProbeCount(int chr, int start, int end){
		int count = 0;
		for(String s: getProbeSet()){
			if(chr(s) == chr){
				if(pos(s) >= start & pos(s) <= end){
					count++;
				}
			}
		}
		return count;
	}
	
	
	public int chr(String id){
		return probe.get(id)[0];
	}
	
	public int pos(String id){
		return probe.get(id)[1];
	}
	
	public Set<String> getProbeSet(){
		return probe.keySet();
		
	}
	public static ProbeInfo getProbeInfoFromTsvFile(String inFile) throws IOException, DataFormatException{
		BufferedReader inputStream = new BufferedReader(new FileReader(inFile));
		String line = inputStream.readLine();
		ProbeInfo PI = new ProbeInfo();
		while((line = inputStream.readLine()) != null){
			List<String> tmp = Arrays.asList((line.split("\t")));
			String id = tmp.get(0);
			int chr = Integer.valueOf(tmp.get(1));
			int pos = Integer.valueOf(tmp.get(2));
			PI.put(id, chr, pos);
		}
		PI.sort();
		return PI;
	}
	
	public void thinDown(int n){
		Map <String, Integer[]> probe2 = new LinkedHashMap<String, Integer[]>();
		int i = 0;
		for(String s:probe.keySet()){
			i++;
			if(i == n){
				Integer[] tmp = {chr(s), pos(s)};
				probe2.put(s, tmp);
				i=0;
			}
		}
		probe = probe2;
	}
	
	public void filter(SegmentContainer S){
		ProbeInfo tmp = new ProbeInfo(this);
		clear();
		for(String s: tmp.getProbeSet()){
			int chr = tmp.chr(s);
			int pos = tmp.pos(s);
			if(S.contains(chr,pos)){
				put(s,chr,pos);
			}
		}
	}
	
	public void filter(int chr){
		ProbeInfo tmp = new ProbeInfo(this);
		clear();
		for(String s: tmp.getProbeSet()){
			int c = tmp.chr(s);
			int p = tmp.pos(s);
			if(chr == c){
				put(s,c,p);
			}
		}
	}
	

	public void filter(int chr, int start, int end){
		ProbeInfo tmp = new ProbeInfo(this);
		clear();
		for(String s: tmp.getProbeSet()){
			int c = tmp.chr(s);
			int p = tmp.pos(s);
			if(chr == c & p >= start & p <= end){
				put(s,c,p);
			}
		}
	}
	
	public void filter(Set<String> filter){
		ProbeInfo tmp = new ProbeInfo(this);
		clear();
		for(String s: tmp.getProbeSet()){
			int c = tmp.chr(s);
			int p = tmp.pos(s);
			if(filter.contains(s)){
				put(s,c,p);
			}
		}
	}
	
	
	public String  toString(){
		List <String> out = new ArrayList<String>();
		out.add("Probe" + "\t" +  "Chrom" + "\t" +  "BasePair");
		for(String p: getProbeSet()){
			out.add(p + "\t" +  chr(p) + "\t" + pos(p));
		}
		return MyFunc.join("\n", out);
	}
	
	public static ProbeInfo generatePsuedoProbeInfoByIntervalLength(int interval){
		ProbeInfo PI = new ProbeInfo();
		for(int i = 1; i <= 24; i++){
			int M = SegmentContainer.chrLength(i);
			int m = 1;
			for(int p = m; p <= M; p+= interval){
				PI.put("chr" + i + "_" + p, i, p);
			}
		}
		
		return PI;
	}
	
	public static ProbeInfo generatePsuedoProbeInfo(int probeCount){
		return generatePsuedoProbeInfoByIntervalLength((int)(3*Math.pow(10, 9)/probeCount));
	}

	
	
}
