package snp;


import java.io.*;
import java.util.*;
import java.util.zip.DataFormatException;

import org.apache.commons.cli.*;

import utility.MyFunc;
import utility.MyMat;

import sun.reflect.Reflection;

public class SegmentContainerMap {

	Map<String,SegmentContainer> SCM;
	
	public SegmentContainerMap(){
		SCM = new LinkedHashMap <String,SegmentContainer> ();
	}
	public SegmentContainerMap(Map<String,SegmentContainer> SCM){
		this.SCM = new LinkedHashMap <String,SegmentContainer> (SCM);
	}
	public SegmentContainerMap(SegmentContainerMap SCM){
		this.SCM = new LinkedHashMap <String,SegmentContainer> (SCM.SCM);
	}
	public SegmentContainerMap(String segFile) throws IOException, DataFormatException{
		BufferedReader inputStream = new BufferedReader(new FileReader(segFile));
		String line = inputStream.readLine();
		String id = "";
		SegmentContainer SC = new SegmentContainer();
		SCM = new LinkedHashMap <String,SegmentContainer> ();
		while((line = inputStream.readLine()) != null){
			List<String> tmp = Arrays.asList((line.split("\t")));
			if(tmp.size() != 6){
				continue;
				//throw new DataFormatException("");
			}
			if(id == ""){
				id = tmp.get(0);
			}else if(!id.equals(tmp.get(0))){
				if(SCM.containsKey(id)){
					SCM.get(id).bind(SC);
				}else{
					SCM.put(id, new SegmentContainer(SC));
				}
				SC = new SegmentContainer();
				id = tmp.get(0);
			}
			int chr;
			if(tmp.get(1).equals("X")){
				chr = 23;
			}else if(tmp.get(1).equals("Y")){
				chr = 24;
			}else{
				chr = (int)Double.parseDouble(tmp.get(1));
			}
			int start = (int)Double.parseDouble(tmp.get(2));
			int end = (int)Double.parseDouble(tmp.get(3));
			if(!(start<end)){
				continue;
			}
			int probeCount = (int)Double.parseDouble(tmp.get(4));
			double value = Double.parseDouble(tmp.get(5));
			SC.add(new Segment(chr, start, end, probeCount, value));
		}
		if(SCM.containsKey(id)){
			SCM.get(id).bind(SC);
		}else{
			SCM.put(id, new SegmentContainer(SC));
		}
	}
	
	public static SegmentContainerMap readFromCNAG(String cscFile) throws IOException{
		BufferedReader inputStream = new BufferedReader(new FileReader(cscFile));
		String line = inputStream.readLine();
		String id = "";
		SegmentContainer SC = new SegmentContainer();
		Map<String,SegmentContainer> SCM = new LinkedHashMap <String,SegmentContainer> ();
		while((line = inputStream.readLine()) != null){
			List<String> tmp = Arrays.asList((line.split("\t")));
			if(tmp.size() < 7){
				continue;
				//throw new DataFormatException("");
			}
			if(id == ""){
				id = tmp.get(0);
			}else if(!id.equals(tmp.get(0))){
				if(SCM.containsKey(id)){
					SCM.get(id).bind(SC);
				}else{
					SCM.put(id, new SegmentContainer(SC));
				}
				SC = new SegmentContainer();
				id = tmp.get(0);
			}
			int chr = (int)Double.parseDouble(tmp.get(2));
			int start = (int)Double.parseDouble(tmp.get(4));
			int end = (int)Double.parseDouble(tmp.get(6));
			if(!(start<end)){
				continue;
			}
			int probeCount = 0;
			double value = Double.parseDouble(tmp.get(1));
			SC.add(new Segment(chr, start, end, probeCount, value));
		}
		if(SCM.containsKey(id)){
			SCM.get(id).bind(SC);
		}else{
			SCM.put(id, new SegmentContainer(SC));
		}
		return new SegmentContainerMap(SCM);
	}
	
	public SegmentContainerMap(MyMat M, ProbeInfo PI){
		SCM = new LinkedHashMap <String,SegmentContainer> ();
		PI.filter(new HashSet<String>(M.getRowNames()));
		PI.sort();
		List <String> probe = new ArrayList<String>(PI.getProbeSet());
		for(String s: M.getColNames()){
			SegmentContainer SC = new SegmentContainer();
			List<Double> v = new ArrayList<Double>();
			for(String p: probe){
				v.add(M.get(p,s));
			}
			int start = PI.pos(probe.get(0));
			int end =  PI.pos(probe.get(0));
			int chr = PI.chr(probe.get(0));
			int probeCount = 1;
			double value 	= v.get(0);
 			for(int i=1; i <  v.size(); i++){
 				String p = probe.get(i);
 				if(chr < PI.chr(p) || value != v.get(i)){
 					SC.add(new Segment(chr, start, end, probeCount, value));
 					chr = PI.chr(probe.get(i));
 					start = PI.pos(probe.get(i));
 					probeCount = 0;
 				}
 				end = PI.pos(probe.get(i));
 				probeCount++;
 				value 	= v.get(i);
 			}
 			SC.add(new Segment(chr, start, end, probeCount, value));
 			SCM.put(s, SC);
		}	
	}
	
		
	public boolean checkSegmentContainerConsistency(){
		for(String s: SCM.keySet()){
			SegmentContainer S =  SCM.get(s);
			if(!S.checkSegmentConsistency()){
				System.err.println("WARN: "  + s + " is inconsistent");
				return false;
			}
		}
		return true;
	}
	
	public void resolveSegmentContainerInconsistency(){
		for(String s: SCM.keySet()){
			SegmentContainer S =  SCM.get(s);
			if(!S.checkSegmentConsistency()){
				System.err.println("WARN: "  + s + " is inconsistent");
				S.cleanSegmentContainer();
			}
		}
	}
	
	public void removeInconsistentSegmentContainer(){
		for(String s: SCM.keySet()){
			SegmentContainer S =  SCM.get(s);
			if(!S.checkSegmentConsistency()){
				System.err.println("WARN: "  + s + " is inconsistent");
				SCM.remove(s);
			}
		}
	}
	
	public void bind(SegmentContainerMap S){
		for(String s:S.keySet()){
			if(SCM.containsKey(s)){
				SCM.get(s).bind(S.get(s));
			}else{
				put(s, S.get(s));
			}	
		}
	}
	
	public SegmentContainerMap getSubMap(List <String> sample){
		SegmentContainerMap S = new SegmentContainerMap(); 
		for(String s: sample){
			if(containsKey(s)){
				S.put(s, get(s));
			}
		}
		return S;
	}
	
	
	public SegmentContainer get(String s){
		return SCM.get(s);
	}
	
	public void put(String s, SegmentContainer SC){
		SCM.put(s, SC);
	}
	
	public Set<String> keySet(){
		return SCM.keySet();
	}
	
	public Collection<SegmentContainer> values(){
		return SCM.values();
	}
	
	public boolean containsKey(String s){
		return SCM.containsKey(s);
	}
	
	public void clear(){
		SCM.clear();
	}
	
	public void removeXY(){
		for(SegmentContainer s: values()){
			s.removeXY();
		}
	}
	public void removeY(){
		for(SegmentContainer s: values()){
			s.removeY();
		}
	}
	
	public void mergeSegments(){
		for(SegmentContainer s: values()){
			s.mergeSegments();
		}
	}
	
	
	public void fillGaps(){
		for(SegmentContainer s: values()){
			s.fillGaps();
		}
	}
	
	public void removeShortSegmentsByProbeCount(int probeCount){
		for(SegmentContainer s: values()){
			s.removeShortSegmentsByProbeCount(probeCount);
		}	
	}
	
	public void removeShortSegmentsByLength(int length){
		for(SegmentContainer s: values()){
			s.removeShortSegmentsByLength(length);
		}	
	}
	
	public void fillGaps(double defaultValue){
		for(SegmentContainer s: values()){
			s.fillGaps(defaultValue);
		}
	}
	
	public boolean checkConsistency(){	
		Map <Integer, Integer> chrCount = new HashMap<Integer, Integer>();
		Map <Integer, Integer> chrMax = new HashMap<Integer, Integer>();
		Map <Integer, Integer> chrMin = new HashMap<Integer, Integer>();
		for(String s: SCM.keySet()){
			SegmentContainer sc = SCM.get(s);
			for(Integer c: sc.chrList()){
				if(chrCount.containsKey(c)){
					chrCount.put(c, chrCount.get(c)+1);
				}else{
					chrCount.put(c,1);
				}
				int m = sc.minPosAtChr(c);
				int M = sc.maxPosAtChr(c);
				if(chrMax.containsKey(c)){
					if(chrMax.get(c) !=  M){
						System.err.println("ERR: inconsistent max pos");
						System.err.println(c + " " + M + " " + chrMax.get(c));
						return false;
					}
				}else{
					chrMax.put(c, M); 	
				}
				if(chrMin.containsKey(c)){
					if(chrMin.get(c) !=  m){
						System.err.println("ERR: inconsistent min pos");
						return false;
					}
				}else{
					chrMin.put(c, m); 	
				}
			}
		}	
		for(Integer i: chrCount.keySet()){
			if(SCM.size()!= chrCount.get(i)){
				System.err.println("ERR: inconsistent chr count " + i);
				return false;	
			}
		}
		return true;
	}
	
	
	public void resolveInconsistency(){	
		Map <Integer, Integer> chrCount = new HashMap<Integer, Integer>();
		Map <Integer, Integer> chrMax = new HashMap<Integer, Integer>();
		Map <Integer, Integer> chrMin = new HashMap<Integer, Integer>();
		
		for(String s: SCM.keySet()){
			SegmentContainer sc = SCM.get(s);
			for(Integer c: sc.chrList()){
				if(chrCount.containsKey(c)){
					chrCount.put(c, chrCount.get(c)+1);
				}else{
					chrCount.put(c,1);
				}
				int m = sc.minPosAtChr(c);
				int M = sc.maxPosAtChr(c);
				if(chrMax.containsKey(c)){
					if(chrMax.get(c) <  M){
						chrMax.put(c, M);
					}
				}else{
					chrMax.put(c, M); 	
				}
				if(chrMin.containsKey(c)){
					if(chrMin.get(c) >  m){
						chrMin.put(c, m); 
					}
				}else{
					chrMin.put(c, m); 	
				}
			}
		}	
		
		List <Integer> missingChr = new ArrayList<Integer>();
		for(Integer i: chrCount.keySet()){
			if(SCM.size()!= chrCount.get(i)){
				missingChr.add(i);
			}
		}
		
		Map <String, SegmentContainer>  SCM2  = new LinkedHashMap <String,SegmentContainer> (SCM);
			
		for(String s: SCM.keySet()){
			SegmentContainer sc = SCM.get(s);
			for(Integer c: sc.chrList()){
				 if(sc.maxPosAtChr(c) != chrMax.get(c)){
					 if(SCM2.containsKey(s)){
						 System.err.println("Warn: remove " + s);
						 SCM2.remove(s);	
					 }
				 }
				 if(sc.minPosAtChr(c) != chrMin.get(c)){
					 if(SCM2.containsKey(s)){
						 System.err.println("Warn: remove " + s);
						 SCM2.remove(s);	
					 }
				 }	
			}
			for(Integer c: missingChr){
				if(!sc.chrList().contains(c)){
					if(SCM2.containsKey(s)){
						System.err.println("Warn: remove " + s);
						SCM2.remove(s);
					}
				}
			}
		}
		SCM = SCM2;
	}
	public void filter(int chr, int start, int end){
		for(SegmentContainer SC: SCM.values()){
			SC.filter(chr,start,end);
		}
	}
	
	public void filter(int chr){
		for(SegmentContainer SC: SCM.values()){
			SC.filter(chr);
		}
	}

	public Set<Integer> getCommonChrSet(){
		Map<Integer, Integer> chrs = new HashMap <Integer, Integer>();
		for(String s: keySet()){
			for(Integer i: get(s).chrList()){
				if(chrs.containsKey(i)){
					chrs.put(i, chrs.get(i)+1);
				}else{
					chrs.put(i,1);
				}
			}
		}
		Set<Integer> chrs2 = new HashSet<Integer>();
		for(Integer i: chrs.keySet()){
			if(chrs.get(i)==keySet().size()){
				chrs2.add(i);
			}
		}
		return chrs2;
	}
	
	public void sort(){
		for(SegmentContainer SC: SCM.values()){
			SC.sort();
		}
	}
	
	public void shiftMode(double d){
		for(SegmentContainer SC: SCM.values()){
			SC.shiftMode(d);
		}
	}
	
	public void binarize(double d){
		for(SegmentContainer SC: SCM.values()){
			SC.binarize(d);
		}
	}
	
	public void binarizeLess(double d){
		for(SegmentContainer SC: SCM.values()){
			SC.binarizeLess(d);
		}
	}
	
	public void ternarize(double d, double d2){
		for(SegmentContainer SC: SCM.values()){
			SC.ternarize(d, d2);
		}
	}
	
	public double percentile(double d){
		List <Double> tmp = sampleValues(10000);
		return  MyFunc.percentile(tmp, d);
	}
	
	
	
	public void  print(String outfile) throws IOException{
		PrintWriter os = new PrintWriter(new FileWriter(outfile));
		os.println("ID" + "\t" + "chrom" + "\t" + "loc.start" + "\t" + "loc.end" + "\t" +  "num.mark" + "\t" + "seg.mean");
		for(String s:SCM.keySet()){
			SegmentContainer SC  = SCM.get(s);
			for(Segment s2: SC){
				os.println(s + "\t" + s2);
			}
		}
		os.close();
	}
	
	public void  print(){
		System.out.println("ID" + "\t" + "chrom" + "\t" + "loc.start" + "\t" + "loc.end" + "\t" +  "num.mark" + "\t" + "seg.mean");
		for(String s:SCM.keySet()){
			SegmentContainer SC  = SCM.get(s);
			for(Segment s2: SC){
				System.out.println(s + "\t" + s2);
			}
		}
	}
	
	public MyMat toMyMat(ProbeInfo PI){
		List <String> sample = new ArrayList<String>(new ArrayList<String>(SCM.keySet()));
		List <String> probe = new ArrayList<String>(new ArrayList<String>(PI.getProbeSet()));
		MyMat M = new MyMat(probe, sample);		
		for(String s: sample){
			SegmentContainer SC = SCM.get(s);
			SegmentContainer.SegmentIterator SCitr  = SC.iterator();
			Segment S = SCitr.next();
			for(String p: probe){
				int chr = PI.chr(p);
				int pos = PI.pos(p);
				while(chr > S.chr()){
					S = SCitr.next();
				}
				while(pos > S.end()){
				//while(!(pos >= S.start() & pos <= S.end())){
					S = SCitr.next();
				}
				M.set(p, s, S.value());
			}
		}
		return M;
	}
	
	public MyMat toMyMat(GeneInfo GI){
		List <String> sample = new ArrayList<String>(new ArrayList<String>(SCM.keySet()));
		List <String> gene = new ArrayList<String>(new ArrayList<String>(GI.getGeneSet()));
		MyMat M = new MyMat(gene, sample);		
		for(String s: sample){
			SegmentContainer SC = SCM.get(s);
			for(String g: gene){
				double v = SC.getMeanValue(GI.chr(g), GI.start(g), GI.end(g));
				M.set(g, s, v);
			}
		}
		return M;
	}
	
	
	public double min(){
		double D = Double.MAX_VALUE;
		for(SegmentContainer SC: SCM.values()){
			double d = SC.min();
			if(d < D){
				D = d;
			}
		}
		return D;
	}
	
	public double max(){
		double D = -Double.MAX_VALUE;
		for(SegmentContainer SC: SCM.values()){
			double d = SC.max();
			if(d > D){
				D = d;
			}
		}
		return D;
	}
	
	public long totalLength(){
		long L = 0;
		for(SegmentContainer SC: SCM.values()){
			L += SC.totalLength();
		}
		return L;
	}
	
	
	public double mean(){
		double d = 0;
		long L = totalLength();
		for(SegmentContainer SC: SCM.values()){
			long l = SC.totalLength();
			d += SC.mean()*l/L;
		}
		return d;
	}
	
	public double var(){
		double d = 0;
		long L = totalLength();
		for(SegmentContainer SC: SCM.values()){
			long l = SC.totalLength();
			d += SC.var()*l/L;
		}
		return d;
	}
	
	public double sd(){
		return Math.pow(var(),0.5);
	}
	
	public void printStatistics(){
		System.out.println("sample:\t" + SCM.size());
		System.out.println("mean:\t" + mean());
		System.out.println("sd:\t" +  sd());
		System.out.println("max:\t" + max());
		System.out.println("min:\t" + min());
		return;
	}
	
	
	public ProbeInfo generatePsuedoProbeInfo(int probeCount){
		for(SegmentContainer S: values()){
			return S.generatePsuedoProbeInfo(probeCount);
		}
		return null;
	}
	
	public void extractSamples(List <String> samples){
		Set <String> tmp = new HashSet<String>(SCM.keySet());
		for(String s: tmp){
			if(!samples.contains(s)){
				SCM.remove(s);
			}
		}
	}
	
	public void changePosition2NearestProbe(ProbeInfo PI){
		for(SegmentContainer S: values()){
			S.changePosition2NearestProbe(PI);
		}
	}
	
	public List <Double> sampleValues(int n){
		int m = Math.round(((float)n)/SCM.size());
		List <Double> tmp = new ArrayList<Double>();
		for(SegmentContainer SC: SCM.values()){
			tmp.addAll(SC.sampleValues(m));
		}
		return tmp;
	}
	public Map<Integer, Double> getPercentiles(){
		int n = 1000000;
		List <Double> s = sampleValues(n);
		Map <Integer, Double> tmp = new LinkedHashMap<Integer, Double>();
		for(int i= 0; i <=100; i +=10){
			tmp.put(i,  MyFunc.percentile(s, i));
		}
		return tmp;
	}
	
	public void printPercentiles(){
		Map <Integer, Double> tmp  = getPercentiles();
		for(Integer i: tmp.keySet()){
			System.out.println(i + "\t" + tmp.get(i));
		}
	}
	
		
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("p", "probe", true, "probe coordinate file");
		//options.addOption("g", "gene", true, "gene coordinate file");
		options.addOption("b", "bed", true, "gene coordinate file from bed file");
		options.addOption("t", "thinpi", true, "thin down probe info");
		options.addOption("S", "rmshort", true, "remove short segments");
		options.addOption("T", "rmshortbypc", true, "remove short segments by Probe Count");
		options.addOption("n", "nump", true, "# of psuedo probes");
		options.addOption("l", "locus", true, "locus (e.g., 3:120000-130000)");
		options.addOption("r", "rmxy", false, "remove XY");
		options.addOption("y", "rmy", false, "remove Y");
		options.addOption("f", "fillgap", false, "fill gap");
		options.addOption("F", "Fillgap", true, "fill gap with a value");
		options.addOption("N", "npro", false, "change positions to nearest probe");
		options.addOption("m", "mode", true, "shift mode");
		options.addOption("R", "rmincsc", false, "remove inconsistent segment container");
		options.addOption("P", "prst", false, "print statistics");
		options.addOption("s", "sort", false, "sort");
		options.addOption("c", "check", false, "check consistency");
		options.addOption("C", "checksc", false, "check segment container consistency");
		options.addOption("i", "rsin", false, "resolve  inconsistency");
		options.addOption("I", "rsinsc", false, "resolve  segment container inconsistency");
		options.addOption("A", "frommat", true, "from matrix and probe info");
		options.addOption("a", "fromcnag", false, "from cnag");
		options.addOption("M", "tomat", false, "to matrix");
		options.addOption("O", "opi", true, "output probe info");
		options.addOption("e", "exsmp", true, "extract samples");
		options.addOption("B", "bin", true, "binarize");
		options.addOption("L", "binless", true, "binarize so that the smaller is 1");
		options.addOption("E", "percentile", false, "print percentiles");
		options.addOption("G", "binper", true, "binarize (cutoff by percentile)");
		options.addOption("K", "binlessper", true, "binarize so that the smaller is 1 (cutoff by percentile)");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine = null;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] SegFile", options);
			System.exit(1);
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() < 1){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] SegFile", options);
			System.exit(1);
		}
		SegmentContainerMap SCM = null;
		
		if(commandLine.hasOption("A")){
			ProbeInfo PI = 	ProbeInfo.getProbeInfoFromTsvFile(commandLine.getOptionValue("A"));
			SCM = new SegmentContainerMap(new MyMat(argList.get(0)), PI); 
		}if(commandLine.hasOption("a")){
			SCM = readFromCNAG(argList.get(0));
			if(argList.size() > 1){
				for(int i = 1; i < argList.size(); i++){
					SCM.bind(readFromCNAG(argList.get(i)));
				}
			}
		}else{
			SCM = new SegmentContainerMap(argList.get(0));
			if(argList.size() > 1){
				for(int i = 1; i < argList.size(); i++){
					SCM.bind(new SegmentContainerMap(argList.get(i)));
				}
			}
		}
		
		if(commandLine.hasOption("c")){
			if(!SCM.checkConsistency()){
				System.exit(1);
			}
		}
		if(commandLine.hasOption("C")){
			if(!SCM.checkSegmentContainerConsistency()){
				System.exit(1);
			}
		}
		
		if(commandLine.hasOption("i")){
			SCM.resolveInconsistency();
		}
		if(commandLine.hasOption("I")){
			SCM.resolveSegmentContainerInconsistency();
		}
		
		
		if(commandLine.hasOption("s")){
			SCM.sort();
		}
		
		if(commandLine.hasOption("S")){
			SCM.removeShortSegmentsByLength(Integer.valueOf(commandLine.getOptionValue("S")));
		}
		if(commandLine.hasOption("T")){
			SCM.removeShortSegmentsByProbeCount(Integer.valueOf(commandLine.getOptionValue("T")));
		}
		
		if(commandLine.hasOption("f")){
			SCM.fillGaps();
		}
	
		if(commandLine.hasOption("F")){
			SCM.fillGaps(Double.valueOf(commandLine.getOptionValue("F")));
		}
		
			
		ProbeInfo PI = null;
		GeneInfo GI = null;
		List<String> sample = new ArrayList<String>(SCM.keySet());
		if(commandLine.hasOption("p")){
			PI = 	ProbeInfo.getProbeInfoFromTsvFile(commandLine.getOptionValue("p"));
			if(commandLine.hasOption("t")){
				PI.thinDown(Integer.valueOf(commandLine.getOptionValue("t")));
			}
			PI.filter(SCM.get(sample.get(0)));
		}
		
		//if(commandLine.hasOption("g")){
			//GI = GeneInfo.getGeneInfoFromTsvFile(commandLine.getOptionValue("g"));
			//GI.filter(SCM.get(sample.get(0)));
		//}
		if(commandLine.hasOption("b")){
			GI = GeneInfo.getGeneInfoFromBedFile(commandLine.getOptionValue("b"));
			GI.filter(SCM.get(sample.get(0)));
		}
		
		if(commandLine.hasOption("l")){
			List <String> tmp = MyFunc.split(":", commandLine.getOptionValue("l"));
			int chr = Integer.valueOf(tmp.get(0));
			if(tmp.size() == 1){
				SCM.filter(chr);
			}else{
				List <String> tmp2 = MyFunc.split("-", tmp.get(1));
				if(tmp2.size()==2){
					int start = Integer.valueOf(tmp2.get(0));
					int end = Integer.valueOf(tmp2.get(1));
					SCM.filter(chr,start,end);
				}
			}
		}
		if(commandLine.hasOption("m")){
			SCM.shiftMode(Double.valueOf(commandLine.getOptionValue("m")));
		}

		if(commandLine.hasOption("r")){
			SCM.removeXY();
		}
		if(commandLine.hasOption("y")){
			SCM.removeY();
		}
		
		if(commandLine.hasOption("R")){
			SCM.removeInconsistentSegmentContainer();
		}
		
		if(commandLine.hasOption("n")){
			int n = Integer.valueOf(commandLine.getOptionValue("n"));
			PI = SCM.generatePsuedoProbeInfo(n);
		}
		
		
		if(commandLine.hasOption("P")){
			SCM.printStatistics();
			return;
		}
		
		if(commandLine.hasOption("P")){
			SCM.printStatistics();
			return;
		}
		
		if(commandLine.hasOption("E")){
			SCM.printPercentiles();
			return;
		}
		
		if(commandLine.hasOption("e")){
			SCM.extractSamples(MyFunc.readStringList2(commandLine.getOptionValue("e")));			
		}
		
		if(commandLine.hasOption("N")){
			SCM.changePosition2NearestProbe(PI);
		}
		
		if(commandLine.hasOption("B") & commandLine.hasOption("L")){
			 SCM.ternarize(Double.valueOf(commandLine.getOptionValue("B")), Double.valueOf(commandLine.getOptionValue("L")));
		}else{
			if(commandLine.hasOption("B")){
				SCM.binarize(Double.valueOf(commandLine.getOptionValue("B")));
			}
			if(commandLine.hasOption("L")){
				SCM.binarizeLess(Double.valueOf(commandLine.getOptionValue("L")));
			}
		}
		
		if(commandLine.hasOption("G") & commandLine.hasOption("K")){
			 SCM.ternarize(SCM.percentile(Double.valueOf(commandLine.getOptionValue("G"))), 
					 SCM.percentile(Double.valueOf(commandLine.getOptionValue("K"))));
		}else{
			if(commandLine.hasOption("G")){
				SCM.binarize(SCM.percentile(Double.valueOf(commandLine.getOptionValue("G"))));
			}
			if(commandLine.hasOption("K")){
				SCM.binarizeLess(SCM.percentile(Double.valueOf(commandLine.getOptionValue("K"))));
			}
		}
		
		if(commandLine.hasOption("L")){
			 SCM.binarizeLess(Double.valueOf(commandLine.getOptionValue("B")));
		}
		
		if(commandLine.hasOption("O") & PI != null){
			PrintWriter os = new PrintWriter(new FileWriter(commandLine.getOptionValue("O")));
			os.print(PI);
			os.close();
		}
		
		if(commandLine.hasOption("M") & (PI != null || GI != null)){
			if(!SCM.checkConsistency()){
				System.err.println("ERR: in SegmentContainerMap consistensiy");
				throw new DataFormatException("");
			}
			if(PI != null){
				PI.filter(SCM.get(sample.get(0)));
				SCM.toMyMat(PI).print();
			}else if(GI != null){
				GI.filter(SCM.get(sample.get(0)));
				SCM.toMyMat(GI).print();
			}
		}else{			
			SCM.print();
		}
	}
	
}