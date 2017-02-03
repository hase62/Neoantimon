package snp;

import java.util.*;
import java.io.*;
import org.apache.commons.cli.*;

import sun.reflect.Reflection;
import utility.*;

public class PACT2  {
	SegmentContainerMap SCM;
	ProbeInfo PI;
	MyMat M;
	int segLength;
	Dist C;
	Map<String, Double> P;
	Map<String, Double> minusLogP;
	List <Double> pp;
	double cutoff = 0;
	boolean countLess = false ;
	double pcutoff = 10;
	
	List<BipartiteGraph> BGs;
	List<BipartiteGraph> BCs;
	
	List<String> result;
	
	public PACT2 (SegmentContainerMap SCM,ProbeInfo PI){
		this.SCM = SCM;
		this.PI = PI;
		PI.sort();
		M = SCM.toMyMat(PI);
	}
	
	
	
	private void setCutoff(double d){
		cutoff = d;
	}
	private void setCutoffByPercentile(double d){
		if(d > 1){
			d /= 100;
		}
		cutoff = MyFunc.percentile(M.asList(), d);
	}
	
	private void calculatePoiBinParameters(){
		pp = new ArrayList<Double>();
		for(int j=0; j<M.colSize();j++){
			pp.add(0.0);
		}
		for(int j=0; j<M.colSize();j++){
			for(int i=0; i<M.rowSize();i++){
				if(M.get(i,j)>0){
					pp.set(j,pp.get(j)+1);
				}
			}
		}
		for(int j=0; j<M.colSize();j++){
			pp.set(j, Math.pow(pp.get(j)/M.rowSize(),2));
		}	
	}
	
	private void getCoaberrationCount(){
		C = new Dist(M.getRowNames());
		
		if(!countLess){
			for(int i=0; i<M.rowSize();i++){
				for(int j=0; j<i;j++){
					int tmp = 0;
					for(int k=0; k<M.colSize();k++){
						if(M.get(i,k) >= cutoff &   M.get(j,k) >= cutoff){
							tmp++;
						}
					}
					C.set(i, j, tmp);
				}
			}
		}else{
			for(int i=0; i<M.rowSize();i++){
				for(int j=0; j<i;j++){
					int tmp = 0;
					for(int k=0; k<M.colSize();k++){
						if(M.get(i,k) <= cutoff &   M.get(j,k) <= cutoff){
							tmp++;
						}
					}
					C.set(i, j, tmp);
				}
			}
			
		}
		//System.err.println(C);
	}
	
	private  void calculatePvalues(){
		P =  new  HashMap<String, Double>();
		double M = MyFunc.max(C.asList());
		double m = MyFunc.min(C.asList());
		List <Integer> kk = new ArrayList<Integer>();
		for(int i=(int) m;i<=(int)M;i++){
				kk.add(i);
		}
		List<Double> pv = PART.getPoissonBinomialPvalue(kk,pp);
		Map<Integer, Double> kk2pv = new HashMap<Integer, Double>();
		for(int i=0; i<kk.size();i++){
			kk2pv.put(kk.get(i),pv.get(i));
		}
		for(int i=0; i< C.size();i++){
			for(int j=i+1; j< C.size();j++){
				double pvalue  = kk2pv.get((int)C.get(i,j));				
				P.put(C.getNames().get(i) + "\t"+  C.getNames().get(j), pvalue);
			}
		}
	}
		
	
	
	private  void calculateMinusLogPvalues(){
		minusLogP = new LinkedHashMap<String, Double>();
		for(String p: P.keySet()){
		double tmp = P.get(p);
		if(tmp==1){
			tmp=0;
		}else{
			tmp = - Math.log10(tmp);
		}
		minusLogP.put(p, tmp);
		}
	}
	
	

	private void getBipartiteGraphs(){
		Map<String, Double> p = minusLogP;
		Set <String> edge = new HashSet <String>();
		Set <String> chr  = new HashSet <String>();
		
		for(String s: p.keySet()){
			List<String> tmp = MyFunc.split("\t", s);
			String p1 = tmp.get(0);
			String p2 = tmp.get(1);
			if(p.get(s) < pcutoff){
				continue;
			}
			if(PI.chr(p1)== PI.chr(p2)){
				continue;
			}

			if(PI.chr(p1) < PI.chr(p2)){
				chr.add(PI.chr(p1) + "\t" + PI.chr(p2));
				edge.add(p1 + "\t" + p2);
			}else{
				chr.add(PI.chr(p2) + "\t" + PI.chr(p1));
				edge.add(p2 + "\t" + p1);
			}
			
		}
		BGs = new ArrayList <BipartiteGraph>();
		for(String c :chr){
			List <String> tmp = MyFunc.split("\t", c);
			int c1 = Integer.valueOf(tmp.get(0));
			int c2 = Integer.valueOf(tmp.get(1));
			BipartiteGraph BG = new BipartiteGraph(PI.getProbeSetBychr(c1),PI.getProbeSetBychr(c2)); 
			for(String p1: BG.uNames()){
				for(String p2: BG.vNames()){
					if(edge.contains(p1 + "\t" + p2)){
						BG.add(p1, p2);
					}
				}
			}
			BGs.add(BG);
		}
	}
	
	private void getBicliques(){
		BCs = new ArrayList <BipartiteGraph>();
		for(BipartiteGraph BG: BGs){
			BCs.addAll(BG.getBicliques(0.7));	
		}
	}
	
	
	private void getResult(){
		result = new ArrayList<String>();
		for(BipartiteGraph bg: BCs){
			double d = bg.getDensity();
			double maxp = 0;
			System.err.println(bg);
			for(String u:bg.uNames()){
				for(String v:bg.vNames()){
					if(minusLogP.get(u + "\t" + v) > maxp){
						maxp = minusLogP.get(u + "\t" + v);
					}
				}
			}
			int chr = PI.chr(bg.uNames().get(0));
			int start = PI.pos(bg.uNames().get(0));
			int end = PI.pos(bg.uNames().get(bg.uSize()-1));
			int chr2 = PI.chr(bg.vNames().get(0));
			int start2 = PI.pos(bg.vNames().get(0));
			int end2 = PI.pos(bg.vNames().get(bg.vSize()-1));
			result.add(chr + "\t" + start + "\t" + end  + "\t" + chr2 + "\t" + start2 + "\t" + end2 + "\t" + maxp + "\t" + d);
		}
	}
	
	
	public void perform(){
		System.err.println("calculate parameters....");
		calculatePoiBinParameters();
		
		System.err.println("count coaberrations....");
		getCoaberrationCount();
		
		
		System.err.println("calculate pvalues....");
		calculatePvalues();
		calculateMinusLogPvalues();
		
		System.err.println("get bipartile graph....");
		getBipartiteGraphs();
		
		System.err.println("get bicleques....");
		getBicliques();
		
		System.err.println("get results....");
		getResult();
		
	}
		
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("c", "cutoff", true, "cutoff for aberration");
		options.addOption("o", "outfile", true, "output file name");
		options.addOption("t", "thinpi", true, "thin down probe info");
		options.addOption("n", "nump", true, "# of psuedo probes");
		options.addOption("l", "less", false, "count less than cutoff");
		options.addOption("p", "pcutoff", true, "pvalue cutoff");
		options.addOption("C", "cutbyper", true, "cutoff by percentile");
		options.addOption("b", "bed", false, "gene coordinate file from bed file");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine = null;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] SegFile  probeTsvFile", options);
			System.exit(1);
		}
		List <String> argList = commandLine.getArgList();
		if(!(argList.size() == 1 | argList.size() == 2)){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] SegFile probeTsvFile", options);
			System.exit(1);
		}
		SegmentContainerMap SCM = new SegmentContainerMap(argList.get(0));
		ProbeInfo PI = null;
		List<String> sample = new ArrayList<String>(SCM.keySet());
	
		if(argList.size() == 2 ){
			if(!commandLine.hasOption("b")){
				PI = 	ProbeInfo.getProbeInfoFromTsvFile(argList.get(1));
			}else{
				PI = GeneInfo.getGeneInfoFromBedFile(argList.get(1)).toProbeInfo();
			}
			if(commandLine.hasOption("t")){
				PI.thinDown(Integer.valueOf(commandLine.getOptionValue("t")));
			}
			PI.filter(SCM.get(sample.get(0)));
		}else{
			int n = 100;
			if(commandLine.hasOption("n")){
				n = Integer.valueOf(commandLine.getOptionValue("n"));
			}	
			PI = SCM.generatePsuedoProbeInfo(n);
		}
		
		PACT2  PACT = new PACT2(SCM, PI);
		
		
		if(commandLine.hasOption("c")){
			PACT.setCutoff(Double.valueOf(commandLine.getOptionValue("c")));
		}
		if(commandLine.hasOption("l")){
			PACT.countLess = true;
		}
		if(commandLine.hasOption("C")){
			PACT.setCutoffByPercentile(Double.valueOf(commandLine.getOptionValue("C")));
		}
		
		if(commandLine.hasOption("p")){
			PACT.pcutoff = Double.valueOf(commandLine.getOptionValue("p"));
		}
		
		PACT.perform();		

		
		Writer os;
		if(commandLine.hasOption("o")){
			 os = new BufferedWriter(new FileWriter(commandLine.getOptionValue("o")));	 
		}else{
			os = new PrintWriter(System.out);
		}
		os.write("chr1" + "\t" +  "start1"  + "\t" + "end1" + "\t" + "chr2" + "\t" +  "start2"  + "\t" + "end2" + "\t" +  "maxp"  + "\t" + "d");
		for(String s:PACT.result){
			os.write(s  + "\n");
		}	
		os.flush();
	}
	
	
	
}
