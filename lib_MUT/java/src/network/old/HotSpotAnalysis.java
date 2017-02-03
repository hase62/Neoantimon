package network.old;

import java.io.IOException;
import java.util.*;
import java.util.zip.DataFormatException;

import network.Link;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

import sun.reflect.Reflection;
import utility.MyFunc;


public class HotSpotAnalysis {
	
	private Link link;
	private Map <String, Double> geneWeight;
	private List <String> genes;
	private List <String> allGenes;	
	private int resamplingNumber = 1000;
	private  Map <String,Double> hotSpot;
	private  List <Double> hotSpotNull;
	private  Map <String,Double> hotSpotP;
	private HotSpotFunc hotSpotFunc;
	private String pivot = "";	
	private  Map <String,Integer> gene2degree;
	private  Map <Integer, List<String>> degree2gene;
	private  Map <Integer, List<String>> degree2geneAll;
	private double hotSpotCenterCutoff = 0;
	
	interface HotSpotFunc {
		public Map <String, Double> get (Map <String, Double> geneWeight);
	}
	
	public HotSpotAnalysis(){
		hotSpot = new HashMap<String, Double>();
		hotSpotNull = new ArrayList<Double>();
		hotSpotP = new HashMap<String, Double>();
		hotSpotFunc = new HotPairs();
		gene2degree = new HashMap<String, Integer>();
		degree2gene = new HashMap<Integer, List<String>>();
		degree2geneAll = new HashMap<Integer, List<String>>();
	}	
	
	public void setLink(Link link){
		this.link = link;
		allGenes = link.getNodeName();
	}
	public void setLink(String infile) throws IOException, DataFormatException{
		link = new Link(infile); 
		allGenes = link.getNodeName();
	}
	public void setResamplingNumber(int resamplingNumber){
		this.resamplingNumber = resamplingNumber;
	}
	public void setGeneWeight(String infile) throws IOException, DataFormatException{
		Map<String, String> tmp = MyFunc.readStringStringMap(infile);
		geneWeight = new HashMap<String, Double>();
		for(String k: tmp.keySet()){
			if(link.containsNode(k)){
				geneWeight.put(k, Double.valueOf(tmp.get(k)));
			}
		}
		genes = new ArrayList<String>(geneWeight.keySet());
	}
	
	public void setPivot(String pivot){
		if(link.containsNode(pivot)){
			this.pivot = pivot;
		}
	}
	
	public void setHotSpotFunc(int i){
		switch(i){
		case 0:
			hotSpotFunc = new HotPairs();
			break;
		case 1:
			hotSpotFunc = new HotSpots();
			break;
		case 3:
			hotSpotFunc = new HotSpotsExact();
			break;
		default:	
		}
	}
	
	public void setHotSoptCenterCutoff(double d){
		hotSpotCenterCutoff = d;
	}
	
	private  void setDegree(){
		for(String g1:allGenes){
			int tmp = 1;
			for(String g2:allGenes){
				if(!g1.equals(g2) & link.get(g1,g2)){
					tmp++;
				}
			}
			gene2degree.put(g1,tmp);
			if(degree2geneAll.containsKey(tmp)){
				degree2geneAll.get(tmp).add(g1);
			}else{
				List<String> tmp2 = new ArrayList<String>();
				tmp2.add(g1);
				degree2geneAll.put(tmp,tmp2);
			}
			if(genes.contains(g1)){
				if(degree2gene.containsKey(tmp)){
					degree2gene.get(tmp).add(g1);
				}else{
					List<String> tmp2 = new ArrayList<String>();
					tmp2.add(g1);
					degree2gene.put(tmp,tmp2);
				}
			}
			
		}
	}
	private class HotPairs implements HotSpotFunc {
		public String toString(){
			return "HotPairs";
		}
		public  Map <String, Double> get(Map <String, Double> geneWeight){
			if(pivot.equals("")){
				return getWithoutPivot( geneWeight);
			}else{
				return getWithPivot( geneWeight);
			}
		}
		
		private Map <String, Double> getWithoutPivot(Map <String, Double> geneWeight){
			List <String> genes = new ArrayList<String>(geneWeight.keySet());
			Map <String,Double> hotPairs = new LinkedHashMap<String, Double>();
			for(String g1: genes){
				for(String g2: genes){
					if(g1.compareTo(g2)<0){
						if(link.get(g1,g2)){
							hotPairs.put(g1 + ":" + g2, geneWeight.get(g1)+geneWeight.get(g2));
						}
					}
				}
			}
			return hotPairs;
		}
		
		private  Map <String, Double> getWithPivot(Map <String, Double> geneWeight){
			List <String> genes = new ArrayList<String>(geneWeight.keySet());
			Map <String, Double> hotPairs = new LinkedHashMap<String,  Double>();
			for(String g1: genes){
				if(!g1.equals(pivot) && link.get(g1,pivot)){
					hotPairs.put(g1 + ":" + pivot, geneWeight.get(g1));
				}
			}
			return hotPairs;
		}
		
	}
	
	private class HotSpots implements HotSpotFunc{
		public String toString(){
			return "HotSpots";
		}
		public  Map <String, Double> get(Map <String, Double> geneWeight){
			if(pivot.equals("")){
				return getWithoutPivot( geneWeight);
			}else{
				return getWithPivot( geneWeight);
			}	
		}
		private Map <String, Double> getWithoutPivot(Map <String, Double> geneWeight){
			List <String> genes = new ArrayList<String>(geneWeight.keySet());
			Map <String,Double> hotSpots = new HashMap<String, Double>();
			Set <String> tmp = new HashSet<String>();
			for(String g1: genes){
				if(geneWeight.get(g1) >= hotSpotCenterCutoff){
					tmp.add(g1);
				}
			}
			for(String g1: tmp){
				Set <String> hot = new HashSet<String>();
				double w = 0;
				hot.add(g1);
				w += geneWeight.get(g1);
				for(String g2: link.getNeighbors(g1)){
					if(geneWeight.containsKey(g2) & !hot.contains(g2)){
						hot.add(g2);
						w += geneWeight.get(g2);
					}
				}
				if(hot.size() == 1){
					continue;
				}
				List <String> hot2 = new ArrayList<String>(hot);
				Collections.sort(hot2);
				String tmp2 = MyFunc.join(":", hot2);
				hotSpots.put(tmp2, w);
			}
			return hotSpots;
		}
		
		private Map <String, Double> getWithPivot(Map <String, Double> geneWeight){
			Map <String,Double> hotSpots = new HashMap<String, Double>();
			Set <String> tmp = new HashSet<String>();
			tmp.add(pivot);
			for(String g2: link.getNeighbors(pivot)){
				tmp.add(g2);
			}
			for(String g1: tmp){
				Set <String> hot = new HashSet<String>();
				hot.add(pivot);
				double w = 0;
				if(geneWeight.containsKey(g1)){
					hot.add(g1);
					if(!g1.equals(pivot)){
						w += geneWeight.get(g1);
					}
				}
				for(String g2: link.getNeighbors(g1)){
					if(geneWeight.containsKey(g2) & !hot.contains(g2)){
						hot.add(g2);
						if(!g2.equals(pivot)){
							w += geneWeight.get(g2);
						}
					}
				}
				
				if(hot.size() == 1){
					continue;
				}
				List <String> hot2 = new ArrayList<String>(hot);
				Collections.sort(hot2);
				String tmp2 = MyFunc.join(":", hot2);
				hotSpots.put(tmp2, w);
			}
			return hotSpots;
		}
		
	}
	
	private class HotSpotsExact implements HotSpotFunc{
		public String toString(){
			return "HotSpots";
		}
		public  Map <String, Double> get(Map <String, Double> geneWeight){
			if(pivot.equals("")){
				return getWithoutPivot( geneWeight);
			}else{
				return getWithPivot(geneWeight);
			}	
		}
		private Map <String, Double> getWithoutPivot(Map <String, Double> geneWeight){
			List <String> genes = new ArrayList<String>(geneWeight.keySet());
			Map <String,Double> hotSpots = new HashMap<String, Double>();
			Set <String> tmp = new HashSet<String>(genes);
			for(String g1: genes){
				for(String g2: link.getNeighbors(g1)){
					tmp.add(g2);
				}
			}
			for(String g1: tmp){
				Set <String> hot = new HashSet<String>();
				double w = 0;
				if(geneWeight.containsKey(g1)){
					hot.add(g1);
					w += geneWeight.get(g1);
				}
				for(String g2: link.getNeighbors(g1)){
					if(geneWeight.containsKey(g2) & !hot.contains(g2)){
						hot.add(g2);
						w += geneWeight.get(g2);
					}
				}
				if(hot.size() == 1){
					continue;
				}
				List <String> hot2 = new ArrayList<String>(hot);
				Collections.sort(hot2);
				String tmp2 = MyFunc.join(":", hot2);
				hotSpots.put(tmp2, w);
			}
			return hotSpots;
		}
		
		private Map <String, Double> getWithPivot(Map <String, Double> geneWeight){
			Map <String,Double> hotSpots = new HashMap<String, Double>();
			Set <String> tmp = new HashSet<String>();
			tmp.add(pivot);
			for(String g2: link.getNeighbors(pivot)){
				tmp.add(g2);
			}
			for(String g1: tmp){
				Set <String> hot = new HashSet<String>();
				hot.add(pivot);
				double w = 0;
				if(geneWeight.containsKey(g1)){
					hot.add(g1);
					if(!g1.equals(pivot)){
						w += geneWeight.get(g1);
					}
				}
				for(String g2: link.getNeighbors(g1)){
					if(geneWeight.containsKey(g2) & !hot.contains(g2)){
						hot.add(g2);
						if(!g2.equals(pivot)){
							w += geneWeight.get(g2);
						}
					}
				}
				
				if(hot.size() == 1){
					continue;
				}
				List <String> hot2 = new ArrayList<String>(hot);
				Collections.sort(hot2);
				String tmp2 = MyFunc.join(":", hot2);
				hotSpots.put(tmp2, w);
			}
			return hotSpots;
		}
		
	}
	
	private void gethotSpots(){
		hotSpot = hotSpotFunc.get(geneWeight);
		sortResult();
	}
	
	private void sortResult(){
		hotSpot = sortMap(hotSpot);
	}
	private  static Map <String,Double> sortMap(Map <String,Double> map){
		Map<Double,List<String>> newMap = new HashMap<Double,List<String>>();
		for(Map.Entry<String,Double> e : map.entrySet()){
			if(newMap.containsKey(e.getValue())){
				(newMap.get(e.getValue())).add(e.getKey());
			}else{
				List<String> tmp = new ArrayList<String>();
				tmp.add(e.getKey());
				newMap.put(e.getValue(), tmp);
			}
		}
		List <Double> tmp = new ArrayList<Double>(newMap.keySet());
		Collections.sort(tmp);
		Collections.reverse(tmp);
		Map <String, Double> sortedMap = new LinkedHashMap<String, Double>();
		for(Double i: tmp){
			for(String s: newMap.get(i)){
				sortedMap.put(s, i);
			}
		}
		return sortedMap;
	}
	public String toString(){
		List <String> tmp = new ArrayList<String>();
		for(String s: hotSpot .keySet()){
			tmp.add(s + "\t" + hotSpot .get(s) + "\t" + hotSpotP.get(s));
		}
		if(tmp.size()==0){
			return "no gene pair!\n";
		}else{
			return MyFunc.join("\n", tmp) + "\n";
		}
	}
	public void perform(){
		gethotSpots();
		generateNullDist();
		getPvalues();
	}
	
	private void generateNullDist(){
		setDegree();
		for(int i = 0; i < resamplingNumber; i++){
			System.err.println(i);
			Map <String, Double>  nullGeneWeight =getNullGeneWeight();
			Map <String, Double>  nullResult = hotSpotFunc.get(nullGeneWeight);
			double nullMaxValue = 0.0;
			if(!nullResult.isEmpty()){
				nullMaxValue = max(nullResult.values());
			}
			hotSpotNull.add(nullMaxValue);
		}
	}

	private double max(Collection <Double> v ){
		double  m = 0;
		for(Double i : v){
			if(i > m){
				m = i;
			}
			
		}
		return m;
	}
	
	private Map <String, Double> getNullGeneWeight(){
		 Map <String, Double> nullGeneWeight = new HashMap<String, Double>();
		 for(Integer i: degree2gene.keySet()){
			 List <String> g = degree2gene.get(i);
			int n = g.size();
			List <String> nullg = MyFunc.sample(degree2geneAll.get(i), n);
			for(int j = 0; j < n; j++){
				if(pivot.equals(g.get(j))){
					nullGeneWeight.put(pivot, geneWeight.get(g.get(j)));
				}else{
					nullGeneWeight.put(nullg.get(j), geneWeight.get(g.get(j)));
				}
			}
		}
		return nullGeneWeight;
	}
	
	private void getPvalues(){
		for(String s: hotSpot.keySet()){
			 double count = hotSpot.get(s);
			 double p = 0.0;
			 for(Double i: hotSpotNull){
				 if(i >= count){
					 p++;
				 }
			 }
			 if(p==0.0){
				 p=1.0;
			 }
			p /= hotSpotNull.size();
			hotSpotP.put(s, p);
		}
	}
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("p", "pivot", true,  "geneName");
		options.addOption("r", "resamp", true,  "# of resamplings for P value calculation (default:1000)");
		options.addOption("0", "nbr", false,  "search HotPairs (default)");
		options.addOption("1", "snbr", false,  "search HotSpot");
		options.addOption("2", "dmtr2", false,  "search HotSpot exactly");
		
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options]  linkFile geneCountFile", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 2){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] linkFile geneCountFile", options);
			return;
		}
		HotSpotAnalysis M = new HotSpotAnalysis();
		M.setLink(argList.get(0));
		M.setGeneWeight(argList.get(1));
		
		if(commandLine.hasOption("0")){
			M.setHotSpotFunc(0);
		}
		if(commandLine.hasOption("1")){
			M.setHotSpotFunc(1);
		}
		if(commandLine.hasOption("2")){
			M.setHotSpotFunc(2);
		}
		if(commandLine.hasOption("p")){
			M.setPivot(commandLine.getOptionValue("p"));
		}
		if(commandLine.hasOption("r")){
			M.setResamplingNumber(Integer.valueOf(commandLine.getOptionValue("r")));
		}
		M.perform();
		System.out.print(M);
	
	}
	
	
	
}
