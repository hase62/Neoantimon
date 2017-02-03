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


public class MutationHotSpotAnalysis {
	private Link link;
	private Map <String, Integer> geneCount;
	private List <String> genes;
	private List <String> allGenes;	
	private int resamplingNumber = 1000;
	private  Map <String,Integer> hotSpot;
	private  List <Integer> hotSpotNull;
	private  Map <String,Double> hotSpotP;
	private HotSpotFunc hotSpotFunc;
	private String pivot = "";	
	
	public MutationHotSpotAnalysis(){
		hotSpot = new HashMap<String, Integer>();
		hotSpotNull = new ArrayList<Integer>();
		hotSpotP = new HashMap<String, Double>();
		hotSpotFunc = new HotNeighborPairs();
	}	
	
	public void setLink(Link link){
		this.link = link;
		allGenes = link.getNodeName();
	}
	public void setLink(String infile) throws IOException, DataFormatException{
		link = new Link(infile); 
		allGenes = link.getNodeName();
	}
	public void setGeneCount(String infile) throws IOException, DataFormatException{
		Map<String, String> tmp = MyFunc.readStringStringMap(infile);
		geneCount = new HashMap<String, Integer>();
		for(String k: tmp.keySet()){
			if(link.containsNode(k)){
				geneCount.put(k, Integer.valueOf(tmp.get(k)));
			}
		}
		genes = new ArrayList<String>(geneCount.keySet());
	}
	
	public void setPivot(String pivot){
		if(link.containsNode(pivot)){
			this.pivot = pivot;
		}
	}
	
	public void setHotSpotFunc(int i){
		switch(i){
		case 0:
			hotSpotFunc = new HotNeighborPairs();
			break;
		case 1:
			hotSpotFunc = new HotSemiNeighborPairs();
			break;
		case 2:
			hotSpotFunc = new HotSpotsOfDiameter2();
			break;
		case 3:
			hotSpotFunc = new HotSpotsOfDiameter3();
			break;
		default:	
		}
	}
	
	public void setResamplingNumber(int resamplingNumber){
		this.resamplingNumber = resamplingNumber;
	}
	
	interface HotSpotFunc {
		public Map <String, Integer> get (Map <String, Integer> geneCount);
		public Map <String, Integer> get (Map <String, Integer> geneCount, String pivot);
	}
	
	private class HotNeighborPairs implements HotSpotFunc {
		public String toString(){
			return "HotNeighborPairs";
		}
		public  Map <String, Integer> get(Map <String, Integer> geneCount){
			List <String> genes = new ArrayList<String>(geneCount.keySet());
			Map <String,Integer> hotPairs = new LinkedHashMap<String, Integer>();
			for(String g1: genes){
				for(String g2: genes){
					if(g1.compareTo(g2)<0){
						if(link.get(g1,g2)){
							hotPairs.put(g1 + ":" + g2, geneCount.get(g1)+geneCount.get(g2));
						}
					}
				}
			}
			return hotPairs;
		}
		
		public  Map <String, Integer> get(Map <String, Integer> geneCount, String pivot){
			
			List <String> genes = new ArrayList<String>(geneCount.keySet());
			Map <String,Integer> hotPairs = new LinkedHashMap<String, Integer>();
			for(String g1: genes){
				if(!g1.equals(pivot) && link.get(g1,pivot)){
					hotPairs.put(g1 + ":" + pivot, geneCount.get(g1));
				}
			}
			return hotPairs;
		}
	}
	
	private class HotSemiNeighborPairs implements HotSpotFunc {
		public String toString(){
			return "HotSemiNeighborPairs";
		}
		public  Map <String, Integer> get(Map <String, Integer> geneCount){
			List <String> genes = new ArrayList<String>(geneCount.keySet());
			Map <String,Integer> hotPairs = new LinkedHashMap<String, Integer>();
			for(String g1: genes){
				for(String g2: genes){
					if(g1.compareTo(g2)<0){
						if(link.havePath(g1,g2,2)){
							hotPairs.put(g1 + ":" + g2, geneCount.get(g1)+geneCount.get(g2));
						}
					}
				}
			}
			return hotPairs;
		}
		
		public  Map <String, Integer> get(Map <String, Integer> geneCount, String pivot){
			List <String> genes = new ArrayList<String>(geneCount.keySet());
			Map <String,Integer> hotPairs = new LinkedHashMap<String, Integer>();
			for(String g1: genes){
				if(!g1.equals(pivot) && (link.get(g1,pivot) || link.havePath(g1,pivot,2))){
					hotPairs.put(g1 + ":" + pivot, geneCount.get(g1));
				}
			}
			return hotPairs;
		}
	}
	
	private class HotSpotsOfDiameter2 implements HotSpotFunc{
		public String toString(){
			return "HotSpotsOfDiameter2";
		}
		public  Map <String, Integer> get(Map <String, Integer> geneCount){
			List <String> genes = new ArrayList<String>(geneCount.keySet());
			Map <String,Integer> hotSpotsOfDiameter2 = new HashMap<String, Integer>();
			Set <String> tmp = new HashSet<String>(genes);
			for(String g1: genes){
				for(String g2: link.getNeighbors(g1)){
					tmp.add(g2);
				}
			}
			for(String g1: tmp){
				Set <String> hot = new HashSet<String>();
				int count = 0;
				if(geneCount.containsKey(g1)){
					hot.add(g1);
					count += geneCount.get(g1);
				}
				for(String g2: link.getNeighbors(g1)){
					if(geneCount.containsKey(g2) & !hot.contains(g2)){
						hot.add(g2);
						count += geneCount.get(g2);
					}
				}
				if(hot.size() == 1){
					continue;
				}
				List <String> hot2 = new ArrayList<String>(hot);
				Collections.sort(hot2);
				String tmp2 = MyFunc.join(":", hot2);
				hotSpotsOfDiameter2.put(tmp2, count);
			}
			return hotSpotsOfDiameter2;
			
		}
		
		
		public Map <String, Integer> get(Map <String, Integer> geneCount, String pivot){
			Map <String,Integer> hotSpotsOfDiameter2 = new HashMap<String, Integer>();
			Set <String> tmp = new HashSet<String>();
			tmp.add(pivot);
			for(String g2: link.getNeighbors(pivot)){
				tmp.add(g2);
			}
			for(String g1: tmp){
				Set <String> hot = new HashSet<String>();
				hot.add(pivot);
				int count = 0;
				if(geneCount.containsKey(g1)){
					hot.add(g1);
					if(!g1.equals(pivot)){
						count += geneCount.get(g1);
					}
				}
				for(String g2: link.getNeighbors(g1)){
					if(geneCount.containsKey(g2) & !hot.contains(g2)){
						hot.add(g2);
						if(!g2.equals(pivot)){
							count += geneCount.get(g2);
						}
					}
				}
				
				if(hot.size() == 1){
					continue;
				}
				List <String> hot2 = new ArrayList<String>(hot);
				Collections.sort(hot2);
				String tmp2 = MyFunc.join(":", hot2);
				hotSpotsOfDiameter2.put(tmp2, count);
			}
			return hotSpotsOfDiameter2;
			
		}
		
	}
	
	private class HotSpotsOfDiameter3 implements HotSpotFunc{
		public String toString(){
			return "HotSpotsOfDiameter3";
		}
		public  Map <String, Integer> get(Map <String, Integer> geneCount){
			List <String> genes = new ArrayList<String>(geneCount.keySet());
			Map <String,Integer> hotSpotsOfDiameter3 = new HashMap<String, Integer>();
			Set <String> tmp = new HashSet<String>(genes);
			for(String g1: genes){
				tmp.addAll(link.getNeighbors(g1));
			}
			
			Set <String[]> corePair = new HashSet<String[]>();
			for(String g1: tmp){
				for(String g2: tmp){
					if(g1.compareTo(g2)<0){
						if(link.get(g1, g2)){
							String[] tmp2 = {g1, g2}; 
							corePair.add(tmp2);
						}
					}
				}
			}
			for(String[] cp: corePair){
				Set <String> hot = new HashSet<String>();
				int count = 0;
				if(geneCount.containsKey(cp[0])){
					hot.add(cp[0]);
					count += geneCount.get(cp[0]);
				}
				if(geneCount.containsKey(cp[1])){
					hot.add(cp[1]);
					count += geneCount.get(cp[1]);
				}
				for(String g1: link.getNeighbors(cp[0])){
					if(geneCount.containsKey(g1) & !hot.contains(g1)){
						hot.add(g1);
						count += geneCount.get(g1);
					}
				}
				for(String g2: link.getNeighbors(cp[1])){
					if(geneCount.containsKey(g2) & !hot.contains(g2)){
						hot.add(g2);
						count += geneCount.get(g2);
					}
				}
				List <String> hot2 = new ArrayList<String>(hot);
				Collections.sort(hot2);
				String tmp2 = MyFunc.join(":", hot2);
				if(hot.size() == 1){
					continue;
				}
				hotSpotsOfDiameter3.put(tmp2, count);
			}
			return hotSpotsOfDiameter3;
		}
		

		public  Map <String, Integer> get(Map <String, Integer> geneCount, String pivot){
			List <String> genes = new ArrayList<String>(geneCount.keySet());
			Map <String,Integer> hotSpotsOfDiameter3 = new HashMap<String, Integer>();
			Set <String> tmp = new HashSet<String>(genes);
			tmp.add(pivot);
			tmp.addAll(link.getNeighbors(pivot));
			
			Set <String[]> corePair = new HashSet<String[]>();
			for(String g1: tmp){
				for(String g2: tmp){
					if(g1.compareTo(g2)<0){
						if(link.get(g1, g2)){
							String[] tmp2 = {g1, g2}; 
							corePair.add(tmp2);
						}
					}
				}
			}
			for(String[] cp: corePair){
				Set <String> hot = new HashSet<String>();
				int count = 0;
				
				if(cp[0].equals(pivot)){
					hot.add(cp[0]);
				}else if(geneCount.containsKey(cp[0])){
					hot.add(cp[0]);
					count += geneCount.get(cp[0]);
				}
				
				if(cp[1].equals(pivot)){
					hot.add(cp[1]);
				}else if(geneCount.containsKey(cp[1])){
					hot.add(cp[1]);
					count += geneCount.get(cp[1]);
				}
				
				for(String g1: link.getNeighbors(cp[0])){
					
					if(g1.equals(pivot)){
						hot.add(g1);
					}else if(geneCount.containsKey(g1) & !hot.contains(g1)){
						hot.add(g1);
						count += geneCount.get(g1);
					}
				}
				for(String g2: link.getNeighbors(cp[1])){
					if(g2.equals(pivot)){
						hot.add(g2);
					}else if(geneCount.containsKey(g2) & !hot.contains(g2)){
						hot.add(g2);
						count += geneCount.get(g2);
					}
				}
				if(!hot.contains(pivot)){
					continue;
				}
				List <String> hot2 = new ArrayList<String>(hot);
				Collections.sort(hot2);
				String tmp2 = MyFunc.join(":", hot2);
				if(hot.size() == 1){
					continue;
				}
				hotSpotsOfDiameter3.put(tmp2, count);
			}
			return hotSpotsOfDiameter3;
		}
		
	}
	
	
	
	private void gethotSpots(){
		hotSpot = hotSpotFunc.get(geneCount);
		sortResult();
	}
	
	private void gethotSpots(String pivot){
		hotSpot = hotSpotFunc.get(geneCount, pivot);
		sortResult();
	}
	
	private Map <String, Integer> getNullGeneCount(){
	   List <String> tmp =	MyFunc.sampleWithReplacement(allGenes, genes.size());
	   Map <String, Integer> tmp2 = new HashMap<String, Integer>();
	   for(String s: tmp){
		   if(tmp2.containsKey(s)){
			   tmp2.put(s, tmp2.get(s)+1);
		   }else{
			   tmp2.put(s, 1);
		   }
	   }
	   return tmp2;
	}
	
	private Map <String, Integer> getNullGeneCount(String pivot){
			int k;
			if(geneCount.containsKey(pivot)){
			k = geneCount.get(pivot);
			}else{
			k = 0;	
			}
		   List <String> tmp =	MyFunc.sampleWithReplacement(allGenes, genes.size() - k);
		   Map <String, Integer> tmp2 = new HashMap<String, Integer>();
		   for(String s: tmp){
			   if(tmp2.containsKey(s)){
				   tmp2.put(s, tmp2.get(s)+1);
			   }else{
				   tmp2.put(s, 1);
			   }
		   }
		   if(geneCount.containsKey(pivot)){
		   tmp2.put(pivot, k);
		   }
		   return tmp2;
		}
	
	private void generateNullDist(){
		for(int i = 0; i < resamplingNumber; i++){
			Map <String, Integer>  nullgeneCount = getNullGeneCount();
			hotSpotNull.add(max(hotSpotFunc.get(nullgeneCount).values()));
		}
	}
	private void generateNullDist(String pivot){
		for(int i = 0; i < resamplingNumber; i++){
			Map <String, Integer>  nullgeneCount = getNullGeneCount(pivot);
			hotSpotNull.add(max(hotSpotFunc.get(nullgeneCount,pivot).values()));
		}
	}
	
	private int max(Collection <Integer> v ){
		int  m = 0;
		for(Integer i : v){
			if(i > m){
				m = i;
			}
			
		}
		return m;
	}
	
	
	private void getPvalues(){
		for(String s: hotSpot.keySet()){
			 int count = hotSpot.get(s);
			 double p = 0.0;
			 for(Integer i: hotSpotNull){
				 if(i >= count){
					 p++;
				 }
			 }
			p /= hotSpotNull.size();
			hotSpotP.put(s, p);
		}
	}
	
	
	private void sortResult(){
		hotSpot = sortMap(hotSpot);
	}
	
	
	private  static Map <String,Integer> sortMap(Map <String,Integer> map){
		Map<Integer,List<String>> newMap = new HashMap<Integer,List<String>>();
		for(Map.Entry<String,Integer> e : map.entrySet()){
			if(newMap.containsKey(e.getValue())){
				(newMap.get(e.getValue())).add(e.getKey());
			}else{
				List<String> tmp = new ArrayList<String>();
				tmp.add(e.getKey());
				newMap.put(e.getValue(), tmp);
			}
		}
		List <Integer> tmp = new ArrayList<Integer>(newMap.keySet());
		Collections.sort(tmp);
		Collections.reverse(tmp);
		Map <String, Integer> sortedMap = new LinkedHashMap<String, Integer>();
		for(Integer i: tmp){
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
		if(pivot == ""){
			gethotSpots();
			generateNullDist();
		}else{
			gethotSpots(pivot);
			generateNullDist(pivot);
		}
		getPvalues();
	}
	
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("p", "pivot", true,  "geneName");
		options.addOption("r", "resamp", true,  "# of resamplings for P value calculation (default:1000)");
		options.addOption("0", "nbr", false,  "search HotNeighborPairs (default)");
		options.addOption("1", "snbr", false,  "search HotSemNeighborPairs");
		options.addOption("2", "dmtr2", false,  "search HotSpotsOfDiameter2");
		options.addOption("3", "dmtr3", false,  "search HotSpotsOfDiameter3");
		
		
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
		MutationHotSpotAnalysis M = new MutationHotSpotAnalysis();
		M.setLink(argList.get(0));
		M.setGeneCount(argList.get(1));
		
		if(commandLine.hasOption("0")){
			M.setHotSpotFunc(0);
		}
		if(commandLine.hasOption("1")){
			M.setHotSpotFunc(1);
		}
		if(commandLine.hasOption("2")){
			M.setHotSpotFunc(2);
		}
		if(commandLine.hasOption("3")){
			M.setHotSpotFunc(3);
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
