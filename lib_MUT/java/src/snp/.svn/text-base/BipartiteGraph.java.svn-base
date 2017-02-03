package snp;

import java.util.*;

import org.apache.commons.cli.*;

import sun.reflect.Reflection;
import utility.*;

public class BipartiteGraph {
	MyMat G;
	
	BipartiteGraph(MyMat G){
		this.G = G;
	}
	
	BipartiteGraph(BipartiteGraph BG){
		G = new MyMat(BG.G);
	}
	
	BipartiteGraph(List<String> U, List<String> V){
		G = new MyMat(U, V);
	}
	
	public BipartiteGraph getSubGraph(List<String> U, List<String> V){
		return new BipartiteGraph(G.getSubMatrix(U, V)); 
	}
	
	public int uSize(){
		return G.rowSize();
	}
	
	public int vSize(){
		return G.colSize();
	}
	
	public List<String> uNames(){
		return G.getRowNames();
	}
	
	public List<String> vNames(){
		return G.getColNames();
	}
	
	public void add(String u, String v){
		G.set(u,v,1);
	}
	
	public void remove(String u, String v){
		G.set(u,v,0);
	}
	
	boolean get(String u, String v){
		return G.get(u,v)==1?true:false;
	}
	
	void add(int i, int j){
		G.set(i,j,1);
	}
	
	void remove(int i, int j){
		G.set(i,j,0);
	}
	
	void remove(List<String> U, List<String> V){
		for(String u: U){
			for(String v: V){
				G.set(u,v,0);
			}
		}
	}
	
	
	public boolean get(int i, int j){
		return G.get(i,j)==1?true:false;
	}
	
	public String toString(){
		return G.toString();
	}
	
	public double getDensity(){
		return MyFunc.sum(G.asList())/(uSize()*vSize());
	}
	
	private double getDensitiy(int i_start, int i_end, int j_start, int j_end){
		int n = (i_end -i_start+1)*(j_end-j_start+1);
		double k = 0;
		for(int i = i_start; i<=i_end; i++){
			for(int j = j_start; j<=j_end; j++){
				if(get(i,j)){
					k++;
				}
			}
		}
			return k/n;
	}
	
	private double getDensitiy(List<String> U, List<String>V){
		int n = U.size()*V.size();
		double k = 0;
		for(String u:U){
			for(String v:V){
				if(get(u,v)){
					k++;
				}
			}
		}
			return k/n;
	}
		
	private  BipartiteGraph getBicliqueGreedily(double cutoff){
		int i_start = 0, i_end = uSize()-1, j_start = 0, j_end = vSize()-1;
		double  d = getDensitiy(i_start,  i_end, j_start,j_end);
		
		while(true){
			boolean update = false; 
			while(true){
				double nextd = getDensitiy(i_start+1,  i_end, j_start, j_end);
				if(nextd > d & i_start != i_end){
					i_start++;
					d = nextd;
					update = true;
					if(d >= cutoff){
						break;
					}
				}else{
					break;
				}
			}
			while(true){
				double nextd = getDensitiy(i_start,  i_end-1, j_start,j_end);
				if(nextd > d & i_start != i_end){
					i_end--;
					d = nextd;
					update = true;
					if(d >= cutoff){
						break;
					}
				}else{
					break;
				}
			}
			while(true){
				double nextd = getDensitiy(i_start,  i_end, j_start+1,j_end);
				if(nextd > d & j_start != j_end){
					j_start++;
					d = nextd;
					update = true;
					if(d >= cutoff){
						break;
					}
				}else{
					break;
				}
			}
			while(true){
				double nextd = getDensitiy(i_start,  i_end, j_start,j_end-1);
				if(nextd > d & j_start != j_end){
					j_end--;
					d = nextd;
					update = true;
					if(d >= cutoff){
						break;
					}
				}else{
					break;
				}
			}
			if(!update){
			 break;
			}
		}
		return getSubGraph(uNames().subList(i_start, i_end+1), vNames().subList(j_start, j_end+1));
	}
	
	private  BipartiteGraph getBiclique(double cutoff){
		
		int i_start = 0, i_end = uSize()-1, j_start = 0, j_end = vSize()-1;
		
		while(true){
			double nextd = getDensitiy(0,  i_start, 0, vSize()-1);
			if(nextd==0 & i_start < uSize()-1 & i_end > i_start){
				i_start++;
			}else{
				break;
			}
		}
		
		while(true){
			double nextd = getDensitiy(i_end,  uSize()-1, 0, vSize()-1);
			if(nextd==0 & i_end > 0 & i_end > i_start){
				i_end--;
			}else{
				break;
			}
		}
		
		while(true){
			double nextd = getDensitiy(0,  uSize()-1, 0, j_start);
			if(nextd==0 & j_start < vSize()-1 &  j_end > j_start){
				j_start++;
			}else{
				break;
			}
		}
		
		while(true){
			double nextd = getDensitiy(0, uSize()-1, j_end, vSize()-1);
			if(nextd==0 & j_end >0 & j_end > j_start){
				j_end--;
			}else{
				break;
			}
		}
		
		List<List<Integer>> I = new ArrayList<List<Integer>>();
		for(int i = i_start; i <= i_end; i++){
			for(int j = i; j <= i_end; j++){
				List<Integer> tmp = new ArrayList<Integer>();
				tmp.add(i);
				tmp.add(j);
				I.add(tmp);
			}
		}	
		List<List<Integer>> J = new ArrayList<List<Integer>>();
		for(int i = j_start; i <= j_end; i++){
			for(int j = i; j <= j_end; j++){
				List<Integer> tmp = new ArrayList<Integer>();
				tmp.add(i);
				tmp.add(j);
				J.add(tmp);
			}
		}	
		
		double  d = -1;
		double size = -1;
		
		for(List<Integer> i:I){
			for(List<Integer> j:J){
				if(d < -1){
					d = getDensitiy(i.get(0),  i.get(1), j.get(0), j.get(1)); 
					size = (i.get(1)-i.get(0)+1)*(j.get(1)-j.get(0)+1); 
					i_start = i.get(0);
					i_end = i.get(1); 
					j_start = j.get(0);
					j_end = j.get(1);
				}else{
					double nextd = getDensitiy(i.get(0),  i.get(1), j.get(0), j.get(1));
					double nextsize = (i.get(1)-i.get(0)+1)*(j.get(1)-j.get(0)+1);
					if(nextd >= cutoff  & nextsize > size){
						d = nextd;
						size = nextsize; 
						i_start = i.get(0);
						i_end = i.get(1); 
						j_start = j.get(0);
						j_end = j.get(1);
					}
				}
			}
		}
		
		return getSubGraph(uNames().subList(i_start, i_end+1), vNames().subList(j_start, j_end+1));	
		
		
		
	}
	
	
	public  List <BipartiteGraph> getBicliques(double cutoff){
		BipartiteGraph BP = new BipartiteGraph(this);
		List<BipartiteGraph> list = new ArrayList<BipartiteGraph>(); 
		while(true){
			BipartiteGraph BC = BP.getBicliqueGreedily(cutoff);
			//BipartiteGraph BC = BP.getBiclique(cutoff);
			if(BC.getDensity() < cutoff){
				return list;
			}
			list.add(BC);
			BP.remove(BC.uNames(), BC.vNames());
		}
	}
	
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine = null;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] tabFile", options);
			System.exit(1);
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() < 1){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] tabFile", options);
			System.exit(1);
		}
	
		BipartiteGraph BG  = new BipartiteGraph(new MyMat(argList.get(0)));
		
		System.out.println(BG);
		for(BipartiteGraph bg: BG.getBicliques(0.8)){
			System.out.println(bg);
		}
	}
	

}
