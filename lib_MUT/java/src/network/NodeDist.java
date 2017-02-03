package network;

import java.io.*;
import java.util.*;
import java.util.zip.DataFormatException;

import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;
import utility.*;

public class NodeDist implements Serializable{
	private static final long serialVersionUID = 1L;
	private int M[][];
	private int n;
	private Map<String, Integer> name2index;
	private List<String> name;
	public NodeDist(List <String> name){
		n = name.size();
		this.name = new ArrayList<String>(name);
		name2index = new HashMap<String, Integer>();
		int i;
		for(i =0; i< n; i++){
			name2index.put(name.get(i),i);
		}
		M = new int[n][];
		for(i = 0; i < n; i++){
			M[i] = new int[i];
		}
		for(i=0;i<n;i++){
			for(int j=0;j<i;j++){
		      M[i][j] = Integer.MAX_VALUE;
		    }	
		}
	}
	public NodeDist(NodeDist D){
		n = D.n;
		name = new ArrayList<String>(D.name);
		name2index = new HashMap<String, Integer>(D.name2index);
		int i;
		M = new int[n][];
		for(i = 0; i < n; i++){
			M[i] = new int[i];
		}
		for(i=0;i<n;i++){
			for(int j=0;j<i;j++){
		      M[i][j] = D.M[i][j];
		    }	
		}
	}
	
	public NodeDist(String infile) throws IOException, DataFormatException{
		name2index = new HashMap<String, Integer>();
		BufferedReader inputStream = new BufferedReader(new FileReader(infile));
		String line;
		Set <String>  name_set = new HashSet<String>();
		String[] tmp = new String[3];
 		while((line = inputStream.readLine()) != null){
			if(line.charAt(0) == '#'){
				continue;				
			}
			 tmp = line.split("\t");
			if(tmp.length != 3){
				throw new DataFormatException("Dist: file format is wrong!");
			}
			name_set.add(tmp[0]);
			name_set.add(tmp[1]);
		}
		name  = new ArrayList<String>(name_set);
		Collections.sort(name);
		int i,j;
		n = name.size();
		for(i=0;i<n;i++){
			name2index.put(name.get(i),i);
		}
		M = new int[n][];
		boolean seen[][] = new boolean[n][];
		for(i = 0; i < n; i++){
			M[i] = new int[i];
			seen[i] = new boolean[i];
			
		}
		inputStream = new BufferedReader(new FileReader(infile));
		while((line = inputStream.readLine()) != null){
			if(line.charAt(0) == '#'){
				continue;				
			}
			tmp = line.split("\t");
			i = name2index.get(tmp[0]);
			j = name2index.get(tmp[1]);
			if(i > j){
				M[i][j] = (Integer.valueOf(tmp[2]));
				seen[i][j] = true;
			}
			if(j > i){
				M[j][i] = (Integer.valueOf(tmp[2]));
				seen[j][i] = true;
			}
		}
		for(i=0;i<n;i++){
		    for(j=0;j<i;j++){
		      if(seen[i][j] == false){
		    	  throw new DataFormatException("Dist: file format is wrong!");
		      }
		    }
		  }
	   inputStream.close();
	}
	public Collection<String> getNames(){
		return name;
	}
	public int size(){
		return n;
	}
	public boolean containsName(String name){
		return name2index.containsKey(name);
	}
	public int get(int i, int j){
		if(i >= n || j >= n){
			throw new IndexOutOfBoundsException();
		}
		if(i > j){
			return M[i][j];
		}
		if(j > i){
			return M[j][i];
		}
		return 0;
	}
	public int get(String s, String t){
		int i = name2index.get(s);
		int j = name2index.get(t);
		if(i > j){
			return M[i][j];
		}
		if(j > i){
			return M[j][i];
		}
		return 0;
	}
	public void set(String s, String t, int d){
		int i = name2index.get(s);
		int j = name2index.get(t);
		if(i > j){
			M[i][j] = d;
		}
		if(j > i){
			M[j][i] = d;
		}
	}
	public void set(int i, int j, int d){
		if(i >= n || j >= n){
			throw new IndexOutOfBoundsException();
		}
		if(i > j){
			M[i][j] = d;
		}
		if(j > i){
			M[j][i] = d;
		}
	}
	public NodeDist getSubNodeDist(List <String> names){
		NodeDist D = new NodeDist(names);
		for(String s: names){
			for(String t: names){
				if(s.compareTo(t) > 0){
					D.set(s,t,get(s,t));
				}
			}
		}
		return D;
	}
	
	
	
	
}
