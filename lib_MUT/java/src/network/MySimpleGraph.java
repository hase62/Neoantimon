package network;

import java.util.*;
import java.util.zip.DataFormatException;
import java.io.*;

import org.jgrapht.*;
import org.jgrapht.alg.BellmanFordShortestPath;
import org.jgrapht.alg.NeighborIndex;
import org.jgrapht.graph.*;

import utility.*;

public class MySimpleGraph  extends  SimpleGraph<String, DefaultEdge> {
	private static final long serialVersionUID = -3065542008831673647L;
	NodeDist D = null;
	NeighborIndex<String, DefaultEdge> NI  = new NeighborIndex<String, DefaultEdge>(this);
	
	
	public MySimpleGraph(MySimpleGraph G){
		super(DefaultEdge.class);
		for(String n: G.vertexSet()){
			addVertex(n);
		}
		for(DefaultEdge e: G.edgeSet()){
			addEdge(G.getEdgeSource(e), G.getEdgeTarget(e));
		}
		D = new NodeDist(D);
	}
	
	public MySimpleGraph() {
		super(DefaultEdge.class);
	}
	
	public MySimpleGraph(List <String> node){
		super(DefaultEdge.class);
		for(String n: node){
			addVertex(n);
		}
	}
	
	
	public static MySimpleGraph getInstanceFromSifFile(String infile) throws IOException, DataFormatException{
		MySimpleGraph  g = new MySimpleGraph();
		BufferedReader inputStream = new BufferedReader(new FileReader(infile));
		String line;
		Set <String>  nameSet = new HashSet<String>();
		String[] tmp = new String[3];
 		while((line = inputStream.readLine()) != null){
			if(line.charAt(0) == '#'){
				continue;				
			}
			tmp = line.split("\t");
			if(!nameSet.contains(tmp[0])){
				g.addVertex(tmp[0]);
				nameSet.add(tmp[0]);
			}
			if(nameSet.contains(tmp[2])){
				g.addVertex(tmp[2]);
				nameSet.add(tmp[2]);
			}
			g.addEdge(tmp[0], tmp[2]);
		}
 		inputStream.close();
 		return g;
	}
	
	
	public static MySimpleGraph getInstanceFromNodePairList(List <String[]> nodePairs){
		MySimpleGraph g = new MySimpleGraph();
		for(String[] s: nodePairs){
			if(!g.containsVertex(s[0])){
				g.addVertex(s[0]);
			}
			if(!g.containsVertex(s[1])){
				g.addVertex(s[1]);
			}
		}
		for(String[] s: nodePairs){
			g.addEdge(s[0],s[1]);
		}
		return g;
	}
	
	public NodeDist calculateNodeDistanceByBellmanFordShortestPath(){
		NodeDist D =  new NodeDist(new ArrayList <String>(this.vertexSet()));
		int i = 0,  m  = (int) (0.5*D.size() *(D.size()-1));
		for(String s: this.vertexSet()){
			BellmanFordShortestPath<String, DefaultEdge>  pathFinder = new BellmanFordShortestPath<String, DefaultEdge>(this, s);
			for(String t: this.vertexSet()){
				if(s.compareTo(t) > 0){
					i++;
					try{
						int d = pathFinder.getPathEdgeList(t).size();
						D.set(s, t, d);
						System.err.println(i + "/" + m + "\t" + s + "\t" + t + "\t" + d);
					}
					catch(NullPointerException e){
						System.err.println(i + "/" + m + "\t" + s + "\t" + t + "\t" + "no path!");
					}
				}	
			}
		}
		return D;
	}
	public void setNodeDistance(){
		D = calculateNodeDistanceByBellmanFordShortestPath();	
	}
	public void setNodeDistance(String infile) throws IOException, DataFormatException{
		D = new NodeDist(infile);
	}
	public int getDistanceBetweenTowNodes(String node1, String node2){
		if(D != null){
			return D.get(node1, node2);
		}else{
			try{
				return BellmanFordShortestPath.findPathBetween(this, node1, node2).size();
			}
			catch(NullPointerException e){
				return Integer.MAX_VALUE;
			}
		}
	}
	public List<DefaultEdge> getPath(String node1, String node2){
		return BellmanFordShortestPath.findPathBetween(this, node1, node2);
	}
	
	public MySimpleGraph getSubGraph(List <String> nodes){
		MySimpleGraph g = new MySimpleGraph();
		for(String s: nodes){
			g.addVertex(s);
		}
		for(String s: nodes){
			for(String t: nodes){
				if(s.compareTo(t) > 0 && this.containsEdge(s, t)){
					g.addEdge(s, t);
				}
			}
		}
		return g;
	}
	
	public boolean hasPath(String node1, String node2, int maxPathLength){
		if(D != null){
			return D.get(node1, node2) > maxPathLength?false:true;
		}else{
			BellmanFordShortestPath<String, DefaultEdge>  pathFinder = new BellmanFordShortestPath<String, DefaultEdge>(this, node1, maxPathLength);
			if(pathFinder.getPathEdgeList(node2)== null){
				return false;
			}else{
				return true;
			}
			
		}	
	}
	
	public List<String> getNeighbors(String node){
		return NI.neighborListOf(node);		
	}
	public MySimpleGraph  getSubLinkContainingTargetNodesAndNeighbors(List <String> targetNodes, int maxPathLength){
		 List <String> neighbor = new ArrayList<String>(targetNodes);
		  List <String> member = new ArrayList<String>(targetNodes);
		  List < String[] >linkList = new ArrayList<String[]>();
		  int k =0;
		  int i,j;
		  while(k <maxPathLength){
		    List <String> tmp =new ArrayList<String>();
		    for(i=0;i<neighbor.size();i++){
		      List <String> tmp2 = getNeighbors(neighbor.get(i));
		      for(j=0;j<tmp2.size();j++){
		    	  member.add(tmp2.get(j));
		    	  tmp.add(tmp2.get(j));
		    	  String[]  tmp3  = new String[2];
		    	  tmp3[0] = neighbor.get(i);
		    	  tmp3[1] = tmp2.get(j);
		    	  linkList.add(tmp3);
		      	}
		    }
		    neighbor = tmp;
		    k++;
		  }
		  
		  member = MyFunc.uniq(member);
		  MySimpleGraph  subG = new  MySimpleGraph(member);
		  for(i=0;i<linkList.size();i++){
		    subG.addEdge((linkList.get(i))[0],(linkList.get(i))[1]);
		  }
		  return subG;
	 }
	
	public  MySimpleGraph getSubLinkContainingNodesithMinDegreeOfLinks(int minDegree){
		 List <String> member = new ArrayList<String>(vertexSet());
		 minDegree = Math.max(minDegree, 1);
		 MySimpleGraph subG;
		  while(true){
			 subG = getSubGraph(member);
		    List <String> member2 = new ArrayList<String>();;
		    for(int i=0;i<member.size();i++){
		      if(subG.getNeighbors(member.get(i)).size() >= minDegree){
		    	  member2.add(member.get(i));
		      }
		    }
		    if(member.size() == member2.size()){
		      break;
		    }
		    member = member2;
		  }
		  return subG;
	 }
	

}
