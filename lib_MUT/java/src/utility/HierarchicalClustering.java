package utility;

import java.io.*;
import java.util.*;


public class HierarchicalClustering implements Serializable{
	private static final long serialVersionUID = 4183693325140498634L;
	private Dist dist;
	private MyMat M;
	private VariableDist clusterDist;
	
	
	private Set <String> terminalNodes;
	private List<String> sortedTerminalNodes;
	private Map <String, Set<String>> node2terminalNodes;
	private Map <String, String[]> node2daughters;
	private Map <String, String> node2parent;
	private Map <String, Double> node2dist;
	
	
	
	private int itr;
	private int N;
	
	public static class ClusterDistFuncType  implements Serializable {
		private static final long serialVersionUID = -1140169245215004688L;
		private final String name;
		private ClusterDistFuncType(String name){this.name = name;}
		public ClusterDistFuncType(ClusterDistFuncType b){this.name = b.name;}
		@Override
		public String toString(){return name;}
		public boolean equals(ClusterDistFuncType b){return name.equals(b.toString());};
		public static final ClusterDistFuncType SINGLE = new ClusterDistFuncType("Single");
		public static final ClusterDistFuncType COMPLETE = new ClusterDistFuncType("Complete");
		public static final ClusterDistFuncType AVERAGE = new ClusterDistFuncType("Average");
		public static final ClusterDistFuncType WARD = new ClusterDistFuncType("Ward");
	}
	public static class DistFuncType  implements Serializable {
		private static final long serialVersionUID = 5734173277281865009L;
		private final String name;
		private DistFuncType(String name){this.name = name;}
		public DistFuncType(DistFuncType b){this.name = b.name;}
		@Override
		public String toString(){return name;}
		public boolean equals(DistFuncType b){return name.equals(b.toString());};
		public static final DistFuncType DISSIMILARITY = new DistFuncType("Dissimilality");
		public static final DistFuncType SIMILARITY = new DistFuncType("Similality");
		public static final DistFuncType CORRELATION = new DistFuncType("Correlation");
		public static final DistFuncType ABSOLUTE_CORRELATION = new DistFuncType("Absolute_Correlation");
		public static final DistFuncType WARD = new DistFuncType("Ward");
	}
	
	interface DistFunc {
		double get(String a, String b);
	}
	interface ClusterDistFunc{
		double get(Set<String> a, Set<String> b);
	}
	transient DistFuncType distFuncType;
	transient DistFunc distFunc;
	transient ClusterDistFuncType clusterDistFuncType;
	transient ClusterDistFunc clusterDistFunc;
	private void setDistFunc2Dissimilarity(){
		distFuncType = DistFuncType.DISSIMILARITY;
		distFunc =new DistFunc (){
			public double get(String a, String b){
				return dist.get(a,b);
			}
		};
	}
	private void setDistFunc2Similarity(){
		distFuncType = DistFuncType.SIMILARITY;
		distFunc =new DistFunc(){
			public double get(String a, String b){
				return -dist.get(a,b);
			}
		};
	}
	private void setDistFunc2Correlation(){
		distFuncType = DistFuncType.CORRELATION;
		distFunc =new DistFunc(){
			public double get(String a, String b){
				return 1-dist.get(a,b);
			}
		};
	}
	private void setDistFunc2AbsoluteCorrelation(){
		distFuncType = DistFuncType.ABSOLUTE_CORRELATION;
		distFunc =new DistFunc(){
			public double get(String a, String b){
				return 1-Math.abs(dist.get(a,b));
			}
		};
	}
	private void setDistFunc2Ward(){
		distFuncType = DistFuncType.WARD;
		distFunc =new DistFunc(){
			public double get(String a, String b){
				//return dist.get(a,b);
				return Math.pow(dist.get(a,b),2)/2;
			}
		};
	}
	public void setDistFunc(DistFuncType distFuncType){
		if(distFuncType.equals(DistFuncType.DISSIMILARITY)){
			setDistFunc2Dissimilarity();
			return;
		}
		if(distFuncType.equals(DistFuncType.SIMILARITY)){
			setDistFunc2Similarity();
			return;
		}	
		if(distFuncType.equals(DistFuncType.CORRELATION)){
			setDistFunc2Correlation();
			return;
		}	
		if(distFuncType.equals(DistFuncType.ABSOLUTE_CORRELATION)){
			setDistFunc2AbsoluteCorrelation();
			return;
		}
	}
	private void setClusterDistFunc2Single(){
		clusterDistFuncType  = ClusterDistFuncType.SINGLE;
		clusterDistFunc = new ClusterDistFunc(){
			public double get(Set<String> a, Set<String> b){
				double min = Double.MAX_VALUE; 
				for(String e1 : a){
					for(String e2 : b){
						double d = distFunc.get(e1, e2);
						if(d < min){
							min = d;
						}
					}
				}
				return min;
			}
		};
	}
	private void setClusterDistFunc2Complete(){
		clusterDistFuncType  = ClusterDistFuncType.COMPLETE;
		clusterDistFunc = new ClusterDistFunc(){
			public double get(Set<String> a, Set<String> b){
				double max = -Double.MAX_VALUE; 
				for(String e1 : a){
					for(String e2 : b){
						double d = distFunc.get(e1, e2);
						if(d > max){
							max = d;
						}
					}
				}
				return max;
			}
		};
	}
	private void setClusterDistFunc2Average(){
		clusterDistFuncType  = ClusterDistFuncType.AVERAGE;
		clusterDistFunc = new ClusterDistFunc(){
			public double get(Set<String> a, Set<String> b){
				double ave = 0; 
				for(String e1 : a){
					for(String e2 : b){
						ave += distFunc.get(e1, e2);
					}
				}
				return ave/(a.size()*b.size());
			}
		};
	}
	
	private void setClusterDistFunc2Ward(){
		clusterDistFuncType  = ClusterDistFuncType.WARD;
		clusterDistFunc = new ClusterDistFunc(){
			public double get(Set<String> a, Set<String> b){
				List <Double> meanA =  M.getRowMeans(new ArrayList<String>(a));
				List <Double> meanB =  M.getRowMeans(new ArrayList<String>(b));
				double tmp = 0;
				for(int i = 0; i < meanA.size(); i++){
					tmp += Math.pow(meanA.get(i)-meanB.get(i),2);
				}
				//double tmp2 = 2*(a.size()*b.size())/(a.size()+b.size());	
				//return  Math.pow(tmp*tmp2, 0.5);
				return  (tmp*a.size()*b.size())/(a.size()+b.size());	
				
			}
		};
	}
	
	public void setClusterDistFunc(ClusterDistFuncType clusterDistFuncType){
		if(clusterDistFuncType.equals(ClusterDistFuncType.AVERAGE)){
			setClusterDistFunc2Average();
			return;
		}
		if(clusterDistFuncType.equals(ClusterDistFuncType.COMPLETE)){
			setClusterDistFunc2Complete();
			return;
		}
		if(clusterDistFuncType.equals(ClusterDistFuncType.SINGLE)){
			setClusterDistFunc2Single();
			return;
		}
	}

	
    public HierarchicalClustering(){};
	
    public HierarchicalClustering(HierarchicalClustering H){
    	dist = H.dist;
    	M = H.M;
    	clusterDist = H.clusterDist;
    	
    	distFuncType = H.distFuncType;
    	distFunc = H.distFunc;
    	clusterDistFuncType = H.clusterDistFuncType;
    	clusterDistFunc = H.clusterDistFunc;
    
    	terminalNodes = H.terminalNodes!=null?new HashSet<String>(H.terminalNodes):null;
    	sortedTerminalNodes = H.sortedTerminalNodes!=null?new ArrayList<String>(H.sortedTerminalNodes):null;
    	node2terminalNodes = H.node2terminalNodes!=null?new HashMap<String, Set<String>>(H.node2terminalNodes):null;
    	node2daughters = H.node2daughters!=null?new HashMap<String, String[]>(H.node2daughters):null;
    	node2parent = H.node2parent!=null?new HashMap<String, String>(H.node2parent):null;
    	node2dist = H.node2dist!=null?new HashMap<String, Double>(H.node2dist):null;
    	
    	itr = H.itr;
    	N = H.N;
    };
    
    
	public HierarchicalClustering(Dist dist){
		this.dist = dist;
		terminalNodes = new HashSet<String>(dist.getNames());
		node2terminalNodes  = new HashMap<String, Set<String>>(); 
		for(String e: terminalNodes){
			Set <String> tmp = new HashSet<String>();
			tmp.add(e);
			node2terminalNodes.put(e, tmp);
		}
		
		clusterDist = new VariableDist(new ArrayList<String>(terminalNodes));		
		itr = 0;
		N = terminalNodes.size();
		node2dist = new HashMap<String, Double>();
		node2daughters = new HashMap<String, String[]>();
		node2parent = new HashMap<String, String>();
		
		setClusterDistFunc2Average();
		setDistFunc2Dissimilarity();
	}
		
	public static HierarchicalClustering getHierarchicalClusteringBasedWardMethod(MyMat M){
		HierarchicalClustering H = new HierarchicalClustering(MyMat.getDistBetweenRows(M, 'e'));
		H.M = M;
		H.setClusterDistFunc2Ward();
		H.setDistFunc2Ward();
		return H;
	}
	
	
	private void mergePairWithMinDist(){
		itr++;
		List <String> currentCluster =  clusterDist.getNames();
 		double minDist = Double.MAX_VALUE;
 		String[] minPair  = new String[2]; 
		for(String s: currentCluster){
 			for(String t: currentCluster){
 				if(s.compareTo(t) > 0){
 					double d = clusterDist.get(s, t);
 					if(d < minDist){
 						minDist = d;
 						minPair[0] = s;
 						minPair[1] = t;
 					}
 				}
 			}
 		}
		Set <String> merged = new LinkedHashSet<String>();
		merged.addAll(node2terminalNodes.get(minPair[0]));
		merged.addAll(node2terminalNodes.get(minPair[1]));
		String nodeName = "node" + itr;
		node2dist.put(nodeName, minDist);
		node2daughters.put(nodeName, minPair);
		node2parent.put(minPair[0], nodeName);
		node2parent.put(minPair[1], nodeName);
		node2terminalNodes.put(nodeName, merged);
		clusterDist.remove(minPair[0]);
		clusterDist.remove(minPair[1]);
		Map <String, Double> dist = new HashMap<String, Double>();
		for(String e: clusterDist.getNames()){
			dist.put(e, clusterDistFunc.get(merged, node2terminalNodes.get(e)));
		}
		clusterDist.add(nodeName, dist);
	}
	public void  perform(){
		for(String s: terminalNodes){
 			for(String t: terminalNodes){
 				if(s.compareTo(t) > 0){
 					clusterDist.set(s, t, distFunc.get(s,t));	
 				}
 			}
		}		
		while(itr < N-1){
			mergePairWithMinDist();
		}
		List <Double> d = new ArrayList<Double>(node2dist.values());
		double max = MyFunc.max(d);
		double min  =  MyFunc.min(d);
		for(Map.Entry<String, Double> e:node2dist.entrySet()){
			e.setValue((e.getValue()-min)/(max-min));
		}
		sortedTerminalNodes = sortDescendants("node" + (N-1));		
	}
	
	//get k clusters based on the dendrogram 
	public List<List <String>> cutTree(int k){
		if(k < 2){
			return null;
		}
		Set <String> ids = new HashSet<String>();
		String nextNode = "node" + (N-1);
		ids.add(node2daughters.get(nextNode)[0]);
		ids.add(node2daughters.get(nextNode)[1]);
		for(int i = 1; i<k-1 ; i++){
			nextNode = "node" + (N-1-i);
			ids.remove(nextNode);
			ids.add(node2daughters.get(nextNode)[0]);
			ids.add(node2daughters.get(nextNode)[1]);
		}
		List < List<String> > clusters = new ArrayList< List<String> >();
		for(String s: ids){
			clusters.add(new ArrayList<String>(node2terminalNodes.get(s)));
		}
		return clusters;
	}
	
	
	
	
	public Map <String, String> getCutTreeMap(int k){
		List <List<String>> tmp = cutTree(k);
		Collections.sort(tmp, new MyComparator()); 
		Map <String, String> tmp2 = new LinkedHashMap<String, String>();
		for(int i =0; i < tmp.size(); i++){
			for(String s: tmp.get(i)){
				tmp2.put(s, "cluster" + (i+1));
			}
		}
		return tmp2;
	}
	
	protected class MyComparator implements Comparator { 
		public int compare(Object o1, Object o2) {
			if( ((List)o2).size() >(((List)o1).size())){
				return 1;
			}else if(((List)o2).size()  < (((List)o1).size())){
				return -1;
			}else{
				return 0;
			}
		}  
	}  
		
	public void writeDendrogramData(String outfile) throws IOException{
		PrintWriter os = new PrintWriter(new FileWriter(outfile));
		for(int i = 1; i <= N-1; i++){
			String id = "node" + i;
			os.println(id + "\t" +  node2daughters.get(id)[0] + "\t" + node2daughters.get(id)[1] + "\t" + node2dist.get(id));  
		}
		os.close();
	}
	
	private LinkedList <String> sortDescendants(String node){
		LinkedList <String> sorted;
		String first = node2daughters.get(node)[0];
		String second = node2daughters.get(node)[1];
		if(terminalNodes.contains(first)&& terminalNodes.contains(second)){
			sorted  = new LinkedList<String>();
			sorted.add(first);
			sorted.add(second);
			return sorted;
		}
		if(!terminalNodes.contains(first) && terminalNodes.contains(second)){
			sorted  = new LinkedList<String>(sortDescendants(first));
			String head  = sorted.getFirst();
			String tail = sorted.getLast();
			if(dist.get(head, second) > dist.get(tail, second)){
				sorted.addLast(second);
			}else{
				sorted.addFirst(second);
			}
			return sorted;
		}
		if(terminalNodes.contains(first) && !terminalNodes.contains(second)){
			sorted = new LinkedList<String>(sortDescendants(second));
			String head  = sorted.getFirst();
			String tail = sorted.getLast();
			if(distFunc.get(head, first) > distFunc.get(tail, first)){
				sorted.addLast(first);
			}else{
				sorted.addFirst(first);
			}
			return sorted;
		}
		sorted = new LinkedList<String>(sortDescendants(first));
		LinkedList <String> sorted2 =  new LinkedList<String>(sortDescendants(second));
		String head  = sorted.getFirst();
		String tail = sorted.getLast();
		String head2  = sorted2.getFirst();
		String tail2 = sorted2.getLast();
		double d[] = new double[4];  
		d[0] = distFunc.get(head, head2);
		d[1] = distFunc.get(head, tail2);
		d[2] = distFunc.get(tail, head2);
		d[3] = distFunc.get(tail, tail2);
		double minD = d[0];
		int minI = 0;
		for(int i = 1; i < 4; i++){
			if(d[i] <  minD){
				minD = d[i];
				minI = i;
			}
		}
		switch(minI){
		case 0:
			Collections.reverse(sorted);
			sorted.addAll(sorted2);
			return sorted;
		case 1:
			Collections.reverse(sorted);
			Collections.reverse(sorted2);
			sorted.addAll(sorted2);
			return sorted;
		case 2:
			sorted.addAll(sorted2);
			return sorted;
		case 3:
			Collections.reverse(sorted2);
			sorted.addAll(sorted2);
			return sorted;
		default:
			return null;
		}
	}
	
	public List <String> getSortedTerminalNodes(){
		return sortedTerminalNodes;
	}
	
	public List <String> getTerminalNodes(String node){
		return new ArrayList <String> (node2terminalNodes.get(node));
	}
	
	
	public List <String> getDownstreamNodes(String node){
		List <String> tmp = new ArrayList <String>();
		if(node2daughters.containsKey(node)){
			String d1 = node2daughters.get(node)[0];
			String d2 = node2daughters.get(node)[1];
			tmp.add(d1);
			tmp.add(d2);
			tmp.addAll(getDownstreamNodes(d1));
			tmp.addAll(getDownstreamNodes(d2));
			return tmp;
		}else{
			return tmp;
		}
	}
	
	
	
	private  List <String> getUpstreamNodes(List <String> node){
		Set <String> tmp = new HashSet<String>();
		boolean update = false;
		for(String s: node){
			if(node2parent.containsKey(s) ){
				String parent = node2parent.get(s);	
				if(!tmp.contains(parent) & node.contains(node2daughters.get(parent)[0]) & node.contains(node2daughters.get(parent)[1])){
					tmp.add(parent);
					update = true;
				}
				}
			}
			List <String> tmp2 = new ArrayList<String>(tmp);
		if(update){
			tmp2.addAll(getUpstreamNodes(tmp2));
		}
		return tmp2;
	}
	
	private  String getUpstreamCommonTopNodes(List<String> node){
		Set <String> tmp  = new HashSet <String>(getUpstreamNodes(node));
		String top = null;
		for(int i = N-1; i> 0; i--){
			if(tmp.contains("node" + i)){
				top  = "node" + i;
			}
		}
		return top;
	}
	
 	public Map <String, String[]> getNode2daughterMap(){
		return node2daughters;
	}
	
	public Map <String, Double> getNode2distMap(){
		return node2dist;
	}
	public Map <String, String> getNode2parentMap(){
		return node2parent;
	}
	
	public HierarchicalClustering getSubHierarchicalClustering(String node){
		HierarchicalClustering H = new HierarchicalClustering();
		H.terminalNodes = node2terminalNodes.get(node);
		H.sortedTerminalNodes = new ArrayList <String>();
		for(String n:sortedTerminalNodes){
			if(H.terminalNodes.contains(n)){
				H.sortedTerminalNodes.add(n);
			}
		}
		
		Set <String> middle = new  HashSet<String>(getDownstreamNodes(node));
		middle.removeAll(H.terminalNodes);
		Set <String> topAndMiddle = new  HashSet<String>(middle);
		topAndMiddle.add(node);
		
		
		Map <String, String> old2new = new HashMap<String,String>();
		int j = 1;
		for(int i = 1; i < terminalNodes.size(); i++){
			if(topAndMiddle.contains("node" + i)){
				old2new.put("node" + i, "node" + j);
				j++;
			}
		}
		
		
		H.node2terminalNodes = new HashMap<String, Set<String>>();
		H.node2daughters = new HashMap<String, String[]>();
		for(String s: topAndMiddle){
			String s2 = old2new.get(s); 
			H.node2terminalNodes.put(s2, node2terminalNodes.get(s));
			String[] tmp = node2daughters.get(s);
			String[] tmp2 = new String[2];
			if(H.terminalNodes.contains(tmp[0])){
				tmp2[0] = tmp[0];
			}else{
				tmp2[0] = old2new.get(tmp[0]);
			}
			if(H.terminalNodes.contains(tmp[1])){
				tmp2[1] = tmp[1];
			}else{
				tmp2[1] = old2new.get(tmp[1]);
			}
			H.node2daughters.put(s2, tmp2);
		}
		
		
		
		H.node2parent = new HashMap<String, String>();
		for(String s: terminalNodes){
			H.node2parent.put(s, old2new.get(node2parent.get(s)));
		}
		for(String s: middle){
			String s2 = old2new.get(s);
			H.node2parent.put(s2, old2new.get(node2parent.get(s)));
		}
		
		
		
		H.node2dist = new HashMap<String, Double>();
		for(String s: topAndMiddle){
			String s2 = old2new.get(s);
			H.node2dist.put(s2, node2dist.get(s));
		}
		List <Double> d = new ArrayList<Double>(H.node2dist.values());
		double max = MyFunc.max(d);
		double min  =  MyFunc.min(d);
		for(Map.Entry<String, Double> e:H.node2dist.entrySet()){
			e.setValue((e.getValue()-min)/(max-min));
		}
		
		

		/*H.dist = dist.getSubDist(H.sortedTerminalNodes);
		if(M != null){
			H.M = M.getSubMatByRow(sortedTerminalNodes);
		}
		H.clusterDist = clusterDist.getSubVariableDist(sortedTerminalNodes);*/
		return H;
	}
	
	
	public HierarchicalClustering getSubHierarchicalClusteringFromCuttingTreeMap( Map <String, String> gene2cluster){
		HierarchicalClustering H = new HierarchicalClustering(); 
		
		H.terminalNodes = new HashSet <String>(gene2cluster.values());
		Set <String> middle = new  HashSet<String>(getUpstreamNodes(new ArrayList<String>(gene2cluster.keySet())));
		Set <String> topAndMiddle = new  HashSet<String>(middle);
		middle.remove("node" + (N-1));
	
		Map <String, List <String>> cluster2gene = new HashMap<String, List<String>>();
		for(String c: gene2cluster.values()){
				cluster2gene.put(c, new ArrayList<String>());
		}
		for(String g: gene2cluster.keySet()){
			cluster2gene.get(gene2cluster.get(g)).add(g);			
		}
		
		Map <String, String> old2new = new HashMap<String,String>();
		int j = 1;
		for(int i = 1; i < terminalNodes.size(); i++){
			if(topAndMiddle.contains("node" + i)){
				old2new.put("node" + i, "node" + j);
				j++;
			}
		}
		for(String c: cluster2gene.keySet()){
			String tmp = getUpstreamCommonTopNodes(cluster2gene.get(c));
			old2new.put(tmp,c);
		}
		
		H.node2terminalNodes = new HashMap<String, Set<String>>();
		H.node2daughters = new HashMap<String, String[]>();
		for(String s: topAndMiddle){
			String s2 = old2new.get(s); 
			List<String> tmp3 = getDownstreamNodes(s);
			tmp3.retainAll(gene2cluster.keySet());
			Set<String> tmp4 = new HashSet<String>();
			for(String t: tmp3){
				tmp4.add(old2new.get(t));
			}
			H.node2terminalNodes.put(s2, tmp4);
			String[] tmp = node2daughters.get(s);
			String[] tmp2 = new String[2];
			tmp2[0] = old2new.get(tmp[0]);
			tmp2[1] = old2new.get(tmp[1]);
			H.node2daughters.put(s2, tmp2);
		}
		
		
		H.node2parent = new HashMap<String, String>();
		for(String s: gene2cluster.keySet()){
			H.node2parent.put(old2new.get(s), old2new.get(node2parent.get(s)));
		}
		for(String s: middle){
			String s2 = old2new.get(s);
			H.node2parent.put(s2, old2new.get(node2parent.get(s)));
		}
		
		

		H.node2dist = new HashMap<String, Double>();
		for(String s: topAndMiddle){
			String s2 = old2new.get(s);
			H.node2dist.put(s2, node2dist.get(s));
		}
		List <Double> d = new ArrayList<Double>(H.node2dist.values());
		double max = MyFunc.max(d);
		double min  =  MyFunc.min(d);
		for(Map.Entry<String, Double> e:H.node2dist.entrySet()){
			e.setValue((e.getValue()-min)/(max-min));
		}
		
		
		H.sortedTerminalNodes = new ArrayList<String>();
		for(String s: sortedTerminalNodes){
			H.sortedTerminalNodes.add(gene2cluster.get(s));			
		}
		H.sortedTerminalNodes = MyFunc.uniq(H.sortedTerminalNodes);
		
		return H;
	}
	
}
