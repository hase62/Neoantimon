package tensor;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.util.zip.DataFormatException;

import utility.Dist;
import utility.HierarchicalClustering;
import utility.HierarchicalClustering.ClusterDistFuncType;
import utility.HierarchicalClustering.DistFuncType;
import utility.MyMat;

public class ClusteredOrder3Tensor extends Order3Tensor {	
	private static final long serialVersionUID = 1L;
	protected  HierarchicalClustering order1Clustering  = null;
	protected  HierarchicalClustering order2Clustering  = null;
	protected  HierarchicalClustering order3Clustering  = null;
	protected DistFuncType order1DistType = DistFuncType.DISSIMILARITY;
	protected DistFuncType order2DistType = DistFuncType.DISSIMILARITY;
	protected DistFuncType order3DistType = DistFuncType.DISSIMILARITY;
	protected ClusterDistFuncType order1ClusteringType = ClusterDistFuncType.AVERAGE;
	protected ClusterDistFuncType order2ClusteringType = ClusterDistFuncType.AVERAGE;
	protected ClusterDistFuncType order3ClusteringType = ClusterDistFuncType.AVERAGE;
	protected boolean doOrder1Clustering = true;
	protected boolean doOrder2Clustering = true;
	protected boolean doOrder3Clustering = true;
	
	public ClusteredOrder3Tensor(int n1, int n2, int n3) {
		super(n1, n2, n3);
	}
	public ClusteredOrder3Tensor(List<String> S1, List<String> S2, List<String> S3) {
		super(S1, S2, S3);
	}
	public ClusteredOrder3Tensor(Order3Tensor T) {
		 super(T);
	}
	public ClusteredOrder3Tensor(List<MyMat> Mlist){
		super(Mlist);
	}
	public ClusteredOrder3Tensor(String infile) throws IOException, DataFormatException{
		super(infile);
	}
	//shallow copy
	public ClusteredOrder3Tensor(ClusteredOrder3Tensor T) {
		 super(T);
		 order1Clustering  = (T.order1Clustering!=null)?new HierarchicalClustering(T.order1Clustering):null;
		 order2Clustering  = (T.order2Clustering!=null)?new HierarchicalClustering(T.order2Clustering):null;
		 order3Clustering  = (T.order3Clustering!=null)?new HierarchicalClustering(T.order3Clustering):null;
		 order1DistType = T.order1DistType;
		 order2DistType = T.order2DistType;
		 order3DistType = T.order3DistType;
		 order1ClusteringType = T.order1ClusteringType;
		 order2ClusteringType = T.order2ClusteringType;
		 order3ClusteringType = T.order3ClusteringType;
		 doOrder1Clustering = T.doOrder1Clustering;
		 doOrder2Clustering = T.doOrder2Clustering;
		 doOrder3Clustering = T.doOrder3Clustering;
	}
	
	public static ClusteredOrder3Tensor readClusteredOrder3TensorFromTxet(List<String> infiles) throws IOException, DataFormatException{
		List <MyMat> Mlist = new ArrayList<MyMat>();
		List <String> name3 = new ArrayList<String>();
		Pattern p = Pattern.compile("^.*/(.*)");
		Pattern p2 = Pattern.compile("([a-zA-Z0-9]+?)[.].*");
		for(int i = 0;i < infiles.size();i++){
			Mlist.add(MyMat.readMyMatFromText(infiles.get(i)));
			 Matcher m = p.matcher(infiles.get(i));
			 if(m.matches()){
				 m = p2.matcher(m.group(1));
			 }else{
				 m = p2.matcher(infiles.get(i));
			 }
			 if(m.matches()){
				 name3.add(m.group(1));
			 }else{
				 name3.add(infiles.get(i));
			 }
		}
		ClusteredOrder3Tensor T = new ClusteredOrder3Tensor(Mlist);
		T.setName3(name3);
		return T;
	}
	
	public void setOrder1Clustering(HierarchicalClustering H){
		reorderOrder1(H.getSortedTerminalNodes());
		order1Clustering = H;
	}
	
	public void setOrder2Clustering(HierarchicalClustering H){
		reorderOrder2(H.getSortedTerminalNodes());
		order2Clustering = H;
	}
	
	public void setOrder3Clustering(HierarchicalClustering H){
		reorderOrder3(H.getSortedTerminalNodes());
		order3Clustering = H;
	}
	
	
	
	public void supressOrder1Clustering(){
		doOrder1Clustering = false;
	}
	public HierarchicalClustering getOrder1Clustering(){
		return order1Clustering;
	}
	public boolean isOrder1Clustered(){
		return (order1Clustering == null)?false:true;
	}
	public void setOrder1DistFunc(DistFuncType distFuncType){
		order1DistType  = distFuncType;
	}
	public void setOrder1ClusterDistFunc(ClusterDistFuncType clusterDistFuncType){
		 order1ClusteringType = clusterDistFuncType;
	}
	
	public void supressOrder2Clustering(){
		doOrder2Clustering = false;
	}
	public HierarchicalClustering getOrder2Clustering(){
		return order2Clustering;
	}
	public boolean isOrder2Clustered(){
		return (order2Clustering == null)?false:true;
	}
	public void setOrder2DistFunc(DistFuncType distFuncType){
		order2DistType  = distFuncType;
	}
	public void setOrder2ClusterDistFunc(ClusterDistFuncType clusterDistFuncType){
		 order2ClusteringType = clusterDistFuncType;
	}
	
	public void supressOrder3Clustering(){
		doOrder3Clustering = false;
	}
	public HierarchicalClustering getOrder3Clustering(){
		return order3Clustering;
	}
	public boolean isOrder3Clustered(){
		return (order3Clustering == null)?false:true;
	}
	public void setOrder3DistFunc(DistFuncType distFuncType){
		order3DistType  = distFuncType;
	}
	public void setOrder3ClusterDistFunc(ClusterDistFuncType clusterDistFuncType){
		 order3ClusteringType = clusterDistFuncType;
	}

	public void setClusterDistFunc(ClusterDistFuncType clusterDistFuncType){
		 order1ClusteringType = clusterDistFuncType;
		 order2ClusteringType = clusterDistFuncType;
		 order3ClusteringType = clusterDistFuncType;
	}
	
	public void performClustering(){
		if(doOrder1Clustering){
			if(order1ClusteringType.equals(ClusterDistFuncType.WARD)){
				order1Clustering = HierarchicalClustering.getHierarchicalClusteringBasedWardMethod(Order3Tensor.collapseOrder2and3(this));
			}else{
				Dist d;	
				if(order1DistType.equals(DistFuncType.DISSIMILARITY)){
					d = getDistBetweenOrder1Slices(this);
				}else{
					d = getDistBetweenOrder1Slices(this);
				}
			
				order1Clustering = new HierarchicalClustering(d);
				order1Clustering.setDistFunc(order1DistType);
				order1Clustering.setClusterDistFunc(order1ClusteringType);
			}
			order1Clustering.perform();
			reorderOrder1(order1Clustering.getSortedTerminalNodes());
		}
		if(doOrder2Clustering){
			if(order2ClusteringType.equals(ClusterDistFuncType.WARD)){
				order2Clustering = HierarchicalClustering.getHierarchicalClusteringBasedWardMethod(Order3Tensor.collapseOrder3and1(this));
			}else{
				Dist d;
				if(order2DistType.equals(DistFuncType.DISSIMILARITY)){
					d = getDistBetweenOrder2Slices(this);
				}else{
					d = getDistBetweenOrder2Slices(this);
				}
			
				order2Clustering = new HierarchicalClustering(d);
				order2Clustering.setDistFunc(order2DistType);
				order2Clustering.setClusterDistFunc(order2ClusteringType);
			}
			order2Clustering.perform();
			reorderOrder2(order2Clustering.getSortedTerminalNodes());
		}
		if(doOrder3Clustering){
			if(order3ClusteringType.equals(ClusterDistFuncType.WARD)){
				order3Clustering = HierarchicalClustering.getHierarchicalClusteringBasedWardMethod(Order3Tensor.collapseOrder1and2(this));
			}else{
			Dist d;
				if(order3DistType.equals(DistFuncType.DISSIMILARITY)){
					d = getDistBetweenOrder3Slices(this);
				}else{
					d = getDistBetweenOrder3Slices(this);
				}
			
				order3Clustering = new HierarchicalClustering(d);
				order3Clustering.setDistFunc(order3DistType);
				order3Clustering.setClusterDistFunc(order3ClusteringType);
			}
			order3Clustering.perform();
			reorderOrder3(order3Clustering.getSortedTerminalNodes());
		}
		
	}
	
	public void printAsBinary(String outfile) throws FileNotFoundException, IOException{
		ObjectOutputStream out = new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream(outfile)));
		out.writeObject(this);
		out.close();
	}
	public static ClusteredOrder3Tensor readFromBinary(String infile) throws FileNotFoundException, IOException, ClassNotFoundException{
		ObjectInputStream in = new ObjectInputStream(new BufferedInputStream(new FileInputStream(infile)));
		ClusteredOrder3Tensor T =  (ClusteredOrder3Tensor)in.readObject();
		in.close();
		return T;
	}
	
	public void copyClustering(int from, int to){
		 HierarchicalClustering tmp = null;
		if(from == 1){
			tmp = order1Clustering;
		}else if(from == 2){
			tmp = order2Clustering;
		}else if(from == 3){
			tmp = order3Clustering;
		}else{
			return;
		}
		
		if(to == 1){
			setOrder1Clustering(tmp);
		}else if(to == 2){
			setOrder2Clustering(tmp);
		}else if(to == 3){
			setOrder3Clustering(tmp);
		}else{
			return;
		}
		
	}
	
}
