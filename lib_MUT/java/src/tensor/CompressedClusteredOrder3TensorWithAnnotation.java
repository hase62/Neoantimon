package tensor;

import java.io.IOException;
import java.util.*;
import java.util.zip.DataFormatException;

import utility.HierarchicalClustering;


public class CompressedClusteredOrder3TensorWithAnnotation extends ClusteredOrder3TensorWithAnnotation {

	private Map <String, List <String>> order1cluster2gene = null;
	private Map <String, List <String>> order2cluster2gene = null;
	private Map <String, List <String>> order3cluster2gene = null;
	
	private Map <String, String> order1gene2cluster = null;
	private Map <String, String> order2gene2cluster = null;
	private Map <String, String> order3gene2cluster = null;
	
	
	public CompressedClusteredOrder3TensorWithAnnotation(ClusteredOrder3TensorWithAnnotation T) {
			super(T);
	}
	
	public CompressedClusteredOrder3TensorWithAnnotation(ClusteredOrder3Tensor T) {
		super(T);
}
	
	public CompressedClusteredOrder3TensorWithAnnotation(String string) throws IOException, DataFormatException {
		super(string);
	}

	public void compressOrder1(int k){
		order1gene2cluster = getOrder1Clustering().getCutTreeMap(k);
		order1cluster2gene = new HashMap<String, List<String>>();
		for(String c: order1gene2cluster.values()){
				order1cluster2gene.put(c, new ArrayList<String>());
		}
		for(String g: order1gene2cluster.keySet()){
			order1cluster2gene.get(order1gene2cluster.get(g)).add(g);			
		}
		order1Annotation = null;
		order1Clustering = order1Clustering.getSubHierarchicalClusteringFromCuttingTreeMap(order1gene2cluster);		
	
	
	
	
	
	}
	
	






}
