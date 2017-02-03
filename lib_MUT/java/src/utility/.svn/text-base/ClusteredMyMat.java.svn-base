package utility;

import java.io.*;
import java.io.IOException;
import java.util.List;
import java.util.zip.DataFormatException;

import org.apache.commons.cli.Options;

import tensor.ClusteredOrder3Tensor;
import utility.HierarchicalClustering.*;

public class ClusteredMyMat extends MyMat {
	private static final long serialVersionUID = 7178378743611753498L;
	protected  HierarchicalClustering rowClustering  = null;
	protected  HierarchicalClustering colClustering = null;
	private DistFuncType rowDistType = DistFuncType.DISSIMILARITY;
	private DistFuncType colDistType = DistFuncType.DISSIMILARITY;
	private ClusterDistFuncType rowClusteringType = ClusterDistFuncType.AVERAGE;
	private ClusterDistFuncType colClusteringType = ClusterDistFuncType.AVERAGE;
	private boolean doRowClustering = true;
	private boolean doColClustering = true;
	public ClusteredMyMat(){
		super();
	}
	public ClusteredMyMat(String infile) throws IOException, DataFormatException{
		super(infile);
	}
	public ClusteredMyMat(MyMat m){
		super(m);
	}
	public ClusteredMyMat(List <String> row, List <String> col){
		super(row, col);
	}
	public void supressRowClustering(){
		doRowClustering = false;
	}
	public void supressColClustering(){
		doColClustering = false;
	}
	public HierarchicalClustering getRowClustering(){
		return rowClustering;
	}
	public HierarchicalClustering getColClustering(){
		return colClustering;
		
	}
	
	
	public ClusteredMyMat(ClusteredMyMat M) {
		 super(M);
		 rowClustering  = (M.rowClustering==null)?null:new HierarchicalClustering(M.rowClustering);
		 colClustering  = (M.colClustering==null)?null:new HierarchicalClustering(M.colClustering);
		
		 rowDistType = M.rowDistType;
		 colDistType = M.colDistType;
		 rowClusteringType = M.rowClusteringType;
		 colClusteringType = M.colClusteringType;
		 doRowClustering = M.doRowClustering;
		 doColClustering = M.doColClustering;
		
	}
	public boolean isRowClustered(){
		return (rowClustering == null)?false:true;
	}
	public boolean isColClustered(){
		return (colClustering == null)?false:true;
	}
	
	public void setRowDistFunc(DistFuncType distFuncType){
		rowDistType  = distFuncType;
	}
	public void setColDistFunc(DistFuncType distFuncType){
		colDistType = distFuncType;
	}
	public void setRowClusterDistFunc(ClusterDistFuncType clusterDistFuncType){
		 rowClusteringType = clusterDistFuncType;
	}
	public void setColClusterDistFunc(ClusterDistFuncType clusterDistFuncType){
		colClusteringType = clusterDistFuncType;
	}
	
	public void performClustering(){
		if(doRowClustering){
			if(rowClusteringType.equals(ClusterDistFuncType.WARD)){
				rowClustering = HierarchicalClustering.getHierarchicalClusteringBasedWardMethod(this);
			}else{
				Dist d;
				if(rowDistType.equals(DistFuncType.DISSIMILARITY)){
					d = getDistBetweenRows(this, 'e');
				}else{
					d = getDistBetweenRows(this, 'c');
				}
			
				rowClustering = new HierarchicalClustering(d);
				rowClustering.setDistFunc(rowDistType);
				rowClustering.setClusterDistFunc(rowClusteringType);
			}
			rowClustering.perform();
			reorderRows(rowClustering.getSortedTerminalNodes());
		}
		if(doColClustering){
			if(colClusteringType.equals(ClusterDistFuncType.WARD)){
				MyMat tmp = new MyMat(this);
				tmp.transpose();
				colClustering = HierarchicalClustering.getHierarchicalClusteringBasedWardMethod(tmp);
			}else{
				Dist d;	
				if(colDistType.equals(DistFuncType.DISSIMILARITY)){
					d = getDistBetweenCols(this, 'e');
				}else{
					d = getDistBetweenCols(this, 'c');
				}
				colClustering = new HierarchicalClustering(d);
				colClustering.setDistFunc(colDistType);
				colClustering.setClusterDistFunc(colClusteringType);
			}
			colClustering.perform();
			reorderCols(colClustering.getSortedTerminalNodes());
		}
	}

	public ClusteredMyMatViewer getViewer(){
		return new ClusteredMyMatViewer(this);
	}
	
	

	public void setRowClustering(HierarchicalClustering H){
		reorderRows(H.getSortedTerminalNodes());
		rowClustering = H;
	}
	
	public void setColClustering(HierarchicalClustering H){
		reorderCols(H.getSortedTerminalNodes());
		colClustering = H;
	}
	
	public void printAsBinary(String outfile) throws FileNotFoundException, IOException{
		ObjectOutputStream out = new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream(outfile)));
		out.writeObject(this);
		out.close();
	}
	public static ClusteredMyMat  readFromBinary(String infile) throws FileNotFoundException, IOException, ClassNotFoundException{
		ObjectInputStream in = new ObjectInputStream(new BufferedInputStream(new FileInputStream(infile)));
		ClusteredMyMat M =  (ClusteredMyMat)in.readObject();
		in.close();
		return M;
	}
	
	
}
