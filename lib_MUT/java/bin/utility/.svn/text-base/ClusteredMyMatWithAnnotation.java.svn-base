package utility;

import java.io.*;
import java.io.ObjectOutputStream;
import java.util.*;
import java.util.zip.DataFormatException;


public class ClusteredMyMatWithAnnotation extends ClusteredMyMat {
	private static final long serialVersionUID = 7504010036081659458L;
	private StringMat rowAnnotation = null;  // row:rowname  col:annotationType
	private StringMat colAnnotation = null; //  row:colname  col:annotationType
	
	
	public ClusteredMyMatWithAnnotation(ClusteredMyMatWithAnnotation M) {
		super(M);
		if(M.rowAnnotation != null){
			rowAnnotation = new StringMat(M.rowAnnotation);
		}
		if(M.colAnnotation != null){
			colAnnotation = new StringMat(M.colAnnotation);
		}
	}
	
	public String getRowAnnotation(String annotationType, String rowname){
		return rowAnnotation.get(rowname, annotationType);
		
	}
	public List <String> getRowAnnotation(String annotationType){
		return  rowAnnotation.getCol(annotationType);
	}
	public Map <String, String> getRowAnnotationMap(String annotationType){
		return  rowAnnotation.getColMap(annotationType);
	}
	public void setRowAnnotation(String annotationType, String rowname, String value){
		rowAnnotation.set(rowname, annotationType, value);
	}
	public String getColAnnotation(String annotationType, String colname){
		return colAnnotation.get(colname, annotationType);
		
	}
	public List <String> getColAnnotation(String annotationType){
		return  colAnnotation.getCol(annotationType);
	}
	public Map <String, String> getColAnnotationMap(String annotationType){
		return  colAnnotation.getColMap(annotationType);
	}
	public StringMat getColAnnotationMatrix(){
		return colAnnotation;
	}
	public StringMat getRowAnnotationMatrix(){
		return rowAnnotation;
	}	
	public void setColAnnotation(String annotationType, String colname, String value){
		colAnnotation.set(colname, annotationType, value);
	}
	
	public List<String> getRowAnnotationTypes(){
		return rowAnnotation.getColNames();
	}
	public List<String> getColAnnotationTypes(){
		return colAnnotation.getColNames();
	}
	public ClusteredMyMatWithAnnotation(){
		super();
	}
	public ClusteredMyMatWithAnnotation(String infile) throws IOException, DataFormatException{
		super(infile);
	}
	public ClusteredMyMatWithAnnotation(MyMat m){
		super(m);
	}
	public ClusteredMyMatWithAnnotation(List <String> row, List <String> col){
		super(row, col);
	}
	
	public void setAnnotation(StringMat annotation){
		if(!MyFunc.isect(rowname, annotation.getRowNames()).isEmpty() || !MyFunc.isect(rowname, annotation.getColNames()).isEmpty()){
			setRowAnnotation(annotation);
		}
		if(!MyFunc.isect(colname, annotation.getRowNames()).isEmpty() || !MyFunc.isect(colname, annotation.getColNames()).isEmpty()){
			setColAnnotation(annotation);
		}
		return; 
	}
	
	public void addAnnotation(StringMat annotation){
		if(!MyFunc.isect(rowname, annotation.getRowNames()).isEmpty() || !MyFunc.isect(rowname, annotation.getColNames()).isEmpty()){
			addRowAnnotation(annotation);
		}
		if(!MyFunc.isect(colname, annotation.getRowNames()).isEmpty() || !MyFunc.isect(colname, annotation.getColNames()).isEmpty()){
			addColAnnotation(annotation);
		}
		return; 
	}

	
	public void setRowAnnotation(StringMat annotation){
		if(!MyFunc.isect(rowname, annotation.getRowNames()).isEmpty()){
			rowAnnotation = new StringMat(new ArrayList<String>(rowname), annotation.getColNames());
			for(String s: rowname){
				for(String t: annotation.getColNames()){
					if(annotation.containsRowName(s)){
						rowAnnotation.set(s, t, annotation.get(s, t));
					}
				}
			}
			return;
		}
		if(!MyFunc.isect(rowname, annotation.getColNames()).isEmpty()){
			rowAnnotation = new StringMat(new ArrayList<String>(rowname), annotation.getRowNames());
			for(String s: rowname){
				for(String t: annotation.getRowNames()){
					if(annotation.containsColName(s)){
						rowAnnotation.set(s, t, annotation.get(t, s));
					}
				}
			}
			return;
		}
	}
	public void addRowAnnotation(StringMat annotation){
		if(!hasRowAnnotation()){
			setRowAnnotation(annotation);
		}else{
			if(!MyFunc.isect(rowname, annotation.getRowNames()).isEmpty()){
				StringMat tmp = new StringMat(new ArrayList<String>(rowname), annotation.getColNames());
				for(String s: rowname){
					for(String t: annotation.getColNames()){
						if(annotation.containsRowName(s)){
							tmp.set(s, t, annotation.get(s, t));
						}
					}
				}
				rowAnnotation = rowAnnotation.bindCol(tmp);
				return;
			}
			if(!MyFunc.isect(rowname, annotation.getColNames()).isEmpty()){
				StringMat tmp = new StringMat(new ArrayList<String>(rowname), annotation.getRowNames());
				for(String s: rowname){
					for(String t: annotation.getRowNames()){
						if(annotation.containsColName(s)){
							tmp.set(s, t, annotation.get(t, s));
						}
					}
				}
				rowAnnotation = rowAnnotation.bindCol(tmp);
				return;
	}
		}
	}
	public void setColAnnotation(StringMat annotation){
		if(!MyFunc.isect(colname, annotation.getRowNames()).isEmpty()){
			colAnnotation = new StringMat(new ArrayList<String>(colname), annotation.getColNames());
			for(String s: colname){
				for(String t: annotation.getColNames()){
					if(annotation.containsRowName(s)){
						colAnnotation.set(s, t, annotation.get(s, t));
					}
				}
			}
			return;
		}
		if(!MyFunc.isect(colname, annotation.getColNames()).isEmpty()){
			colAnnotation = new StringMat(new ArrayList<String>(colname), annotation.getRowNames());
			for(String s: colname){
				for(String t: annotation.getRowNames()){
					if(annotation.containsColName(s)){
						colAnnotation.set(s, t, annotation.get(t, s));
					}
				}
			}
			return;
		}
	}
	public void addColAnnotation(StringMat annotation){
		if(!hasColAnnotation()){
			setColAnnotation(annotation);
		}else{
			if(!MyFunc.isect(colname, annotation.getRowNames()).isEmpty()){
				StringMat tmp = new StringMat(new ArrayList<String>(colname), annotation.getColNames());
				for(String s: colname){
					for(String t: annotation.getColNames()){
						if(annotation.containsRowName(s)){
							tmp.set(s, t, annotation.get(s, t));
						}
					}
				}
				colAnnotation = colAnnotation.bindCol(tmp);
				return;
			}
			if(!MyFunc.isect(colname, annotation.getColNames()).isEmpty()){
				StringMat tmp  = new StringMat(new ArrayList<String>(colname), annotation.getRowNames());
				for(String s: colname){
					for(String t: annotation.getRowNames()){
						if(annotation.containsColName(s)){
							tmp.set(s, t, annotation.get(t, s));
						}
					}
				}
				colAnnotation = colAnnotation.bindCol(tmp);
				return;
			}
		}
	}
	
	public boolean hasRowAnnotation(){
		return (rowAnnotation==null)?false:true;
	}
	public boolean hasColAnnotation(){
		return (colAnnotation==null)?false:true;
	}
	@SuppressWarnings("unchecked")
	public void reorderRows(List row){
		super.reorderRows(row);
		if(hasRowAnnotation()){
			rowAnnotation.reorderRows(rowname);
		}
	}
	@SuppressWarnings("unchecked")
	public void reorderCols(List col){
		super.reorderCols(col);
		if(hasColAnnotation()){
			colAnnotation.reorderRows(colname);
		}
	}
	
	public void sortColsByValue(String rowname){
		if(hasColAnnotation() && colAnnotation.containsColName(rowname)){
			Map <String,Double> tmp = MyFunc.StirngStringMap2StringDoubleMap(colAnnotation.getColMap(rowname));
			List <String> tmp2 = MyFunc.sortKeysByAscendingOrderOfValues(tmp);
			reorderCols(tmp2);
		}else{
			Map <String,Double> tmp = getRowMap(rowname);
			List <String> tmp2 = MyFunc.sortKeysByAscendingOrderOfValues(tmp);
			reorderCols(tmp2);
		}
	}
	public void sortRowsByValue(String colname){
		if(hasRowAnnotation() && rowAnnotation.containsColName(colname)){
			Map <String,Double> tmp = MyFunc.StirngStringMap2StringDoubleMap(rowAnnotation.getColMap(colname));
			List <String> tmp2 = MyFunc.sortKeysByAscendingOrderOfValues(tmp);
			reorderCols(tmp2);
		}else{
			Map <String,Double> tmp = getColMap(colname);
			List <String> tmp2 = MyFunc.sortKeysByAscendingOrderOfValues(tmp);
			reorderRows(tmp2);
		}
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
	public static ClusteredMyMatWithAnnotation  readFromBinary(String infile) throws FileNotFoundException, IOException, ClassNotFoundException{
		ObjectInputStream in = new ObjectInputStream(new BufferedInputStream(new FileInputStream(infile)));
		ClusteredMyMatWithAnnotation M =  (ClusteredMyMatWithAnnotation)in.readObject();
		in.close();
		return M;
	}
	
	
	
}
