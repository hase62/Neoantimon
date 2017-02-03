package tensor;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.util.zip.DataFormatException;
import utility.*;

public class ClusteredOrder3TensorWithAnnotation extends ClusteredOrder3Tensor{
	
	protected StringMat order1Annotation = null;  // row:name1  col:annotationType
	protected StringMat order2Annotation = null; //  row:name2  col:annotationType
	protected StringMat order3Annotation = null; //  row:name3  col:annotationType
	
	
	private static final long serialVersionUID = -3533225747002429577L;
	public ClusteredOrder3TensorWithAnnotation(int n1, int n2, int n3) {
		super(n1, n2, n3);
	}
	public ClusteredOrder3TensorWithAnnotation(List<String> S1, List<String> S2, List<String> S3) {
		super(S1, S2, S3);
	}
	public ClusteredOrder3TensorWithAnnotation(Order3Tensor T) {
		 super(T);
	}
	public ClusteredOrder3TensorWithAnnotation(List<MyMat> Mlist){
		super(Mlist);
	}
	public ClusteredOrder3TensorWithAnnotation(String infile) throws IOException, DataFormatException{
		super(infile);
	}
    
	public ClusteredOrder3TensorWithAnnotation(ClusteredOrder3Tensor T) {
		super(T);
	}
	
	public static ClusteredOrder3TensorWithAnnotation getBigClusteredOrder3TensorWithAnnotation(String infile) throws IOException, DataFormatException{
		BufferedReader inputStream = new BufferedReader(new FileReader(infile));
		String line;
		int  i,j,k;
		
		List <String> name1 = new ArrayList<String>();
		List <String> name2 = new ArrayList<String>();
		List <String> name3 = new ArrayList<String>();
		
		Pattern p = Pattern.compile(">(.+)");
		
		while((line = inputStream.readLine()) != null){
			if(line.charAt(0) == '#'){
				continue;				
			}
			if(line.charAt(0) == '>'){
				List<String> str = Arrays.asList((line.split("\t")));
				Matcher m = p.matcher(str.get(0));
				if(!m.matches()){
					System.err.println(line);
					throw new DataFormatException("Order3Tensor: file format is wrong!");
				}
				name3.add(m.group(1));
						
				if(name3.size() == 1){
					for(i = 1; i < str.size(); i++){
						name2.add(str.get(i));	
					}
				}
			}else{
				List<String>  str = Arrays.asList(line.split("\t"));
				if(str.size() != name2.size()+1){
					throw new DataFormatException("Order3Tensor: file format is wrong!");
				}
				if(name3.size() == 1){
					name1.add(str.get(0));
				}
			}
		}
		inputStream = new BufferedReader(new FileReader(infile));
		
		ClusteredOrder3TensorWithAnnotation T = new ClusteredOrder3TensorWithAnnotation(name1, name2, name3);
		
		k=-1;
		i=0;
		
		
		while((line = inputStream.readLine()) != null){
			if(line.charAt(0) == '#'){
				continue;				
			}
			if(line.charAt(0) == '>'){
				k++;
				i=0;
			}else{
				List<String>  str = Arrays.asList(line.split("\t"));
				if(str.size() != name2.size()+1){
					throw new DataFormatException("Order3Tensor: file format is wrong!");
				}
				for(j=1; j<str.size(); j++){
					T.set(i, j-1, k,  Double.valueOf(str.get(j)));		
				 }			
				i++;
			}
		}
		
		return T;
	}
	
	
	public ClusteredOrder3TensorWithAnnotation(ClusteredOrder3TensorWithAnnotation T) {
		super(T);
		if(T.order1Annotation != null){
			order1Annotation = new StringMat(T.order1Annotation);
		}
		if(T.order2Annotation != null){
			order2Annotation = new StringMat(T.order2Annotation);
		}
		if(T.order3Annotation != null){
			order3Annotation = new StringMat(T.order3Annotation);
		}
	}

	public String getOrder1Annotation(String annotationType, String name1){
		return order1Annotation.get(name1, annotationType);
		
	}
	public List <String> getOrder1Annotation(String annotationType){
		return  order1Annotation.getCol(annotationType);
	}
	public Map <String, String> getOrder1AnnotationMap(String annotationType){
		return  order1Annotation.getColMap(annotationType);
	}
	public void setOrder1Annotation(String annotationType, String name1, String value){
		order1Annotation.set(name1, annotationType, value);
	}
	public StringMat getOrder1AnnotationMatrix(){
		return order1Annotation;
	}	
	public List<String> getOrder1AnnotationTypes(){
		return order1Annotation.getColNames();
	}
	
	public void setAnnotation(StringMat annotation){
		if(!MyFunc.isect(name1, annotation.getRowNames()).isEmpty() || !MyFunc.isect(name1, annotation.getColNames()).isEmpty()){
			setOrder1Annotation(annotation);
		}
		if(!MyFunc.isect(name2, annotation.getRowNames()).isEmpty() || !MyFunc.isect(name2, annotation.getColNames()).isEmpty()){
			setOrder2Annotation(annotation);
		}
		if(!MyFunc.isect(name3, annotation.getRowNames()).isEmpty() || !MyFunc.isect(name3, annotation.getColNames()).isEmpty()){
			setOrder3Annotation(annotation);
		}
		return; 
	}
	
	public void addAnnotation(StringMat annotation){
		if(!MyFunc.isect(name1, annotation.getRowNames()).isEmpty() || !MyFunc.isect(name1, annotation.getColNames()).isEmpty()){
			addOrder1Annotation(annotation);
		}
		if(!MyFunc.isect(name2, annotation.getRowNames()).isEmpty() || !MyFunc.isect(name2, annotation.getColNames()).isEmpty()){
			addOrder2Annotation(annotation);
		}
		if(!MyFunc.isect(name3, annotation.getRowNames()).isEmpty() || !MyFunc.isect(name3, annotation.getColNames()).isEmpty()){
			addOrder3Annotation(annotation);
		}
		return; 
	}
	
	public void setOrder1Annotation(StringMat annotation){
		if(!MyFunc.isect(name1, annotation.getRowNames()).isEmpty()){
			order1Annotation = new StringMat(new ArrayList<String>(name1), annotation.getColNames());
			for(String s: name1){
				for(String t: annotation.getColNames()){
					if(annotation.containsRowName(s)){
						order1Annotation.set(s, t, annotation.get(s, t));
					}
				}
			}
			return;
		}
		if(!MyFunc.isect(name1, annotation.getColNames()).isEmpty()){
			order1Annotation = new StringMat(new ArrayList<String>(name1), annotation.getRowNames());
			for(String s: name1){
				for(String t: annotation.getRowNames()){
					if(annotation.containsColName(s)){
						order1Annotation.set(s, t, annotation.get(t, s));
					}
				}
			}
			return;
		}
	}
	
	public void addOrder1Annotation(StringMat annotation){
		if(!hasOrder1Annotation()){
			setOrder1Annotation(annotation);
		}else{
			if(!MyFunc.isect(name1, annotation.getRowNames()).isEmpty()){
				StringMat tmp = new StringMat(new ArrayList<String>(name1), annotation.getColNames());
				for(String s: name1){
					for(String t: annotation.getColNames()){
						if(annotation.containsRowName(s)){
							tmp.set(s, t, annotation.get(s, t));
						}
					}
				}
				order1Annotation = order1Annotation.bindCol(tmp);
				return;
			}
			if(!MyFunc.isect(name1, annotation.getColNames()).isEmpty()){
				StringMat tmp = new StringMat(new ArrayList<String>(name1), annotation.getRowNames());
				for(String s: name1){
					for(String t: annotation.getRowNames()){
						if(annotation.containsColName(s)){
							tmp.set(s, t, annotation.get(t, s));
						}
					}
				}
				order1Annotation = order1Annotation.bindCol(tmp);
				return;
			}
		}	
	}
			
	public boolean hasOrder1Annotation(){
		return (order1Annotation==null)?false:true;
	}

	@SuppressWarnings("unchecked")
	public void reorderOrder1(List name1){
		super.reorderOrder1(name1);
		if(hasOrder1Annotation()){
			order1Annotation.reorderRows(name1);
		}
	}

	

	public String getOrder2Annotation(String annotationType, String name2){
		return order2Annotation.get(name2, annotationType);
		
	}
	public List <String> getOrder2Annotation(String annotationType){
		return  order2Annotation.getCol(annotationType);
	}
	public Map <String, String> getOrder2AnnotationMap(String annotationType){
		return  order2Annotation.getColMap(annotationType);
	}
	public void setOrder2Annotation(String annotationType, String name2, String value){
		order2Annotation.set(name2, annotationType, value);
	}
	public StringMat getOrder2AnnotationMatrix(){
		return order2Annotation;
	}	
	public List<String> getOrder2AnnotationTypes(){
		return order2Annotation.getColNames();
	}
	
	
	public void setOrder2Annotation(StringMat annotation){
		if(!MyFunc.isect(name2, annotation.getRowNames()).isEmpty()){
			order2Annotation = new StringMat(new ArrayList<String>(name2), annotation.getColNames());
			for(String s: name2){
				for(String t: annotation.getColNames()){
					if(annotation.containsRowName(s)){
						order2Annotation.set(s, t, annotation.get(s, t));
					}
				}
			}
			return;
		}
		if(!MyFunc.isect(name2, annotation.getColNames()).isEmpty()){
			order2Annotation = new StringMat(new ArrayList<String>(name2), annotation.getRowNames());
			for(String s: name2){
				for(String t: annotation.getRowNames()){
					if(annotation.containsColName(s)){
						order2Annotation.set(s, t, annotation.get(t, s));
					}
				}
			}
			return;
		}
	}
	
	public void addOrder2Annotation(StringMat annotation){
		if(!hasOrder2Annotation()){
			setOrder2Annotation(annotation);
		}else{
			if(!MyFunc.isect(name2, annotation.getRowNames()).isEmpty()){
				StringMat tmp = new StringMat(new ArrayList<String>(name2), annotation.getColNames());
				for(String s: name2){
					for(String t: annotation.getColNames()){
						if(annotation.containsRowName(s)){
							tmp.set(s, t, annotation.get(s, t));
						}
					}
				}
				order2Annotation = order2Annotation.bindCol(tmp);
				return;
			}
			if(!MyFunc.isect(name2, annotation.getColNames()).isEmpty()){
				StringMat tmp = new StringMat(new ArrayList<String>(name2), annotation.getRowNames());
				for(String s: name2){
					for(String t: annotation.getRowNames()){
						if(annotation.containsColName(s)){
							tmp.set(s, t, annotation.get(t, s));
						}
					}
				}
				order2Annotation = order2Annotation.bindCol(tmp);
				return;
			}
		}	
	}
			
	public boolean hasOrder2Annotation(){
		return (order2Annotation==null)?false:true;
	}

	@SuppressWarnings("unchecked")
	public void reorderOrder2(List name2){
		super.reorderOrder2(name2);
		if(hasOrder2Annotation()){
			order2Annotation.reorderRows(name2);
		}
	}

	

	public String getOrder3Annotation(String annotationType, String name3){
		return order3Annotation.get(name3, annotationType);
		
	}
	public List <String> getOrder3Annotation(String annotationType){
		return  order3Annotation.getCol(annotationType);
	}
	public Map <String, String> getOrder3AnnotationMap(String annotationType){
		return  order3Annotation.getColMap(annotationType);
	}
	public void setOrder3Annotation(String annotationType, String name3, String value){
		order3Annotation.set(name3, annotationType, value);
	}
	public StringMat getOrder3AnnotationMatrix(){
		return order3Annotation;
	}	
	public List<String> getOrder3AnnotationTypes(){
		return order3Annotation.getColNames();
	}
	
	
	public void setOrder3Annotation(StringMat annotation){
		if(!MyFunc.isect(name3, annotation.getRowNames()).isEmpty()){
			order3Annotation = new StringMat(new ArrayList<String>(name3), annotation.getColNames());
			for(String s: name3){
				for(String t: annotation.getColNames()){
					if(annotation.containsRowName(s)){
						order3Annotation.set(s, t, annotation.get(s, t));
					}
				}
			}
			return;
		}
		if(!MyFunc.isect(name3, annotation.getColNames()).isEmpty()){
			order3Annotation = new StringMat(new ArrayList<String>(name3), annotation.getRowNames());
			for(String s: name3){
				for(String t: annotation.getRowNames()){
					if(annotation.containsColName(s)){
						order3Annotation.set(s, t, annotation.get(t, s));
					}
				}
			}
			return;
		}
	}
	
	public void addOrder3Annotation(StringMat annotation){
		if(!hasOrder3Annotation()){
			setOrder3Annotation(annotation);
		}else{
			if(!MyFunc.isect(name3, annotation.getRowNames()).isEmpty()){
				StringMat tmp = new StringMat(new ArrayList<String>(name3), annotation.getColNames());
				for(String s: name3){
					for(String t: annotation.getColNames()){
						if(annotation.containsRowName(s)){
							tmp.set(s, t, annotation.get(s, t));
						}
					}
				}
				order3Annotation = order3Annotation.bindCol(tmp);
				return;
			}
			if(!MyFunc.isect(name3, annotation.getColNames()).isEmpty()){
				StringMat tmp = new StringMat(new ArrayList<String>(name3), annotation.getRowNames());
				for(String s: name3){
					for(String t: annotation.getRowNames()){
						if(annotation.containsColName(s)){
							tmp.set(s, t, annotation.get(t, s));
						}
					}
				}
				order3Annotation = order3Annotation.bindCol(tmp);
				return;
			}
		}	
	}
			
	public boolean hasOrder3Annotation(){
		return (order3Annotation==null)?false:true;
	}

	@SuppressWarnings("unchecked")
	public void reorderOrder3(List name3){
		super.reorderOrder3(name3);
		if(hasOrder3Annotation()){
			order3Annotation.reorderRows(name3);
		}
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
	
	public ClusteredOrder3TensorViewer getViewer(){
		return new  ClusteredOrder3TensorViewer2(this);
	}
	
	public void printAsBinary(String outfile) throws FileNotFoundException, IOException{
		ObjectOutputStream out = new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream(outfile)));
		out.writeObject(this);
		out.close();
	}
	
	public static ClusteredOrder3TensorWithAnnotation readFromBinary(String infile) throws FileNotFoundException, IOException, ClassNotFoundException{
		ObjectInputStream in = new ObjectInputStream(new BufferedInputStream(new FileInputStream(infile)));
		ClusteredOrder3TensorWithAnnotation M =  (ClusteredOrder3TensorWithAnnotation)in.readObject();
		in.close();
		return M;
	}

}
