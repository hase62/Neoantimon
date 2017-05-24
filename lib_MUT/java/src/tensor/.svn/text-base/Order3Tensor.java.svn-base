package tensor;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.DataFormatException;

import utility.Dist;
import utility.MyFunc;
import utility.MyMat;

public class Order3Tensor implements Serializable{
	private static final long serialVersionUID = 1L;
	protected double M[][][];
	protected int n1;
	protected int n2;
	protected int n3;
	protected List <String> name1;
	protected List <String> name2;
	protected List <String> name3;
	protected Map<String, Integer> name12index;
	protected Map<String, Integer> name22index;
	protected Map<String, Integer> name32index;
	
	public Order3Tensor(int n1, int n2, int n3){
		this.n1 = n1;
		this.n2 = n2;
		this.n3 = n3;
		int i, j, k; 
		List <String> name1 = new ArrayList<String>();
		List <String>  name2 = new ArrayList<String>();
		List <String>  name3 = new ArrayList<String>();
		for(i=1;i<=n1;i++){
			name1.add("1_"+i);
		}
		for(j=1;j<=n2;j++){
			name2.add("2_"+j);
		}
		for(k=1;k<=n3;k++){
			name3.add("3_"+k);
		}
		setName1(name1);
		setName2(name2);
		setName3(name3);
		M = new double[n1][n2][n3];
		for(i=0;i<n1;i++){
			for(j=0;j<n2;j++){
				for(k=0;k<n3;k++){
					M[i][j][k] = 0;	
				}
			}
		}
	}
	
	
	public Order3Tensor(List<String> name1,List<String> name2, List<String> name3){
		n1 =name1.size(); 
		n2 =name2.size(); 
		n3 =name3.size(); 
		setName1(name1);
		setName2(name2);
		setName3(name3);
		int i, j, k; 
		M = new double[n1][n2][n3];
		for(i=0;i<n1;i++){
			for(j=0;j<n2;j++){
				for(k=0;k<n3;k++){
					M[i][j][k] = 0;	
				}
			}
		}
	}
	
	public Order3Tensor(Order3Tensor T) {
		n1 =T.name1.size(); 
		n2 =T.name2.size(); 
		n3 =T.name3.size(); 
		setName1(T.name1);
		setName2(T.name2);
		setName3(T.name3);
		int i, j, k;
		M = new double[n1][n2][n3];
		for(i=0;i<n1;i++){
			for(j=0;j<n2;j++){
				for(k=0;k<n3;k++){
					M[i][j][k] = T.M[i][j][k];	
				}
			}
		}
	}
	
	
	
	public Order3Tensor(List<MyMat> Mlist){
		Set<String> commonRows = new TreeSet<String>(Mlist.get(0).getRowNames());
		for(int i = 1; i < Mlist.size(); i++){
			Set<String> tmp = new HashSet<String>(Mlist.get(i).getRowNames());
			commonRows.retainAll(tmp);
		}
		Set<String> commonCols = new TreeSet<String>(Mlist.get(0).getColNames());
		for(int i = 1; i < Mlist.size(); i++){
			Set<String> tmp = new HashSet<String>(Mlist.get(i).getColNames());
			commonCols.retainAll(tmp);
		}
		
		List <String> name1 = new ArrayList<String>(commonRows);
		List <String> name2 = new ArrayList<String>(commonCols);
		List <String> name3 = new ArrayList<String>();
		for(int k=1;k<=Mlist.size();k++){
			name3.add("3_"+k);
		}
		n1 =name1.size(); 
		n2 =name2.size(); 
		n3 =name3.size(); 
		setName1(name1);
		setName2(name2);
		setName3(name3);
		int i, j, k; 
		M = new double[n1][n2][n3];
		for(i=0;i<n1;i++){
			for(j=0;j<n2;j++){
				for(k=0;k<n3;k++){
					M[i][j][k] = Mlist.get(k).get(name1.get(i), name2.get(j));	
				}
			}
		}
	}
	
	public static Order3Tensor readOrder3TensorFromTxet(List<String> infiles) throws IOException, DataFormatException{
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
		Order3Tensor T = new Order3Tensor(Mlist);
		T.setName3(name3);
		return T;
	}
	
	public Order3Tensor(String infile) throws IOException, DataFormatException{
		BufferedReader inputStream = new BufferedReader(new FileReader(infile));
		List <List<String>> Lines = new ArrayList<List <String>>();
		String line;
		int i = -1;
		while((line = inputStream.readLine()) != null){
			if(line.charAt(0) == '#'){
				continue;				
			}
			if(line.charAt(0) == '>'){
			   Lines.add(new ArrayList<String>());
				i++;
			}
			Lines.get(i).add(line);
		}
		List <MyMat> Mlist = new ArrayList<MyMat>();
		List <String> name3 = new ArrayList<String>();
		Pattern p = Pattern.compile(">(.+)");
		for(int k = 0; k < Lines.size(); k++){
			List<String> lines = Lines.get(k);
			int l = 0;
			line = lines.get(0);
			List<String> str = Arrays.asList((line.split("\t")));
			Matcher m = p.matcher(str.get(0));
			if(!m.matches()){
				continue;
			}
			name3.add(m.group(1));
			int ncol = str.size()-1;
			int nrow = lines.size() -1;
			List <String> colname = new ArrayList<String>();
			List <String> rowname = new ArrayList<String>();
			
			for(i = 1; i < str.size(); i++){
			 colname.add(str.get(i));	
			}
			double[][] M = new double[nrow][ncol];
			for(l=1;l<=nrow;l++){
				line = lines.get(l);
				str = Arrays.asList(line.split("\t"));
				if(str.size() != ncol+1){
					throw new DataFormatException("Order3Tensor: file format is wrong!");
				}
				rowname.add(str.get(0));
				for(i=1; i<str.size(); i++){
					M[l-1][i-1] = Double.valueOf(str.get(i));		
				 }
			}
			Mlist.add(new MyMat(rowname, colname, M));
		}
		Set<String> commonRows = new TreeSet<String>(Mlist.get(0).getRowNames());
		for(i = 1; i < Mlist.size(); i++){
			Set<String> tmp = new HashSet<String>(Mlist.get(i).getRowNames());
			commonRows.retainAll(tmp);
		}
		Set<String> commonCols = new TreeSet<String>(Mlist.get(0).getColNames());
		for(i = 1; i < Mlist.size(); i++){
			Set<String> tmp = new HashSet<String>(Mlist.get(i).getColNames());
			commonCols.retainAll(tmp);
		}
		
		List <String> name1 = new ArrayList<String>(commonRows);
		List <String> name2 = new ArrayList<String>(commonCols);
		n1 =name1.size(); 
		n2 =name2.size(); 
		n3 =name3.size(); 
		setName1(name1);
		setName2(name2);
		setName3(name3);
		int  j, k; 
		M = new double[n1][n2][n3];
		for(i=0;i<n1;i++){
			for(j=0;j<n2;j++){
				for(k=0;k<n3;k++){
					M[i][j][k] = Mlist.get(k).get(name1.get(i), name2.get(j));	
				}
			}
		}
	}
	
	
	public String toString(){
		StringBuffer S = new StringBuffer();
		for(int i=0; i<n3; i++){	
			S.append(">"+getName3().get(i) + getOrder3Slice(i));
		}
		return S.toString();
	}
	
	public boolean containsName1(String name1){
		return name12index.containsKey(name1);
	}
	public int  Order1IndexOf(String name1){
		return name12index.get(name1);
	}
	public boolean containsName2(String name2){
		return name22index.containsKey(name2);
	}
	public int  Order2IndexOf(String name2){
		return name22index.get(name2);
	}public boolean containsName3(String name3){
		return name32index.containsKey(name3);
	}
	public int  Order3IndexOf(String name3){
		return name32index.get(name3);
	}
	
	public double get(int i, int j, int k){
		if(i >= n1 || j >= n2 || k >= n3){
			throw new IndexOutOfBoundsException();
		}
		return M[i][j][k];
	}
	public double get(String i, String j, String k){
		return M[name12index.get(i)][name22index.get(j)][name32index.get(k)];
	}
	public void set(int i, int j, int k, double d){
		M[i][j][k] = d;
	}
	public void set(String i, String j, String k, double d){
		M[name12index.get(i)][name22index.get(j)][name32index.get(k)] = d;	
	}
	public void setName1(List <String> v){
		name1 = new ArrayList<String>(v);
		if(name12index == null){
			name12index = new HashMap<String, Integer>();
		}else{
			name12index.clear();
		}
		int i;
		for(i=0;i<n1;i++){
			name12index.put(name1.get(i),i);
		}
	}
	public void setName2(List <String> v){
		name2 = new ArrayList<String>(v);
		if(name22index == null){
			name22index = new HashMap<String, Integer>();
		}else{
			name22index.clear();
		}
		int i;
		for(i=0;i<n2;i++){
			name22index.put(name2.get(i),i);
		}
	}
	public void setName3(List <String> v){
		name3 = new ArrayList<String>(v);
		if(name32index == null){
			name32index = new HashMap<String, Integer>();
		}else{
			name32index.clear();
		}
		int i;
		for(i=0;i<n3;i++){
			name32index.put(name3.get(i),i);
		}
	}
	public List<String> getName1(){
		return name1;
	}
	public List<String> getName2(){
		return name2;
	}
	public List<String> getName3(){
		return name3;
	}
	public int getDimOfOrder1(){
		return n1;
	}
	public int getDimOfOrder2(){
		return n2;
	}
	public int getDimOfOrder3(){
		return n3;
	}
	public MyMat getOrder1Slice(int i){
		MyMat S = new MyMat(name2, name3);
		int j,k;
		for(j = 0; j < n2; j++){
			for(k = 0; k < n3; k++){
				S.set(j, k, M[i][j][k]);
			}
		}	
		return S;
	}
	public MyMat getOrder2Slice(int j){
		MyMat S = new MyMat(name3, name1);
		int i,k;
		for(k = 0; k < n3; k++){	
				for(i = 0; i < n1; i++){
				S.set(k, i, M[i][j][k]);
			}
		}	
		return S;
	}
	public MyMat getOrder3Slice(int k){
		MyMat S = new MyMat(name1, name2);
		int i,j;
		for(i = 0; i < n1; i++){
			for(j = 0; j < n2; j++){	
				S.set(i, j, M[i][j][k]);
			}
		}	
		return S;
	}
	public MyMat getOrder1Slice(String s){
		return getOrder1Slice(Order1IndexOf(s));
	}
	public MyMat getOrder2Slice(String s){
		return getOrder2Slice(Order2IndexOf(s));
	}
	public MyMat getOrder3Slice(String s){
		return getOrder3Slice(Order3IndexOf(s));
	}
	public void reorderOrder1(List name1){
		List<Integer> v = new ArrayList<Integer>();
		int i,j,k;
		Set<Object>seen = new HashSet<Object>();
		for(i=0;i<name1.size();i++){
			if(name1.get(i) instanceof Integer){
				if((Integer)name1.get(i) < n1 && !seen.contains(name1.get(i)))
				v.add((Integer) name1.get(i));
				seen.add(name1.get(i));
				continue;	
			}
			if(name1.get(i) instanceof String){
				if(name12index.containsKey(name1.get(i)) && !seen.contains(name1.get(i)))
				v.add(name12index.get(name1.get(i)));	
				seen.add(name1.get(i));
				continue;	
			}
			throw new IllegalArgumentException("Need Integer or String");	
		}
		List <String> new_name1= new ArrayList<String>();
		n1 = v.size();
		double[][][] new_M = new double[n1][n2][n3];
		for(i=0;i<n1;i++){
		  new_name1.add(this.name1.get(v.get(i)));
		  for(j = 0; j< n2; j++){
			  for(k = 0; k< n3; k++){
			  new_M[i][j][k] = M[v.get(i)][j][k];
			  }	
		  }	
		}
		setName1(new_name1);
		M = new_M;	
	}
	public void reorderOrder2(List name2){
		List<Integer> v = new ArrayList<Integer>();
		int i,j,k;
		Set<Object>seen = new HashSet<Object>();
		for(i=0;i<name2.size();i++){
			if(name2.get(i) instanceof Integer){
				if((Integer)name2.get(i) < n2 && !seen.contains(name2.get(i)))
					v.add((Integer) name2.get(i));
						seen.add(name2.get(i));
						continue;	
			}
			if(name2.get(i) instanceof String){
				if(name22index.containsKey(name2.get(i)) && !seen.contains(name2.get(i)))
					v.add(name22index.get(name2.get(i)));	
					seen.add(name2.get(i));
					continue;	
				}
			throw new IllegalArgumentException("Need Integer or String");	
		}
		List <String> new_name2= new ArrayList<String>();
		n2 = v.size();
		double[][][] new_M = new double[n1][n2][n3];
			
		for(j=0;j<n2;j++){
			new_name2.add(this.name2.get(v.get(j)));
			for(i = 0; i< n1; i++){
				for(k = 0; k< n3; k++){
					new_M[i][j][k] = M[i][v.get(j)][k];
			   }	
			}
		}
		setName2(new_name2);
		M = new_M;	
	}
	
	public void reorderOrder3(List name3){
		List<Integer> v = new ArrayList<Integer>();
		int i,j,k;
		Set<Object>seen = new HashSet<Object>();
		for(i=0;i<name3.size();i++){
			if(name3.get(i) instanceof Integer){
				if((Integer)name3.get(i) < n3 && !seen.contains(name3.get(i)))
					v.add((Integer) name3.get(i));
						seen.add(name3.get(i));
						continue;	
			}
			if(name3.get(i) instanceof String){
				if(name32index.containsKey(name3.get(i)) && !seen.contains(name3.get(i)))
					v.add(name32index.get(name3.get(i)));	
					seen.add(name3.get(i));
					continue;	
				}
			throw new IllegalArgumentException("Need Integer or String");	
		}
		List <String> new_name3= new ArrayList<String>();
		n3 = v.size();
		double[][][] new_M = new double[n1][n2][n3];
			
		for(k=0;k<n3;k++){
			new_name3.add(this.name3.get(v.get(k)));
			for(i = 0; i< n1; i++){
				for(j = 0; j< n2; j++){
					new_M[i][j][k] = M[i][j][v.get(k)];
			   }	
			}
		}
		setName3(new_name3);
		M = new_M;	
	}
	
	public static Dist  getDistBetweenOrder1Slices(Order3Tensor T){
		Dist dist = new Dist(T.getName1());
		int i,j;
		for(i = 0; i < dist.size(); i++){
			for(j = 0; j < i; j++){
				dist.set(i, j, MyFunc.euclideanDist(T.getOrder1Slice(i).asList(),T.getOrder1Slice(j).asList()));
		    }
		}
		return dist;		
	}
	public static Dist  getDistBetweenOrder2Slices(Order3Tensor T){
		Dist dist = new Dist(T.getName2());
		int i,j;
		for(i = 0; i < dist.size(); i++){
			for(j = 0; j < i; j++){
				dist.set(i, j, MyFunc.euclideanDist(T.getOrder2Slice(i).asList(),T.getOrder2Slice(j).asList()));
		    }
		}
		return dist;	
	}
	public static Dist  getDistBetweenOrder3Slices(Order3Tensor T){
		Dist dist = new Dist(T.getName3());
		int i,j;
		for(i = 0; i < dist.size(); i++){
			for(j = 0; j < i; j++){
				dist.set(i, j, MyFunc.euclideanDist(T.getOrder3Slice(i).asList(),T.getOrder3Slice(j).asList()));
		    }
		}
		return dist;	
	}
	public static Dist  getDistBetweenSlices(Order3Tensor T, int sliceOrder){
		if(sliceOrder == 1){
			return getDistBetweenOrder1Slices(T);
		}else if(sliceOrder == 2){
			return getDistBetweenOrder2Slices(T);
		}else if(sliceOrder == 3){
			return getDistBetweenOrder3Slices(T);
		}else{
			return null;
		}
	}	
	
	public static MyMat collapseOrder2and3(Order3Tensor T){
		MyMat tmp = T.getOrder3Slice(0);
		List<String> newColnames = new ArrayList<String>();
		for(String s: tmp.getColNames()){
			newColnames.add(s + "0");
		}
		tmp.setColNames(newColnames);
		for(int i = 1; i < T.n3; i++){
			MyMat tmp2 = T.getOrder3Slice(i);
			List<String> newColnames2 = new ArrayList<String>();
			for(String s: tmp2.getColNames()){
				newColnames2.add(s + i);
			}
			tmp2.setColNames(newColnames2);
			tmp = tmp.bindCol(tmp2);
		}
		return tmp;
	}
	
	public static MyMat collapseOrder3and1(Order3Tensor T){
		MyMat tmp = T.getOrder1Slice(0);
		List<String> newColnames = new ArrayList<String>();
		for(String s: tmp.getColNames()){
			newColnames.add(s + "0");
		}
		tmp.setColNames(newColnames);
		for(int i = 1; i < T.n1; i++){
			MyMat tmp2 = T.getOrder1Slice(i);
			List<String> newColnames2 = new ArrayList<String>();
			for(String s: tmp2.getColNames()){
				newColnames2.add(s + i);
			}
			tmp2.setColNames(newColnames2);
			tmp = tmp.bindCol(tmp2);
		}
		return tmp;
	}
	public static MyMat collapseOrder1and2(Order3Tensor T){
		MyMat tmp = T.getOrder2Slice(0);
		List<String> newColnames = new ArrayList<String>();
		for(String s: tmp.getColNames()){
			newColnames.add(s + "0");
		}
		tmp.setColNames(newColnames);
		for(int i = 1; i < T.n2; i++){
			MyMat tmp2 = T.getOrder2Slice(i);
			List<String> newColnames2 = new ArrayList<String>();
			for(String s: tmp2.getColNames()){
				newColnames2.add(s + i);
			}
			tmp2.setColNames(newColnames2);
			tmp = tmp.bindCol(tmp2);
		}
		return tmp;
	}
	
	public void filterOrder1ByVariance(double d){
		SortedMap < Double, List<Integer>> sm = new TreeMap<Double, List<Integer>>();
		  int i;
		  for(i=0;i<n1;i++){
			  double s = MyFunc.sd(getOrder1Slice(i).asList()); 
			  if(sm.containsKey(s)){
					 (sm.get(s)).add(i);
				 }else{
					 List<Integer> tmp = new ArrayList<Integer>();
					 tmp.add(i);
					 sm.put(s, tmp);
				 }
		  }
		  
		  List< Double > tmp = new ArrayList<Double>(sm.keySet());
		  Collections.reverse(tmp);
		 List<Integer> tmp2 = new ArrayList<Integer>(); 
		  if(d < 1){
		    d = Math.round(n1*d);
		  } 
		 for(i=0;;i++){
			 List<Integer> tmp3 = sm.get(tmp.get(i));
			 if(tmp3.size()+tmp2.size() <= d){
			  tmp2.addAll(sm.get(tmp.get(i)));
			 }else{
				break;
			 }
		 }
		 reorderOrder1(tmp2);
	}
	public void filterOrder2ByVariance(double d){
		SortedMap < Double, List<Integer>> sm = new TreeMap<Double, List<Integer>>();
		  int i;
		  for(i=0;i<n2;i++){
			  double s = MyFunc.sd(getOrder2Slice(i).asList()); 
			  if(sm.containsKey(s)){
					 (sm.get(s)).add(i);
				 }else{
					 List<Integer> tmp = new ArrayList<Integer>();
					 tmp.add(i);
					 sm.put(s, tmp);
				 }
		  }
		  
		  List< Double > tmp = new ArrayList<Double>(sm.keySet());
		  Collections.reverse(tmp);
		 List<Integer> tmp2 = new ArrayList<Integer>(); 
		  if(d < 1){
		    d = Math.round(n2*d);
		  } 
		 for(i=0;;i++){
			 List<Integer> tmp3 = sm.get(tmp.get(i));
			 if(tmp3.size()+tmp2.size() <= d){
			  tmp2.addAll(sm.get(tmp.get(i)));
			 }else{
				break;
			 }
		 }
		 reorderOrder2(tmp2);
	}
	public void filterOrder3ByVariance(double d){
		SortedMap < Double, List<Integer>> sm = new TreeMap<Double, List<Integer>>();
		  int i;
		  for(i=0;i<n3;i++){
			  double s = MyFunc.sd(getOrder3Slice(i).asList()); 
			  if(sm.containsKey(s)){
					 (sm.get(s)).add(i);
				 }else{
					 List<Integer> tmp = new ArrayList<Integer>();
					 tmp.add(i);
					 sm.put(s, tmp);
				 }
		  }
		  
		  List< Double > tmp = new ArrayList<Double>(sm.keySet());
		  Collections.reverse(tmp);
		 List<Integer> tmp2 = new ArrayList<Integer>(); 
		  if(d < 1){
		    d = Math.round(n3*d);
		  } 
		 for(i=0;;i++){
			 List<Integer> tmp3 = sm.get(tmp.get(i));
			 if(tmp3.size()+tmp2.size() <= d){
			  tmp2.addAll(sm.get(tmp.get(i)));
			 }else{
				break;
			 }
		 }
		 reorderOrder3(tmp2);
	}
	
	public ArrayList<Double> asList(){
		ArrayList<Double> tmp = new ArrayList<Double>();
		for(int i = 0; i < n3; i++){
			tmp.addAll(getOrder3Slice(i).asList());
		}
		return tmp;
	}
	

	public void printDimension(){
		System.out.println("order1:\t" + n1);
		System.out.println("order2:\t" + n2);
		System.out.println("order3:\t" + n3);
	}
	
	
	public void printStatistics(){
		List <Double> l  = asList();
 		System.out.println("mean:\t" + MyFunc.mean(l));
		System.out.println("sd:\t" + MyFunc.sd(l));
		System.out.println("max:\t" + MyFunc.max(l));
		System.out.println("min:\t" + MyFunc.min(l));
	}
	
	public void print(){
		System.out.print(this);
	}
	
	public Order3TensorViewer getViewer(){
		return new  Order3TensorViewer(this);
	}
	
	public MyMat getOrder1MeanSlice(List <String> name1){
		MyMat M = new MyMat(name2,name3);
		for(String s2: name2){
			for(String s3: name3){
				double tmp = 0;
				for(String s1: name1){
					tmp += get(s1,s2,s3);
				}
				M.set(s2,s3,tmp);
			}
		}
		return M;
	}
	
	public MyMat getOrder2MeanSlice(List <String> name2){
		MyMat M = new MyMat(name3,name1);
		for(String s3: name3){
			for(String s1: name1){
				double tmp = 0;
				for(String s2: name2){
					tmp += get(s1,s2,s3);
				}
				M.set(s3,s1,tmp);
			}
		}
		return M;
	}
	
	public MyMat getOrder3MeanSlice(List <String> name3){
		MyMat M = new MyMat(name1,name2);
		for(String s1: name1){
			for(String s2: name2){
				double tmp = 0;
				for(String s3: name3){
					tmp += get(s1,s2,s3);
				}
				M.set(s1,s2,tmp);
			}
		}
		return M;
	}
	
	public double max(){
		double tmp = -Double.MAX_VALUE;
		for(int i = 0; i < n1; i++){
			for(int j = 0; j < n2; j++){
				for(int k = 0; k < n3; k++){
					if(tmp < get(i,j,k)){
						tmp = get(i,j,k);
					}
				}
			}
		}
		return tmp;
	}
	
	public double min(){
		double tmp = Double.MAX_VALUE;
		for(int i = 0; i < n1; i++){
			for(int j = 0; j < n2; j++){
				for(int k = 0; k < n3; k++){
					if(tmp > get(i,j,k)){
						tmp = get(i,j,k);
					}
				}
			}
		}
		return tmp;
	}
	
	
}
