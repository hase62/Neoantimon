package utility;


import java.util.*;
import java.util.zip.DataFormatException;
import java.io.*;

public class DoubleMat extends MyMat  {
	private static final long serialVersionUID = 1L;
	protected Double M[][];
	public DoubleMat(String infile) throws IOException, DataFormatException{
		colname = new ArrayList<String>();
		rowname = new ArrayList<String>();
		colname2index = new HashMap<String, Integer>();
		rowname2index = new HashMap<String, Integer>();
		BufferedReader inputStream = new BufferedReader(new FileReader(infile));
		List<String> lines = new ArrayList<String>();
		String line;
		nrow = -1;
		while((line = inputStream.readLine()) != null){
			if(line.charAt(0) == '#'){
				continue;				
			}
			lines.add(line);
			nrow++;
		}
		
		int i;	
		
		int l = 0;
		line = lines.get(0);
		List<String> str = Arrays.asList((line.split("\t")));
		ncol = str.size()-1;
		for(i = 1; i < str.size(); i++){
		 colname.add(str.get(i));	
		}
			
		M = new Double[nrow][ncol];
		for(l=1;l<=nrow;l++){
			line = lines.get(l);
			str = Arrays.asList(line.split("\t"));
			if(str.size() != ncol+1){
				throw new DataFormatException("DoubleMat: file format is wrong!");
			}
			rowname.add(str.get(0));
			for(i=1; i<str.size(); i++){
				try{
					M[l-1][i-1] = Double.valueOf(str.get(i));
				}catch(NumberFormatException e){
					M[l-1][i-1] = Double.NaN;
				}
			 }
		}
	   for(i=0;i<ncol;i++){
		 colname2index.put(colname.get(i),i);
	   }
	  for(i=0;i<nrow;i++){
		 rowname2index.put(rowname.get(i),i);
	  }
	  inputStream.close();
	}
	public DoubleMat(int i, int j){
		nrow = i;
		ncol = j;
		colname = new ArrayList<String>();
		rowname = new ArrayList<String>();
		colname2index = new HashMap<String, Integer>();
		rowname2index = new HashMap<String, Integer>();
		M = new Double[i][j];
		for(i=0;i<nrow;i++){
			for(j=0;j<ncol;j++){
		      M[i][j] = 0.0;
		    }	
		}
		for(i=1;i<=ncol;i++){
			colname.add("c"+i);
		}
		for(j=1;j<=nrow;j++){
			rowname.add("r"+j);
		}
		for(i=0;i<ncol;i++){
			colname2index.put(colname.get(i),i);
		 }
		for(i=0;i<nrow;i++){
			rowname2index.put(rowname.get(i),i);
	    }
	}
	
	public DoubleMat(List <String> row, List <String> col){
		nrow = row.size();
		ncol = col.size();
		M = new Double[nrow][ncol];
		for(int i=0;i<nrow;i++){
			for(int j=0;j<ncol;j++){
		      M[i][j] = 0.0;
		    }	
		}
		setColNames(col);
		setRowNames(row);		
	}

	public DoubleMat(List <String> row, List <String> col, double[][] M){
		nrow = row.size();
		ncol = col.size();
		this.M = new Double[nrow][ncol];
		for(int i=0;i<nrow;i++){
			for(int j=0;j<ncol;j++){
		      this.M[i][j] = M[i][j];
		    }	
		}
		setColNames(col);
		setRowNames(row);		
	}
	
	
	
	public DoubleMat(DoubleMat m){
		nrow = m.nrow;
		ncol = m.ncol;
		colname2index = new HashMap<String, Integer>(m.colname2index);
		rowname2index = new HashMap<String, Integer>(m.rowname2index);
		colname = new ArrayList<String>(m.colname);
		rowname = new ArrayList<String>(m.rowname);
		M = new Double[nrow][ncol];
		int i,j;
		for(i=0;i<nrow;i++){
			for(j=0;j<ncol;j++){
				M[i][j] = m.M[i][j];
			}	
		}
	}
	public DoubleMat(String colname, Map<String, Double> rowname2value){
		nrow = rowname2value.size();
		ncol = 1;
		this.colname = new ArrayList<String>();
		rowname = new ArrayList<String>();
		colname2index = new HashMap<String, Integer>();
		rowname2index = new HashMap<String, Integer>();
		M = new Double[nrow][ncol];
		for(int i=0;i<nrow;i++){
			for(int j=0;j<ncol;j++){
		      M[i][j] = 0.0;
		    }	
		}
		List <String> tmp = new ArrayList<String>();
		tmp.add(colname);
		setColNames(tmp);
		setRowNames(new ArrayList<String>(rowname2value.keySet()));		
		for(Map.Entry<String, Double> e: rowname2value.entrySet()){
			set(e.getKey(), colname, e.getValue());	
		}	
	}
	
	
	public boolean containsColName(String colname){
		return colname2index.containsKey(colname);
	}
	public boolean containsRowName(String rowname){
		return rowname2index.containsKey(rowname);
	}
	public int colIndexOf(String colname){
		return colname2index.get(colname);
	}
	public int rowIndexOf(String rowname){
		return rowname2index.get(rowname);
	}
	public double get(int i, int j){
		if(i >= nrow || j >= ncol){
			throw new IndexOutOfBoundsException();
		}
		return M[i][j];
	}
	public double get(String i, String j){
		return M[rowname2index.get(i)][colname2index.get(j)];
	}
	public void set(int i, int j, double d){
		M[i][j] = d;
	}
	public void set(String i, String j, double d){
		M[rowname2index.get(i)][colname2index.get(j)] = d;	
	}
	public void setRowNames(List <String> v){
		rowname = new ArrayList<String>(v);
		if(rowname2index == null){
			rowname2index = new HashMap<String, Integer>();
		}else{
			rowname2index.clear();
		}
		int i;
		for(i=0;i<nrow;i++){
			rowname2index.put(rowname.get(i),i);
		}
	}
	public void setColNames(List <String> v){
		colname = new ArrayList<String>(v);
		if(colname2index == null){
			colname2index = new HashMap<String, Integer>();
		}else{
			colname2index.clear();
		}
		int i;
		for(i=0;i<ncol;i++){
			colname2index.put(colname.get(i),i);
		}
	}
	public List<String> getRowNames(){
		return rowname;
	}
	public List<String> getColNames(){
		return colname;
	}
	public int rowSize(){
		return nrow;
	}
	public int colSize(){
		return ncol;
	}
	public void transpose(){
		int i,j;
		Double[][] tmp = new Double[ncol][nrow];
		for(i=0;i<nrow;i++){
			for(j=0;j<ncol;j++){
				tmp[j][i] = M[i][j];
		    }
		}
		M = tmp;
		i = nrow;
		nrow = ncol;
		ncol = i;
		List <String> tmp2 = new ArrayList<String>(rowname);
		setRowNames(colname);
		setColNames(tmp2);
	}
	@SuppressWarnings("unchecked")
	public void reorderRows(List row){
		List<Integer> v = new ArrayList<Integer>();
		int i,j;
		Set<Object>seen = new HashSet<Object>();
		for(i=0;i<row.size();i++){
			if(row.get(i) instanceof Integer){
				if((Integer)row.get(i) < nrow && !seen.contains(row.get(i)))
				v.add((Integer) row.get(i));
				seen.add(row.get(i));
				continue;	
			}
			if(row.get(i) instanceof String){
				if(rowname2index.containsKey(row.get(i)) && !seen.contains(row.get(i)))
				v.add(rowname2index.get(row.get(i)));	
				seen.add(row.get(i));
				continue;	
			}
			throw new IllegalArgumentException("Need Integer or String");	
		}
		List <String> new_rowname= new ArrayList<String>();
		Double[][] new_M = new Double[v.size()][ncol];
		nrow = v.size();
		for(i=0;i<nrow;i++){
		  new_rowname.add(rowname.get(v.get(i)));
		  for(j = 0; j< ncol; j++){
			  new_M[i][j] = M[v.get(i)][j];
		   }	
		}
		setRowNames(new_rowname);
		M = new_M;	
	}
	@SuppressWarnings("unchecked")
	public void reorderCols(List col){
		List<Integer> v = new ArrayList<Integer>();
		int i,j;
		Set<Object>seen = new HashSet<Object>();
		for(i=0;i<col.size();i++){
			if(col.get(i) instanceof Integer){
				if((Integer)col.get(i) < ncol && !seen.contains(col.get(i)))
				v.add((Integer) col.get(i));
				seen.add(col.get(i));
				continue;	
			}
			if(col.get(i) instanceof String){
				if(colname2index.containsKey(col.get(i)) && !seen.contains(col.get(i)))
				v.add(colname2index.get(col.get(i)));
				seen.add(col.get(i));
				continue;	
			}
			throw new IllegalArgumentException("Need Integer or String");	
		}
		List <String> new_colname= new ArrayList<String>();
		Double[][] new_M = new Double[nrow][v.size()];
		ncol = v.size();
		for(j=0;j<ncol;j++){
		  new_colname.add(colname.get(v.get(j)));
		  for(i = 0; i< nrow; i++){
			  new_M[i][j] = M[i][v.get(j)];
		   }	
		}
		setColNames(new_colname);
		M = new_M;	
	}
	public DoubleMat getSubMatrix(List row, List col){
		List<Integer> v = new ArrayList<Integer>();
		int i,j;
		Set<Object>seen = new HashSet<Object>();
		
		for(i=0;i<row.size();i++){
			if(row.get(i) instanceof Integer){
				if((Integer)row.get(i) < nrow && !seen.contains(row.get(i)))
				v.add((Integer) row.get(i));	
				seen.add(row.get(i));
				continue;	
			}
			if(row.get(i) instanceof String){
				if(rowname2index.containsKey(row.get(i)) && !seen.contains(row.get(i)))
				v.add(rowname2index.get(row.get(i)));
				seen.add(row.get(i));
				continue;	
			}
			throw new IllegalArgumentException("Need Integer or String");	
		}
		List <String> new_rowname= new ArrayList<String>();
		int new_nrow = v.size();
		for(i=0;i<new_nrow;i++){
		  new_rowname.add(rowname.get(v.get(i)));
		}
		List<Integer> u = new ArrayList<Integer>();
		seen.clear();
		for(j=0;j<col.size();j++){
			if(col.get(j) instanceof Integer){
				if((Integer)col.get(j) < ncol && !seen.contains(col.get(j)))
				u.add((Integer) col.get(j));
				seen.add(col.get(j));
				continue;	
			}
			if(col.get(j) instanceof String){
				if(colname2index.containsKey(col.get(j)) && !seen.contains(col.get(j)))
				u.add(colname2index.get(col.get(j)));	
				seen.add(col.get(j));
				continue;	
			}
			throw new IllegalArgumentException("Need Integer or String");	
		}
		List <String> new_colname= new ArrayList<String>();
		int new_ncol = u.size();
		for(j=0;j<new_ncol;j++){
		  new_colname.add(colname.get(u.get(j)));
		}
		DoubleMat newMat = new DoubleMat(new_rowname, new_colname);
		for(i = 0; i < new_nrow; i++){
			for(j = 0; j < new_ncol; j++){
				 newMat.M[i][j] = M[v.get(i)][u.get(j)];
				
			}
		}
		return newMat;
	}
	@SuppressWarnings("unchecked")
	public DoubleMat getSubMatByRow(List row){
		List<Integer> v = new ArrayList<Integer>();
		int i;
		Set<Object>seen = new HashSet<Object>();
		for(i=0;i<row.size();i++){
			if(row.get(i) instanceof Integer){
				if((Integer)row.get(i) < nrow && !seen.contains(row.get(i)))
				v.add((Integer) row.get(i));
				seen.add(row.get(i));
				continue;	
			}
			if(row.get(i) instanceof String){
				if(rowname2index.containsKey(row.get(i)) && !seen.contains(row.get(i)))
				v.add(rowname2index.get(row.get(i)));
				seen.add(row.get(i));
				continue;	
			}
			throw new IllegalArgumentException("Need Integer or String");	
		}
		List <String> new_rowname= new ArrayList<String>();
		for(i=0;i<v.size();i++){
		  new_rowname.add(rowname.get(v.get(i)));
		}
		
		DoubleMat new_M = new DoubleMat(new_rowname, colname);
		for(String s: new_rowname){
			for(String t: colname){
				new_M.set(s, t, get(s,t));		
			}
		}
		return new_M;	
	}
	
	
	@SuppressWarnings("unchecked")
	public DoubleMat getSubMatByCol(List col){
		List<Integer> v = new ArrayList<Integer>();
		int i,j;
		Set<Object>seen = new HashSet<Object>();
		for(i=0;i<col.size();i++){
			if(col.get(i) instanceof Integer){
				if((Integer)col.get(i) < ncol && !seen.contains(col.get(i)))
				v.add((Integer) col.get(i));
				seen.add(col.get(i));
				continue;	
			}
			if(col.get(i) instanceof String){
				if(colname2index.containsKey(col.get(i)) && !seen.contains(col.get(i)))
				v.add(colname2index.get(col.get(i)));
				seen.add(col.get(i));
				continue;	
			}
			throw new IllegalArgumentException("Need Integer or String");	
		}
		List <String> new_colname= new ArrayList<String>();
		for(j=0;j<v.size();j++){
		  new_colname.add(colname.get(v.get(j)));
		}
		DoubleMat new_M = new DoubleMat(rowname, new_colname);
		for(String s: rowname){
			for(String t: new_colname){
				new_M.set(s, t, get(s,t));
			}
			
		}
		return new_M;		
	}
	
	public List<Double> getRow(int i){
		List<Double> tmp = new ArrayList<Double>();
		int j;
		for(j = 0; j < ncol; j++){
			tmp.add(M[i][j]);
		}
		return tmp;
	}
	public List<Double> getRow(String i){
		List<Double> tmp = new ArrayList<Double>();
		int j;
		for(j = 0; j < ncol; j++){
			tmp.add(M[rowname2index.get(i)][j]);
		}
		return tmp;
	}
	public Map<String,Double> getRowMap(int i){
		Map<String, Double> tmp = new HashMap<String, Double>();
		int j;
		for(j = 0; j < ncol; j++){
			tmp.put(colname.get(j), M[i][j]);
		}
		return tmp;
	}
	public Map<String,Double> getRowMap(String i){
		Map<String, Double> tmp = new HashMap<String, Double>();
		int j;
		for(j = 0; j < ncol; j++){
			tmp.put(colname.get(j),M[rowname2index.get(i)][j]);
		}
		return tmp;
	}
	public List<Double> getCol(int j){
		List<Double> tmp = new ArrayList<Double>();
		int i;
		for(i = 0; i < nrow; i++){
			tmp.add(M[i][j]);
		}
		return tmp;
	}
	public List<Double> getCol(String j){
		List<Double> tmp = new ArrayList<Double>();
		int i;
		for(i = 0; i < nrow; i++){
			tmp.add(M[i][colname2index.get(j)]);
		}
		return tmp;
	}
	public Map<String,Double> getColMap(int j){
		Map<String, Double> tmp = new HashMap<String, Double>();
		int i;
		for(i = 0; i < nrow; i++){
			tmp.put(rowname.get(i), M[i][j]);
		}
		return tmp;
	}
	public Map<String,Double> getColMap(String j){
		Map<String, Double> tmp = new HashMap<String, Double>();
		int i;
		for(i = 0; i < nrow;i++){
			tmp.put(rowname.get(i),M[i][colname2index.get(j)]);
		}
		return tmp;
	}
	public void filterRowByVariance(double d){
		 SortedMap < Double, List<Integer>> sm = new TreeMap<Double, List<Integer>>();
		  int i;
		  for(i=0;i<nrow;i++){
			  double s = MyFunc.sd(getRow(i)); 
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
		    d = Math.round(nrow*d);
		  } 
		 for(i=0;;i++){
			 List<Integer> tmp3 = sm.get(tmp.get(i));
			 if(tmp3.size()+tmp2.size() <= d){
			  tmp2.addAll(sm.get(tmp.get(i)));
			 }else{
				break;
			 }
		 }
		 reorderRows(tmp2);
	}
	public DoubleMat bindRow(DoubleMat m){
		List <String> newRowname = new ArrayList<String>(rowname);
		newRowname.addAll(m.rowname);
		newRowname = MyFunc.uniq(newRowname);
		List <String> newColname = MyFunc.isect(colname, m.colname);
		DoubleMat tmp = new DoubleMat(newRowname, newColname);
		for(String t: tmp.colname){
			for(String s: tmp.rowname){
				if(containsRowName(s)){	
					tmp.set(s, t, get(s, t));
				}else{
					tmp.set(s, t, m.get(s, t));
				}
			}
		}	
		return tmp;
	}
	public DoubleMat addRow(String rowname, Map<String, Double> value){
		List <String> newRowname = new ArrayList<String>(this.rowname);
		newRowname.add(rowname);
		DoubleMat tmp = new DoubleMat(newRowname, colname);
		for(String t: colname){
			for(String s: this.rowname){
				tmp.set(s, t, get(s, t));
			}
			if(value.containsKey(t)){
				tmp.set(rowname, t, value.get(t));
			}
		}
		return tmp;
	}
	public DoubleMat bindCol(DoubleMat m){
		List <String> newColname = new ArrayList<String>(colname);
		newColname.addAll(m.colname);
		newColname = MyFunc.uniq(newColname);
		List <String> newRowname = MyFunc.isect(rowname, m.rowname);
		DoubleMat tmp = new DoubleMat(newRowname, newColname);
		for(String s: tmp.rowname){
			for(String t: tmp.colname){
				if(containsColName(t)){
				tmp.set(s, t, get(s, t));
			}else{
				tmp.set(s, t, m.get(s, t));
				}
			}
		}
		return tmp;
	}
	public DoubleMat addCol(String colname, Map<String, Double> value){
		List <String> newColname = new ArrayList<String>(this.colname);
		newColname.add(colname);
		DoubleMat tmp = new DoubleMat(rowname, newColname);
		for(String s: rowname){
			for(String t: this.colname){
				tmp.set(s, t, get(s, t));
			}
			if(value.containsKey(s)){
				tmp.set(s, colname, value.get(s));
			}
		}
		return tmp;
	}
	public DoubleMat bind(DoubleMat m){
		if(!MyFunc.isect(getRowNames(), m.getRowNames()).isEmpty()){
			return bindCol(m);
		}
		if(!MyFunc.isect(getRowNames(),m.getColNames()).isEmpty()){
			DoubleMat tmp = new DoubleMat(m);
			tmp.transpose();
			return bindCol(tmp);
		}
		if(!MyFunc.isect(getColNames(), m.getRowNames()).isEmpty()){
			DoubleMat tmp = new DoubleMat(m);
			tmp.transpose();
			return bindRow(tmp);
		}
		if(!MyFunc.isect(getColNames(),m.getColNames()).isEmpty()){
			return bindRow(m);
		}
		return this;		
	}
	
	public DoubleMat add(String name, Map<String, Double> value){
		if(!MyFunc.isect(getRowNames(),new ArrayList<String>(value.keySet())).isEmpty()){
			return addCol(name, value);
		}
		if(!MyFunc.isect(getColNames(), new ArrayList<String>(value.keySet())).isEmpty()){
			return addRow(name, value);
		}
		return this;
	}
	
	
	
	public DoubleMat multiply(DoubleMat m){
		if(ncol != m.nrow){
			throw new RuntimeException("colsize of the first DoubleMat must be equal to rowSize of the second DoubleMat");
		}
		DoubleMat tmp = new DoubleMat(nrow, m.ncol);
		  int i,j,k;
		  for(i = 0; i < nrow; i++){
		    for(j = 0; j < m.ncol; j++){
		      double tmp2 = 0;
		      for(k = 0; k < ncol; k++){
		        tmp2 += M[i][k] * m.get(k,j);
		      }
		      tmp.M[i][j] = tmp2;
		    }
		  }
		 tmp.setRowNames(rowname);
		 tmp.setColNames(m.colname);
		 return tmp;
	}
	public void normalizeRows(){
		int i,j;
		for(i=0;i<nrow;i++){
		  List <Double> tmp = getRow(i);
		  double tmp2 = MyFunc.mean(tmp);
		  double tmp3 = MyFunc.sd(tmp);
		  for(j=0;j<ncol;j++){
		     M[i][j] = (M[i][j] - tmp2)/tmp3;
		  }
		}		
	}
	public void normalizeCols(){
		int i,j;
		for(j=0;j<ncol;j++){
		  List <Double> tmp = getCol(j);
		  double tmp2 = MyFunc.mean(tmp);
		  double tmp3 = MyFunc.sd(tmp);
		  for(i=0;i<nrow;i++){
		     M[i][j] = (M[i][j] - tmp2)/tmp3;
		  }
		}		
	}
	public void shuffleRows(){
	    List<String> new_rowname = new ArrayList<String>(rowname);
	    List<String> old_rowname = new ArrayList<String>(colname);
	    Collections.shuffle(new_rowname);
	    reorderRows(new_rowname);
	    setRowNames(old_rowname);
	}
	public void shuffleCols(){
	    List<String> new_colname = new ArrayList<String>(colname);
	    List<String> old_colname = new ArrayList<String>(colname);
	    Collections.shuffle(new_colname);
	    reorderCols(new_colname);	
	    setColNames(old_colname);
	}
	
	
	public void sortColsByValue(String rowname){
		Map <String,Double> tmp = getRowMap(rowname);
		List <String> tmp2 = MyFunc.sortKeysByAscendingOrderOfValues(tmp);
		reorderCols(tmp2);
	}
	public void sortRowsByValue(String colname){
		Map <String,Double> tmp = getRowMap(colname);
		List <String> tmp2 = MyFunc.sortKeysByAscendingOrderOfValues(tmp);
		reorderRows(tmp2);
	}
	
	public void sortRowsByRowNames(){
	    List<String> new_rowname = new ArrayList<String>(rowname);
	    Collections.sort(new_rowname);
	    reorderRows(new_rowname); 		
	}
	public void sortRowsByVariance(){
		 SortedMap < Double, List<Integer>> sm = new TreeMap<Double, List<Integer>>();
		  int i;
		  for(i=0;i<nrow;i++){
			 double s = MyFunc.sd(getRow(i)); 
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
		 for(i=0;i < tmp.size() ;i++){
		    tmp2.addAll(sm.get(tmp.get(i)));
		  }
		 reorderRows(tmp2);
		
	}
	public void print(String outfile) throws IOException {
		PrintWriter os = new PrintWriter(new FileWriter(outfile));
		os.println("\t"+  MyFunc.join("\t", getColNames())); 
		int i,j;
		for(i=0; i<nrow; i++){
			List<Double> tmp = getRow(i);
			List<String> tmp2 = new ArrayList<String>();
			for(j=0;j<ncol;j++){
				tmp2.add( (tmp.get(j)).toString() );
			}	
			os.println(rowname.get(i) + "\t"+  MyFunc.join("\t",tmp2));
		}
		os.flush();
		os.close();
	}
	public void print() throws IOException {
		PrintWriter os = new PrintWriter(System.out);
		os.println("\t"+  MyFunc.join("\t", getColNames())); 
		int i,j;
		for(i=0; i<nrow; i++){
			List<Double> tmp = getRow(i);
			List<String> tmp2 = new ArrayList<String>();
			for(j=0;j<ncol;j++){
				tmp2.add( (tmp.get(j)).toString() );
			}	
			os.println(rowname.get(i) + "\t"+  MyFunc.join("\t",tmp2));
		}
		os.flush();
		os.close();
	}
	public String toString(){
		StringBuffer S = new StringBuffer("\t"+  MyFunc.join("\t", getColNames()) + "\n");
		for(int i=0; i<nrow; i++){	
			S.append(getRowNames().get(i) + "\t"+  MyFunc.join("\t",MyFunc.toString(getRow(i))) + "\n");
		}
		return S.toString();
	}
	
	public boolean equals(DoubleMat m){
		if(colname != m.colname || rowname != m.rowname){
			return false;	
		}
		int i,j;
		for(i=0;i<nrow;i++){
			for(j=0;j<ncol;j++){
				if(M[i][j] == m.M[i][j]){
				  return false;
				}
			}
		}
		return true;
	}
	public List<Double> asList(){
		List <Double> v = new ArrayList<Double>();
		int i,j;
		for(i=0;i<nrow;i++){
			for(j=0;j<ncol;j++){
				v.add(M[i][j]);
			}
		}	
		return v;
	}
	
	
	public Map <String, Double> asMap(){
		Map <String, Double> m = new HashMap<String, Double>();
		for(String s: rowname){
			for(String t: colname){
				m.put(s + "\t" + t, get(s,t));	
			}
		}
		return m;		
	}
	
	
	public List<Double> getRowMeans(List <String> rows){
		 List <Double> v = new ArrayList<Double>();
		  int i,j;
		  for(j=0; j<ncol; j++){
		    double tmp = 0;
		    for(i=0; i<rows.size(); i++){
		      tmp += get(rowname2index.get(rows.get(i)),j);
		    }
		    v.add(tmp/rows.size()) ;
		  }
		  return v;
	}
	public List<Double> getRowMeans(){
		 List <Double> v = new ArrayList<Double>();
		  int i,j;
		  for(j=0; j<ncol; j++){
		    double tmp = 0;
		    for(i=0; i<nrow; i++){
		      tmp += get(i,j);
		    }
		    v.add(tmp/nrow) ;
		  }
		  return v;
	}
	public List<Double> getColMeans(List <String> cols){
		 List <Double> v = new ArrayList<Double>();
		  int i,j;
		  for(i=0; i<nrow; i++){
		    double tmp = 0;
		    for(j=0; j<cols.size(); j++){
		      tmp += get(rowname2index.get(cols.get(i)),j);
		    }
		    v.add(tmp/cols.size()) ;
		  }
		  return v;
	}
	public List<Double> getColMeans(){
		 List <Double> v = new ArrayList<Double>();
		  int i,j;
		  for(i=0; i<nrow; i++){
		    double tmp = 0;
		    for(j=0; j<ncol; j++){
		      tmp += get(i,j);
		    }
		    v.add(tmp/ncol) ;
		  }
		  return v;
	}
		
	public DoubleMat getCovMatForRow(){
		DoubleMat M = new DoubleMat(this);
		List <Double>  ColMean =  M.getColMeans();
		for(int i = 0; i <  M.rowSize(); i++){
			for(int j = 0; j <  M.colSize(); j++){
				M.set(i, j, M.get(i, j)-ColMean.get(i));
			}	
		}
		DoubleMat Mt = new DoubleMat(M);
		Mt.transpose();
		M = M.multiply(Mt);
		for(int i = 0; i <  M.rowSize(); i++){
			for(int j = 0; j <  M.colSize(); j++){
				M.set(i, j, M.get(i, j)/ncol);
			}
		}
		return M;
		
	}
		
	public DoubleMat getCovMatForCol(){
		DoubleMat M = new DoubleMat(this);
		List <Double>  RowMean = M.getRowMeans();
		for(int i = 0; i <  M.rowSize(); i++){
			for(int j = 0; j <  M.colSize(); j++){
				M.set(i, j, M.get(i, j)-RowMean.get(j));
			}	
		}
		DoubleMat Mt = new DoubleMat(M);
		Mt.transpose();
		M =  Mt.multiply(M);
		for(int i = 0; i <  M.rowSize(); i++){
			for(int j = 0; j <  M.colSize(); j++){
				M.set(i, j, M.get(i, j)/nrow);
			}
		}
		return M;
	}
	
	public static Dist getDistBetweenRows(DoubleMat m, char type){
		 Dist dist = new Dist(m.getRowNames());
		 int i,j;
		 switch(type){
			case 'c':
				DoubleMat copy_m = new DoubleMat(m);
				copy_m.normalizeRows();
				for(i = 0; i < dist.size(); i++){
					for(j = 0; j < i; j++){
						dist.set(i, j, MyFunc.pearsonCorrelationForNormarizedList(copy_m.getRow(i),copy_m.getRow(j)));
				    }
				}
				dist.setDiagonalElement(1);
				break;
			case 'C':
				for(i = 0; i < dist.size(); i++){
					for(j = 0; j < i; j++){
						dist.set(i, j,  MyFunc.pearsonCorrelationForNormarizedList(m.getRow(i),m.getRow(j)));
				    }
				}
				dist.setDiagonalElement(1);
				break;
			case 'e':
				for(i = 0; i < dist.size(); i++){
					for(j = 0; j < i; j++){
						dist.set(i, j, MyFunc.euclideanDist(m.getRow(i),m.getRow(j)));
				    }
				}
				break;
			default:
				throw new IllegalArgumentException("type must be 'c','C', or 'e'");
		}
		return dist;
	}
	
	public static Dist getDistBetweenCols(DoubleMat m, char type){
		 Dist dist = new Dist(m.getColNames());
		 int i,j;
		 switch(type){
			case 'c':
				DoubleMat copy_m = new DoubleMat(m);
				copy_m.normalizeCols();
				for(i = 0; i < dist.size(); i++){
					for(j = 0; j < i; j++){
						dist.set(i, j,  MyFunc.pearsonCorrelationForNormarizedList(copy_m.getCol(i),copy_m.getCol(j)));
				    }
				}
				dist.setDiagonalElement(1);
				break;
			case 'C':
				for(i = 0; i < dist.size(); i++){
					for(j = 0; j < i; j++){
						dist.set(i, j, MyFunc.pearsonCorrelationForNormarizedList(m.getCol(i),m.getCol(j)));
				    }
				}
				dist.setDiagonalElement(1);
				break;
			case 'e':
				for(i = 0; i < dist.size(); i++){
					for(j = 0; j < i; j++){
						dist.set(i, j, MyFunc.euclideanDist(m.getCol(i),m.getCol(j)));
				    }
				}
				break;
			default:
				throw new IllegalArgumentException("type must be 'c','C', or 'e'");
		}
		return dist;
	}
	

	public void printDimension(){
		System.out.println("ncol:\t" + ncol);
		System.out.println("nrow:\t" + nrow);
	}
	
	
	public void printStatistics(){
		List <Double> l  = asList();
 		System.out.println("mean:\t" + MyFunc.mean(l));
		System.out.println("sd:\t" + MyFunc.sd(l));
		System.out.println("max:\t" + MyFunc.max(l));
		System.out.println("min:\t" + MyFunc.min(l));
	}
	












}
