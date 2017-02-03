package utility;
import java.util.*;
import java.util.zip.DataFormatException;
import java.io.*;


public class StringMat implements Serializable{
	private static final long serialVersionUID = 8715882894418912249L;
	protected String M[][];
	protected Map<String, Integer> colname2index;
	protected Map<String, Integer> rowname2index;
	protected int ncol;
	protected int nrow;
	protected List <String> colname;
	protected List <String> rowname;
	public StringMat(String infile) throws IOException, DataFormatException{
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
			
		M = new String[nrow][ncol];
		for(l=1;l<=nrow;l++){
			line = lines.get(l);
			str = new ArrayList<String>(Arrays.asList(line.split("\t")));
			while(str.size()<ncol+1){
				str.add("");
			}
			if(str.size() != ncol+1){
				throw new DataFormatException("StringMat: file format is wrong!");
			}
			rowname.add(str.get(0));
			for(i=1; i<str.size(); i++){
				String s = str.get(i);
				M[l-1][i-1] = s;
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
	public StringMat(StringMat m){
		nrow = m.nrow;
		ncol = m.ncol;
		colname2index = new HashMap<String, Integer>(m.colname2index);
		rowname2index = new HashMap<String, Integer>(m.rowname2index);
		colname = new ArrayList<String>(m.colname);
		rowname = new ArrayList<String>(m.rowname);
		M = new String[nrow][ncol];
		int i,j;
		for(i=0;i<nrow;i++){
			for(j=0;j<ncol;j++){
				M[i][j] = m.M[i][j];
			}	
		}
	}
	public StringMat(int i, int j){
		nrow = i;
		ncol = j;
		colname = new ArrayList<String>();
		rowname = new ArrayList<String>();
		colname2index = new HashMap<String, Integer>();
		rowname2index = new HashMap<String, Integer>();
		M = new String[i][j];
		for(i=0;i<nrow;i++){
			for(j=0;j<ncol;j++){
		      M[i][j] = "";
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
	public StringMat(List <String> row, List <String> col){
		nrow = row.size();
		ncol = col.size();
		colname = new ArrayList<String>();
		rowname = new ArrayList<String>();
		colname2index = new HashMap<String, Integer>();
		rowname2index = new HashMap<String, Integer>();
		M = new String[nrow][ncol];
		for(int i=0;i<nrow;i++){
			for(int j=0;j<ncol;j++){
		      M[i][j] = "";
		    }	
		}
		setColNames(col);
		setRowNames(row);		
	}
	public StringMat(String colname, Map<String, String> rowname2value){
		nrow = rowname2value.size();
		ncol = 1;
		this.colname = new ArrayList<String>();
		rowname = new ArrayList<String>();
		colname2index = new HashMap<String, Integer>();
		rowname2index = new HashMap<String, Integer>();
		M = new String[nrow][ncol];
		List <String> tmp = new ArrayList<String>();
		tmp.add(colname);
		setColNames(tmp);
		setRowNames(new ArrayList<String>(rowname2value.keySet()));		
		for(Map.Entry<String, String> e: rowname2value.entrySet()){
			set(e.getKey(), colname,  e.getValue());	
		}	
	}
	public static StringMat getStringMatFromTxtFile(String gmtFile) throws Exception{
			Map <String, String> annotMap = new HashMap<String, String>();	
		 BufferedReader inputStream = new BufferedReader(new FileReader(gmtFile));
			List<String> str = new ArrayList<String>();
			String line;	
			while((line = inputStream.readLine()) != null){
				if(line.charAt(0) == '#'){
					continue;				
				}
				str = Arrays.asList(line.split("\t"));
				String id = str.get(0);
				annotMap.put(id, "1");
			}
			inputStream.close();	
			if(annotMap.isEmpty()){
				throw new DataFormatException("readGeneSetFromTxTFile: file format is wrong!");
			}
		return new StringMat(gmtFile, annotMap);
	}
	
	 public static StringMat getStringMatFromGmtFile(String gmtFile) throws Exception{
		 Map <String, String> annotMap = new HashMap<String, String>();
		 BufferedReader inputStream = new BufferedReader(new FileReader(gmtFile));
		 List<String> str = new ArrayList<String>();
		 String line;
		 Set <String> seen = new HashSet<String>(); 
		 while((line = inputStream.readLine()) != null){
			 if(line.charAt(0) == '#'){
				 continue;
			 }
			 	str = Arrays.asList(line.split("\t"));
			 if(str.size() < 3 ){
				 	continue;
                 }
			 String annot = str.get(0);
			 str = MyFunc.uniq(str.subList(2, str.size()));
			 for(String s: str){
				 if(seen.contains(s)){
					 
				 }
				 seen.add(s);	
				 annotMap.put(s, annot);
              }	
		 }	
		 inputStream.close();
		 if(annotMap.isEmpty()){
			 throw new DataFormatException("readGeneSetFromGmtFile: file format is wrong!");
		 }
		 return new StringMat(gmtFile, annotMap);
	 }
	 
	 public static StringMat getStringMatFromGmtFile2(String gmtFile) throws Exception{
		 Map <String, List <String>> annotMap = new HashMap<String, List <String>>();
		 BufferedReader inputStream = new BufferedReader(new FileReader(gmtFile));
		 List<String> str = new ArrayList<String>();
		 String line;
		 Set <String> seen = new HashSet<String>();
		 while((line = inputStream.readLine()) != null){
			 if(line.charAt(0) == '#'){
				 continue;
			 }
			 str = Arrays.asList(line.split("\t"));
			 if(str.size() < 3 ){
				 	continue;
                 }
			 String annot = str.get(0);
			 str = MyFunc.uniq(str.subList(2, str.size()));	
			 annotMap.put(annot, str);
			 seen.addAll(str);
		 }	
		 inputStream.close();
		 if(annotMap.isEmpty()){
			 throw new DataFormatException("readGeneSetFromGmtFile: file format is wrong!");
		 }
		 StringMat tmp =  new StringMat(new ArrayList<String>(seen), new ArrayList<String>(annotMap.keySet()));
		 for(String s: annotMap.keySet()){
			 for(String t: annotMap.get(s)){
			 tmp.set(t, s, "1");
			 }
		 }
		 return tmp;
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
	public String get(int i, int j){
		if(i >= nrow || j >= ncol){
			throw new IndexOutOfBoundsException();
		}
		return M[i][j];
	}
	public String get(String i, String j){
		return M[rowname2index.get(i)][colname2index.get(j)];
	}
	public void set(int i, int j, String d){
		M[i][j] = d;
	}
	public void set(String i, String j, String d){
		M[rowname2index.get(i)][colname2index.get(j)] = d;	
	}
	public void setRowNames(List <String> v){
		rowname = new ArrayList<String>(v);
		rowname2index.clear();
		int i;
		for(i=0;i<nrow;i++){
			rowname2index.put(rowname.get(i),i);
		}
	}
	public void setColNames(List <String> v){
		colname = new ArrayList<String>(v);
		colname2index.clear();
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
		String[][] tmp = new String[ncol][nrow];
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
	public List<String> getRow(int i){
		List<String> tmp = new ArrayList<String>();
		int j;
		for(j = 0; j < ncol; j++){
			tmp.add(M[i][j]);
		}
		return tmp;
	}
	public List<String> getRow(String i){
		List<String> tmp = new ArrayList<String>();
		int j;
		for(j = 0; j < ncol; j++){
			tmp.add(M[rowname2index.get(i)][j]);
		}
		return tmp;
	}
	public Map<String,String> getRowMap(int i){
		Map<String, String> tmp = new HashMap<String, String>();
		int j;
		for(j = 0; j < ncol; j++){
			String s = M[i][j];
			tmp.put(colname.get(j), s);
		}
		return tmp;
	}
	public Map<String,String> getRowMap(String i){
		Map<String, String> tmp = new HashMap<String, String>();
		int j;
		for(j = 0; j < ncol; j++){
			String s = M[rowname2index.get(i)][j];
			tmp.put(colname.get(j),s);
		}
		return tmp;
	}
	public List<String> getCol(int j){
		List<String> tmp = new ArrayList<String>();
		int i;
		for(i = 0; i < nrow; i++){
			tmp.add(M[i][j]);
		}
		return tmp;
	}
	public List<String> getCol(String j){
		List<String> tmp = new ArrayList<String>();
		int i;
		for(i = 0; i < nrow; i++){
			tmp.add(M[i][colname2index.get(j)]);
		}
		return tmp;
	}
	public Map<String,String> getColMap(int j){
		Map<String, String> tmp = new HashMap<String, String>();
		int i;
		for(i = 0; i < nrow; i++){
			String s = M[i][j];
			tmp.put(rowname.get(i), s);
		}
		return tmp;
	}
	public Map<String,String> getColMap(String j){
		Map<String, String> tmp = new HashMap<String, String>();
		int i;
		for(i = 0; i < nrow;i++){
			String s = M[i][colname2index.get(j)];
			tmp.put(rowname.get(i),s);
		}
		return tmp;
	}
	
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
		String[][] new_M = new String[v.size()][ncol];
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
		String[][] new_M = new String[nrow][v.size()];
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
	
	public StringMat getSubMatrix(List row, List col){
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
		StringMat newMat = new StringMat(new_rowname, new_colname);
		for(i = 0; i < new_nrow; i++){
			for(j = 0; j < new_ncol; j++){
				 newMat.M[i][j] = M[v.get(i)][u.get(j)];
				
			}
		}
		return newMat;
	}
	
	public StringMat getSubMatByRow(List row){
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
		
		StringMat new_M = new StringMat(new_rowname, colname);
		for(String s: new_rowname){
			for(String t: colname){
				new_M.set(s, t, get(s,t));		
			}
		}
		return new_M;	
	}
	@SuppressWarnings("unchecked")
	public StringMat getSubMatByCol(List col){
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
		StringMat new_M = new StringMat(rowname, new_colname);
		for(String s: rowname){
			for(String t: new_colname){
				new_M.set(s, t, get(s,t));
			}
			
		}
		return new_M;		
	}
	public String toString(){
		StringBuffer S = new StringBuffer("\t"+  MyFunc.join("\t", getColNames()) + "\n");
		for(int i=0; i<nrow; i++){	
			S.append(getRowNames().get(i) + "\t"+  MyFunc.join("\t",MyFunc.toString(getRow(i))) + "\n");
		}
		return S.toString();
	}
	public boolean equals(StringMat m){
		if(colname != m.colname || rowname != m.rowname){
			return false;	
		}
		int i,j;
		for(i=0;i<nrow;i++){
			for(j=0;j<ncol;j++){
				if(M[i][j].equals(m.M[i][j])){
				  return false;
				}
			}
		}
		return true;
	}

	public StringMat bindRow(StringMat m){
		List <String> newRowname = new ArrayList<String>(rowname);
		newRowname.addAll(m.rowname);
		newRowname = MyFunc.uniq(newRowname);
		List <String> newColname = MyFunc.isect(colname, m.colname);
		StringMat tmp = new StringMat(newRowname, newColname);
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
	public StringMat addRow(String rowname, Map<String, String> value){
		List <String> newRowname = new ArrayList<String>(this.rowname);
		newRowname.add(rowname);
		StringMat tmp = new StringMat(newRowname, colname);
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
	public StringMat bindCol(StringMat m){
		List <String> newColname = new ArrayList<String>(colname);
		newColname.addAll(m.colname);
		newColname = MyFunc.uniq(newColname);
		List <String> newRowname = MyFunc.isect(rowname, m.rowname);
		StringMat tmp = new StringMat(newRowname, newColname);
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
	public StringMat addCol(String colname, Map<String,String> value){
		List <String> newColname = new ArrayList<String>(this.colname);
		newColname.add(colname);
		StringMat tmp = new StringMat(rowname, newColname);
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
	
	public StringMat bind(StringMat m){
		if(!MyFunc.isect(getRowNames(), m.getRowNames()).isEmpty()){
			return bindCol(m);
		}
		if(!MyFunc.isect(getRowNames(),m.getColNames()).isEmpty()){
			StringMat tmp = new StringMat(m);
			tmp.transpose();
			return bindCol(m);
		}
		if(!MyFunc.isect(getColNames(), m.getRowNames()).isEmpty()){
			StringMat tmp = new StringMat(m);
			tmp.transpose();
			return bindRow(tmp);
		}
		if(!MyFunc.isect(getColNames(),m.getColNames()).isEmpty()){
			return bindRow(m);
		}
		return this;		
	}
	public StringMat add(String name, Map<String, String> value){
		if(!MyFunc.isect(getRowNames(),new ArrayList<String>(value.keySet())).isEmpty()){
			return addCol(name, value);
		}
		if(!MyFunc.isect(getColNames(), new ArrayList<String>(value.keySet())).isEmpty()){
			return addRow(name, value);
		}
		return this;
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
	
	
}
	
	

