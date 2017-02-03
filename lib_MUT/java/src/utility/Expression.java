package utility;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.*;
import java.util.zip.DataFormatException;
import java.io.*;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;


public class Expression extends MyMat{
	private static final long serialVersionUID = 3504989313988715574L;
	Map<String, String> row_annot;
	public Expression(String infile)throws IOException, DataFormatException{
		super(infile);
		row_annot = new HashMap<String, String>();
	}
	public Expression(){};
	public Expression(MyMat E){
		super(E);
		row_annot = new HashMap<String, String>();
	}
	public Expression(Expression E){
		super(E);
		row_annot = new HashMap<String, String>(E.row_annot);
	}
	public Expression(List <String> genes, List <String> samples){
		super(genes, samples);
		row_annot = new HashMap<String, String>();
	}

	public static Expression getFromMySQL(String tableName){
		Connection con = MyFunc.getMySQLConnection("gene_expression");
		try {
			Statement stmt = con.createStatement();
			 ResultSet rs = stmt.executeQuery("describe " + tableName);
			 rs.next();
			 List <String> sampleNames = new ArrayList<String>();
			 while(rs.next()){
				 sampleNames.add(rs.getString(1));
			 }
			 rs = stmt.executeQuery("select * from " + tableName);
			 List <String> geneNames = new ArrayList<String>();
			 while(rs.next()){
				 geneNames.add(rs.getString(1));
			 }
			 Expression Exp =  new Expression(geneNames, sampleNames);
			 rs = stmt.executeQuery("select * from " + tableName);
			 while(rs.next()){
				 String gene = rs.getString(1);
				 for(String s: sampleNames){
					 Exp.set(gene, s, rs.getDouble(s));
				 }
			 }
			 stmt.close();
		     con.close();	
			 return Exp;	
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}
	
	
	public static Expression readGctFile(String infile) throws IOException, DataFormatException {
		Expression E = new Expression();
		E.colname = new ArrayList<String>();
		E.rowname = new ArrayList<String>();
		E.colname2index = new HashMap<String, Integer>();
		E.rowname2index = new HashMap<String, Integer>();
		E.row_annot = new HashMap<String, String>();
		BufferedReader inputStream = new BufferedReader(new FileReader(infile));
		List<String> lines = new ArrayList<String>();
		String line;
		E.nrow = -1;
		inputStream.readLine();
		inputStream.readLine();
		while((line = inputStream.readLine()) != null){
			if(line.charAt(0) == '#'){
				continue;				
			}
			lines.add(line);
			E.nrow++;
		}
		
		int i;	
		
		int l = 0;
		line = lines.get(0);
		List<String> str = Arrays.asList((line.split("\t")));
		E.ncol = str.size()-2;
		for(i = 2; i < str.size(); i++){
			if(E.colname.contains(str.get(i))){
				throw new DataFormatException("Expression.readFromGctFile: file format is wrong!");
			}
			E.colname.add(str.get(i));	
		}
			
		E.M = new double[E.nrow][E.ncol];
		for(l=1;l<=E.nrow;l++){
			line = lines.get(l);
			str = Arrays.asList(line.split("\t"));
			if(str.size() != E.ncol+2){
				throw new DataFormatException("Expression.readFromGctFile: file format is wrong!");
			}
			if(E.rowname.contains(str.get(0))){
				throw new DataFormatException("Expression.readFromGctFile: file format is wrong!");
			}
			E.rowname.add(str.get(0));
			E.row_annot.put(str.get(0), str.get(1));
			for(i=2; i<str.size(); i++){
				E.M[l-1][i-2] = Double.parseDouble(str.get(i));
			  }
		}
	   for(i=0;i<E.ncol;i++){
		   E.colname2index.put(E.colname.get(i),i);
	   }
	  for(i=0;i<E.nrow;i++){
		  E.rowname2index.put(E.rowname.get(i),i);
	  }
	  inputStream.close();
	  return E;
	}
	public static Expression readResFile(String infile) throws IOException, DataFormatException {
		Expression E = new Expression();
		E.colname = new ArrayList<String>();
		E.rowname = new ArrayList<String>();
		E.colname2index = new HashMap<String, Integer>();
		E.rowname2index = new HashMap<String, Integer>();
		E.row_annot = new HashMap<String, String>();
		BufferedReader inputStream = new BufferedReader(new FileReader(infile));
		List<String> lines = new ArrayList<String>();
		String line = inputStream.readLine();
		
		List<String> str = Arrays.asList((line.split("\t")));
		int i;
		for(i = 2; i < str.size(); i+=2){
			if(E.colname.contains(str.get(i))){
				throw new DataFormatException("Expression.readFromGctFile: file format is wrong!");
			}
			E.colname.add(str.get(i));	
		}
		E.ncol = E.colname.size();
		inputStream.readLine();
		inputStream.readLine();
		E.nrow = 0;
		while((line = inputStream.readLine()) != null){
			lines.add(line);
			E.nrow++;
		}
		E.M = new double[E.nrow][E.ncol];
		for(int l=0;l<E.nrow;l++){
			line = lines.get(l);
			str = Arrays.asList(line.split("\t"));
			if(str.size() != E.ncol*2+2){
				throw new DataFormatException("Expression.readFromGctFile: file format is wrong!");
			}
			if(E.rowname.contains(str.get(1))){
				throw new DataFormatException("Expression.readFromGctFile: file format is wrong!");
			}
			E.rowname.add(str.get(1));
			E.row_annot.put(str.get(0), str.get(0));
			int j = 0;
			for(i=2; i<str.size(); i+=2){
				E.M[l][j] = Double.parseDouble(str.get(i));
				j++;
			}
		}
	   for(i=0;i<E.ncol;i++){
		   E.colname2index.put(E.colname.get(i),i);
	   }
	  for(i=0;i<E.nrow;i++){
		  E.rowname2index.put(E.rowname.get(i),i);
	  }
	  inputStream.close();
	  return E;
	}
	public void printGctFile(){
		PrintWriter os = new PrintWriter(System.out);
		os.println("#1.2");
		os.println( nrow + "\t" + ncol);
		os.println("NAME\tDescription\t" +  MyFunc.join("\t", getColNames())); 
		int i,j;
		for(i=0; i<nrow; i++){
			List<Double> tmp = getRow(i);
			List<String> tmp2 = new ArrayList<String>();
			for(j=0;j<ncol;j++){
				tmp2.add( (tmp.get(j)).toString() );
			}
			if(row_annot.containsKey(rowname.get(i))){
				os.println(rowname.get(i) + "\t" + row_annot.get(rowname.get(i)) +  "\t"+  MyFunc.join("\t",tmp2));
			}else{
				os.println(rowname.get(i) + "\t\t"+  MyFunc.join("\t",tmp2));	
			}
		}	
		os.close();
	}
	public void printGctFile(String outfile) throws IOException {
		PrintWriter os = new PrintWriter(new FileWriter(outfile));
		os.println("#1.2");
		os.println( nrow + "\t" + ncol);
		os.println("NAME\tDescription\t" +  MyFunc.join("\t", getColNames())); 
		int i,j;
		for(i=0; i<nrow; i++){
			List<Double> tmp = getRow(i);
			List<String> tmp2 = new ArrayList<String>();
			for(j=0;j<ncol;j++){
				tmp2.add( (tmp.get(j)).toString() );
			}
			if(row_annot.containsKey(rowname.get(i))){
				os.println(rowname.get(i) + "\t" + row_annot.get(rowname.get(i)) +  "\t"+  MyFunc.join("\t",tmp2));
			}else{
				os.println(rowname.get(i) + "\t\t"+  MyFunc.join("\t",tmp2));	
			}
		}	
		os.close();
	}
	public void readRowAnnot(String infile) throws IOException, DataFormatException{
		BufferedReader inputStream = new BufferedReader(new FileReader(infile));
		List<String> str = new ArrayList<String>();
		row_annot = new HashMap<String, String>();
		String line;	
		while((line = inputStream.readLine()) != null){
			if(line.charAt(0) == '#'){
				continue;				
			}
			str = Arrays.asList(line.split("\t"));
			if(str.size() < 2){
				continue;
				//throw new DataFormatException("readRowAnnot: file format is wrong!");
			}
			if(rowname2index.containsKey(str.get(0))){
				row_annot.put(str.get(0), str.get(1));
			}
		}
		inputStream.close();
	}
	public void convertRowNames(){
		if(row_annot.isEmpty()){
			return;
		}
		Map<String,List<String>>annot2id = new TreeMap<String,List<String>>();
		for(Map.Entry<String,String> e : row_annot.entrySet()){
			if(!e.getValue().equals("") && !e.getValue().equals("---") && !e.getValue().equals("NA") && !e.getValue().equals("na") && !e.getValue().equals("NULL") && !e.getValue().equals("null")){
				if(annot2id.containsKey(e.getValue())){
					(annot2id.get(e.getValue())).add(e.getKey());
				}else{
					List<String> tmp = new ArrayList<String>();
					tmp.add(e.getKey());
					annot2id.put(e.getValue(), tmp);
				}
			}
		}
		List<String> selected_row = new ArrayList<String>();
		List<String> new_rowname = new ArrayList<String>();
		row_annot.clear();
		int i;
		for(Map.Entry<String,List<String>> e2 : annot2id.entrySet()){
			List <String> ids = e2.getValue();
			double max_v = MyFunc.var(getRow(ids.get(0)));
			String id = ids.get(0);
			 for(i = 1; i < ids.size(); i++){
			    double v = MyFunc.var(getRow(ids.get(i)));
			    if(v > max_v){
			    	max_v = v;
			    	id = ids.get(i);
			    }
			 }
			 selected_row.add(id);
			 new_rowname.add(e2.getKey());
			 row_annot.put(e2.getKey(), id);
		}
		
		reorderRows(selected_row);
		setRowNames(new_rowname);
		}
	public void convertRowNames(String infile)throws IOException, DataFormatException{	
		readRowAnnot(infile);
		convertRowNames();
		if(rowSize() == 0){
			throw new DataFormatException("convertRowNames: the chip file  is wrong!");	
		}
	}
	public void extractRows(String infile)throws IOException, DataFormatException{	
		BufferedReader inputStream = new BufferedReader(new FileReader(infile));
		List <String> rows  = new ArrayList<String>();
		String line;	
		while((line = inputStream.readLine()) != null){
			if(line.equals("")){
				continue;
			}
			if(line.charAt(0) == '#'){
				continue;				
			}
			List<String> str = Arrays.asList(line.split("\t"));
			rows.add(str.get(0));
		}
		inputStream.close();
		reorderRows(rows);
	}
	
	public void extractCols(String infile)throws IOException, DataFormatException{	
		BufferedReader inputStream = new BufferedReader(new FileReader(infile));
		List <String> rows  = new ArrayList<String>();
		String line;	
		while((line = inputStream.readLine()) != null){
			if(line.charAt(0) == '#'){
				continue;				
			}
			List<String> str = Arrays.asList(line.split("\t"));
			rows.add(str.get(0));
		}
		inputStream.close();
		reorderCols(rows);
	}
	
	
	public boolean hasRowAnnot(){
		if(row_annot.isEmpty()){
			return false;
		}else{
			return true;
		}
	}
	
	public void ceil(double d){
		int i,j;
		for(i=0;i<nrow;i++){
		    for(j=0;j<ncol;j++){
		      if(get(i,j) > d){
		    	  set(i, j, d);
		      }  
		    }
		}
	}

	public void floor(double d){
		int i,j;
		for(i=0;i<nrow;i++){
			for(j=0;j<ncol;j++){
			  if(get(i,j) < d){
				  set(i, j, d);
			  }  
			}
		}
	}
	
	public void takeLog(double base){
		  int i,j;
		  double b = Math.log(base);
		  for(i=0;i<nrow;i++){
		    for(j=0;j<ncol;j++){
		     double d = get(i,j);
		     if(d <= 0){
		    	 throw  new ArithmeticException();
		     }
		      set(i, j,  Math.log(d)/b);
		    }
		  }
	}
	public void adjustColMean(double m){
		  int i,j;
		  List <Double> mean = getRowMeans();
		  for(j=0;j<ncol;j++){
			  for(i=0;i<nrow;i++){
		      set(i, j,  get(i,j)*m/mean.get(j));
		    }
		  }
	}
	
	
	
	public void sampleRows(double r){
		int n;
		if(r > 1){
			n = (int) Math.round(r);
		}else{
			n =   (int) Math.round(this.rowSize()*r);
		}
		reorderRows(MyFunc.sample(getRowNames(), n));
	}
	
	public void sampleColumns(double r){
		int n;
		if(r > 1){
			n = (int) Math.round(r);
		}else{
			n =   (int) Math.round(this.colSize()*r);
		}
		reorderCols(MyFunc.sample(getColNames(), n));
	}
	
	public Expression getModuleProfile( Map <String, List<String> > geneSets){
		Map <String, List<Double> > profile = new LinkedHashMap<String, List<Double>>();
		for(String s: geneSets.keySet()){
			List <String>  genes = MyFunc.isect(geneSets.get(s), getRowNames());
			if(genes.isEmpty()){
				continue;
			}
			List <Double> v = new ArrayList<Double>();
			for(int i = 0; i < ncol; i++){
				v.add(0.0);
			}
			for(String g: genes){ 
				List <Double> tmp = getRow(g);
				tmp = MyFunc.normalize(tmp);
				for(int i = 0; i < ncol; i++){
					v.set(i, v.get(i)+tmp.get(i));	
			    }
			}
			v = MyFunc.normalize(v);
			profile.put(s, v);			
		}
		Expression E  = new Expression(new ArrayList<String>(profile.keySet()), getColNames());
		for(String s: E.rowname){
			for(String t: E.colname){
				E.set(s, t, profile.get(s).get(E.colname2index.get(t)));
			}
		}
		return E;
	}
	
	
	public Expression getRowMeansAsExpression(){
		List<Double> d = getRowMeans();
		List <String> tmp = new ArrayList<String>();
		tmp.add("mean");
		Expression E = new Expression(tmp, getColNames());
		for(int j = 0; j < d.size(); j++){
			E.set(0,j,d.get(j));
		}
		return E;
	}
	
	
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("c", "ncol", false,  "normalize columns");
		options.addOption("r", "nrow", false,  "normalize rows");
		options.addOption("g", "readgct", false,  "read gct file");
		options.addOption("F", "readres", false,  "read res file");
		options.addOption("G", "writegct", false,  "write gct file");
		options.addOption("l", "log", true, "take log (base value)");
		options.addOption("C", "ceil", true, "ceil (cutoff value)");
		options.addOption("f", "foor", true, "floor (cutoff value)");
		options.addOption("a", "annot", true, "add chip annotation (chip file)");
		options.addOption("R", "convrow", false, "convert rowname");
		options.addOption("v", "varfilt", true, "variation filter (the number or rate of filtered rows)");
		options.addOption("s", "shufrow", false, "shuffle rows");
		options.addOption("t", "transpose", false, "transpose");
		options.addOption("e", "exractrow", true, "extract rows (a tab-delimited file containing row names in the 1st column)");
		options.addOption("E", "exractcol", true, "extract columns (a tab-delimited file containing column names in the 1st column)");
		options.addOption("b", "bind", false, "bind expression matrices");
		options.addOption("S", "smpr", true, "sample rows (the number or rate of sampled rows)");
		options.addOption("T", "smpc", true, "sample columns (the number or rate of sampled columns)");
		options.addOption("m", "colmean", true, "adjust col means (mean value)");
		options.addOption("p", "prop", false,  "show property");
		options.addOption("M", "mod", true, "get module profiles (gmt file)");
		options.addOption("N", "mean", false, "take row means");
		options.addOption("B", "bin", true, "binarize");
		options.addOption("L", "binless", true, "binarize so that the smaller is 1");
		options.addOption("H", "binper", true, "binarize (cutoff by percentile)");
		options.addOption("K", "binlessper", true, "binarize so that the smaller is 1 (cutoff by percentile)");
		options.addOption("A", "nall", false, "normalize all");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp("Expression [options] expfile", options);
			return;
		}
		List <String> argList = commandLine.getArgList();
		if(!commandLine.hasOption("b") && argList.size() != 1){
			formatter.printHelp("Expression [options] expfile ", options);
			return;
		}
		
		Expression E;
		if(commandLine.hasOption("g")){
			E = readGctFile(argList.get(0));
		}else if(commandLine.hasOption("F")){
			E = readResFile(argList.get(0));
		}else {
			E = new Expression(argList.get(0));
		}
		if(commandLine.hasOption("b")){
			for(int i = 1; i < argList.size();i++){
				if(commandLine.hasOption("g")){
					E = new Expression (E.bind(readGctFile(argList.get(i))));
				}else{
					E = new Expression(E.bind(new Expression(argList.get(i))));
				}
			}
		}
	
		if(commandLine.hasOption("C")){
			E.ceil(Double.valueOf(commandLine.getOptionValue("C")));
		}
		if(commandLine.hasOption("f")){
			E.floor(Double.valueOf(commandLine.getOptionValue("f")));
		}

		if(commandLine.hasOption("l")){
			String s = commandLine.getOptionValue("l");
			if(s.toLowerCase().equals("e")){
				E.takeLog(Double.valueOf(Math.E));
			}else{
				E.takeLog(Double.valueOf(s));
			}
		}
		
		if(commandLine.hasOption("a")){ 
			E.readRowAnnot(commandLine.getOptionValue("a"));
		}	
		if(commandLine.hasOption("R")){
			if(E.hasRowAnnot()){
				E.convertRowNames();
			}else{
				System.err.println("specify an annot file !");
				return;
			}
		}		
		if(commandLine.hasOption("c")){
			E.normalizeCols();
		}
		if(commandLine.hasOption("r")){
			E.normalizeRows();
		}
		if(commandLine.hasOption("A")){
			E.normalizeAll();
		}
		
		if(commandLine.hasOption("v")){
			E.filterRowByVariance(Double.valueOf(commandLine.getOptionValue("v")));
		}
		if(commandLine.hasOption("s")){
			E.shuffleRows();
		}
		if(commandLine.hasOption("t")){
			E.transpose();
		}
		if(commandLine.hasOption("e")){
			E.extractRows(commandLine.getOptionValue("e"));
		}
		if(commandLine.hasOption("E")){
			E.extractCols(commandLine.getOptionValue("E"));
		}
		if(commandLine.hasOption("S")){
			E.sampleRows(Double.valueOf(commandLine.getOptionValue("S")));
		}
		if(commandLine.hasOption("T")){
			E.sampleColumns(Double.valueOf(commandLine.getOptionValue("T")));
		}
		if(commandLine.hasOption("m")){
			E.adjustColMean(Double.valueOf(commandLine.getOptionValue("m")));
		}
		if(commandLine.hasOption("B") & commandLine.hasOption("L")){
			 E.ternarize(Double.valueOf(commandLine.getOptionValue("B")), Double.valueOf(commandLine.getOptionValue("L")));
		}else{
			if(commandLine.hasOption("B")){
				E.binarize(Double.valueOf(commandLine.getOptionValue("B")));
			}
			if(commandLine.hasOption("L")){
				E.binarizeLess(Double.valueOf(commandLine.getOptionValue("L")));
			}
		}
		
		if(commandLine.hasOption("H") & commandLine.hasOption("K")){
			E.ternarize(E.percentile(Double.valueOf(commandLine.getOptionValue("H"))), 
					E.percentile(Double.valueOf(commandLine.getOptionValue("K"))));
		}else{
			if(commandLine.hasOption("H")){
				E.binarize(E.percentile(Double.valueOf(commandLine.getOptionValue("H"))));
			}
			if(commandLine.hasOption("K")){
				E.binarizeLess(E.percentile(Double.valueOf(commandLine.getOptionValue("K"))));
			}
		}
		
		
		if(commandLine.hasOption("N")){
			E = E.getRowMeansAsExpression();
		}
		if(commandLine.hasOption("p")){
			E.printDimension();
			E.printStatistics();
			return;
		}
		if(commandLine.hasOption("M")){
			E = E.getModuleProfile(MyFunc.readGeneSetFromGmtFile(commandLine.getOptionValue("M")));
		}
		if(commandLine.hasOption("G")){
			E.printGctFile();
		}else{
			E.print();
		}
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
}
