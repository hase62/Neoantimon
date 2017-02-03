package utility;
import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.util.zip.DataFormatException;


public class CoxRegression {
	private MyMat Exp;
	private Map <String, Double> sample2time;
	private Map <String, Boolean> sample2status; //true: dead, false: sensored
	private Set <String> samples;
	private List <String> genes;
	private String tmpFile;
	private Map <String, Double> Zscore;
	private Map <String, Double> Pvalue;
	
	public CoxRegression(MyMat Exp){
		this.Exp = Exp;
		samples = new HashSet<String>(Exp.colname);
		genes = new ArrayList<String>(Exp.rowname);
		sample2time = new HashMap<String, Double>();
		sample2status = new HashMap<String, Boolean>();
	}
	public void setTime(Map <String, String> time){
		for(Map.Entry<String, String> e: time.entrySet()){
			if(samples.contains(e.getKey()) && !e.getValue().equals("") && MyFunc.canBeDouble(e.getValue())){
				sample2time.put(e.getKey(), Double.valueOf(e.getValue()));
			}
		}
		samples = sample2time.keySet();
	}
	public void setStatus(Map <String, String> status){
		for(Map.Entry<String, String> e: status.entrySet()){
			if(samples.contains(e.getKey())){
				if(e.getValue().equals("1")|| e.getValue().toUpperCase().equals("TRUE")){
					sample2status.put(e.getKey(), true);
				}
				if(e.getValue().equals("0")|| e.getValue().toUpperCase().equals("FALSE")){
					sample2status.put(e.getKey(), false);
				}
			}
		}
		samples = sample2status.keySet();
	}
	
	
	private void printInputFile() throws IOException{
		tmpFile = "tmp" + hashCode();
		PrintWriter os = new PrintWriter(new FileWriter(tmpFile + ".in"));
		for(String s: samples){
			os.println(s + "\t" + sample2time.get(s) + "\t" + (sample2status.get(s)==true?1:0));
		}
		os.close();
		Exp.print(tmpFile + ".exp");
	}
	private void printScriptFile()throws IOException{
		PrintWriter os = new PrintWriter(new FileWriter(tmpFile + ".R"));
		os.println("library(survival)");
		os.println("A<-read.table('"  + tmpFile  +".in')"); 
		os.println("rownames(A)<-A[,1]");
		os.println("A<-A[,2:3]");
		os.println("colnames(A)<-c('time', 'status')"); 
		os.println("E<-as.matrix(read.table(\"" +  tmpFile + ".exp" +  "\",row.names=1,sep=\"\\t\",quote=\"\",header=T))");
		os.println("E<-t(E)");
		os.println("E<-scale(E)");
		os.println("E<-E[rownames(A),]");
		os.println("m<-colnames(E)");
		os.println("X<-cbind(as.data.frame(E),as.data.frame(A))");
		os.println("CR<-NULL");
		os.println("for(i in 1:ncol(E)){");
		os.println("x<-c(X[\"time\"], X[\"status\"], X[m[i]])");
		os.println("names(x)[3]<-\"m\"");
		os.println(" CR[[i]] <- coxph(Surv(time, status) ~ m, data=x)");
		os.println("}");
		os.println("p<-NULL");
		os.println("z<-NULL");
		os.println("for(i in 1:ncol(E)){");
		os.println("z<-c(z, (summary(CR[[i]])$coef)[4])");
		os.println("p<-c(p, (summary(CR[[i]])$coef)[5])");
		os.println("}");
		os.println("outfile<-paste(c('" + tmpFile + "', '.p'), collapse='')"); 
		os.println("write(p, outfile ,ncolumns=1) "); 
		os.println("outfile<-paste(c('" + tmpFile + "', '.z'), collapse='')");
		os.println("write(z, outfile ,ncolumns=1) ");
		os.close();
	}
	
	private  void runRscript() throws Exception{
		MyFunc.runRscript(tmpFile + ".R");
	}
	
	private  void readResultFiles() throws IOException{
	   Pvalue = new HashMap<String, Double>();
	   Zscore = new HashMap<String, Double>();
	   List <Double> tmp = MyFunc.readDoubleList(tmpFile + ".p");
	   for(int i = 0, n = genes.size(); i < n; i++){
		   Pvalue.put(genes.get(i), tmp.get(i));	   
		}
	   tmp = MyFunc.readDoubleList(tmpFile + ".z");
	   for(int i = 0, n = genes.size(); i < n; i++){
		   Zscore.put(genes.get(i), tmp.get(i));	   
		}
	}
	private  void removeTmpFiles() throws Exception{
		List <String> files = new ArrayList<String>();
		files.add(tmpFile + ".in");
		files.add(tmpFile + ".exp");
		files.add(tmpFile + ".R");
		files.add(tmpFile + ".p");
		files.add(tmpFile + ".z");
		MyFunc.removeFiles(files);
	}
	
	public void perform() throws Exception{
		if(sample2status.isEmpty() || sample2time.isEmpty()){
			throw new MyException("Status and time must be set!"); 
		}
		printInputFile();
		printScriptFile();
		runRscript();
		readResultFiles();
		removeTmpFiles();
	}
	
	public Map<String, Double> getPvalue(){
		return Pvalue;
	}
	public Map<String, Double> getZscore(){
		return Zscore;
	}	
}
