package utility;
import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.util.zip.DataFormatException;


public class KaplanMeierAnalysis {
	private Map <String, Double> sample2time;
	private Map <String, Boolean> sample2status; //true: dead, false: sensored
	private Map <String, String> sample2group;
	private Set <String> samples;
	private List <String> groups;
	private String tmpFile;
	private Map <String, List <Double[]>> survivalCurves;
	private Double Pvalue; 
	
	public KaplanMeierAnalysis(Map <String, String> sample2group){
		this.sample2group = new HashMap<String,String>();
		for(Map.Entry<String, String> e: sample2group.entrySet()){
			if(!e.getValue().equals("")){
				this.sample2group.put(e.getKey(), e.getValue());
			}
		}
		samples = new HashSet<String>(this.sample2group.keySet());
		sample2time = new HashMap<String, Double>();
		sample2status = new HashMap<String, Boolean>();
		groups = new ArrayList<String>(new TreeSet<String>(this.sample2group.values()));
	}
	
	public void setTime(Map <String, String> time){
		for(Map.Entry<String, String> e: time.entrySet()){
			if(samples.contains(e.getKey()) && !e.getValue().equals("") && MyFunc.canBeDouble(e.getValue())){
				sample2time.put(e.getKey(), Double.valueOf(e.getValue()));
			}
		}
		samples.retainAll(sample2time.keySet());
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
		samples.retainAll(sample2status.keySet());
	}
	
	private void printInputFile() throws IOException{
		tmpFile = "tmp" + hashCode();
		PrintWriter os = new PrintWriter(new FileWriter(tmpFile + ".in"));
		for(String s: samples){
			os.println(s + "\t" + (groups.indexOf(sample2group.get(s))+1) + "\t" + sample2time.get(s) + "\t" + (sample2status.get(s)==true?1:0));
		}
		os.close();
	}
		
	private  void printScriptFile()throws IOException{
		PrintWriter os = new PrintWriter(new FileWriter(tmpFile + ".R"));
		os.println("library(survival)");
		os.println("A<-read.table('"  + tmpFile  +".in')"); 
		os.println("rownames(A)<-A[,1]");
		os.println("A<-A[,2:4]");
		os.println("colnames(A)<-c('group', 'time', 'status')"); 
		os.println("Sv<-survfit(Surv(time, status)~group, data=A)"); 
		os.println("n<-" + groups.size()); 
		os.println("k<-1"); 
		os.println("for(i in 1:n){"); 
		os.println("time<-Sv$time[k:(k+Sv$ntimes.strata[i]-1)]"); 
		os.println("surv<-Sv$surv[k:(k+Sv$ntimes.strata[i]-1)]"); 
		os.println("k<-k+Sv$ntimes.strata[i]"); 
		os.println("outfile<-paste(c('" + tmpFile + "', '.time', i), collapse='')"); 
		os.println("write(time, outfile ,ncolumns=1) "); 
		os.println("outfile<-paste(c('" + tmpFile + "', '.surv', i), collapse='')");
		os.println("write(surv, outfile ,ncolumns=1) ");
		os.println("}"); 
		os.println("p<-pchisq(survdiff(Surv(time, status)~group, data=A)$chisq,"+ (groups.size()-1 )+",lower.tail=FALSE)"); 
		os.println("outfile<-paste(c('" + tmpFile + "', '.pvalue'), collapse='')");
		os.println("write(p, outfile ,ncolumns=1) ");
		os.close();
	}
	
	private  void readResultFiles() throws IOException{
		survivalCurves = new HashMap<String, List<Double[]>>();
		for(int i = 0, n = groups.size(); i < n; i++){
			List <Double> t = MyFunc.readDoubleList(tmpFile + ".time" + (i+1));
			List <Double> s = MyFunc.readDoubleList(tmpFile + ".surv" + (i+1));
 			List <Double[]> tmp = new ArrayList<Double[]>();
 			for(int j = 0, m = t.size(); j < m; j++){
 				Double[] tmp2 = {t.get(j), s.get(j)}; 
 				tmp.add(tmp2);
 			}
 			
			survivalCurves.put(groups.get(i), tmp);
		}
		BufferedReader inputStream = new BufferedReader(new FileReader(tmpFile + ".pvalue"));
		Pvalue = Double.valueOf(inputStream.readLine());
		inputStream.close();
	}
	
	private  void runRscript() throws Exception{
		MyFunc.runRscript(tmpFile + ".R");
	}
	private  void removeTmpFiles() throws Exception{
		List <String> files = new ArrayList<String>();
		files.add(tmpFile + ".in");
		files.add(tmpFile + ".R");
		files.add(tmpFile + ".pvalue");
		for(int i = 0; i < groups.size(); i++){
			files.add(tmpFile + ".time" + (i+1));
			files.add(tmpFile + ".surv" + (i+1));
		}
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
	
	public SurvivalCurvePlotViewer getSurvivalCorvePlotViewer(){
		SurvivalCurvePlotViewer viewer = new SurvivalCurvePlotViewer(survivalCurves);
		viewer.setPvalue(Pvalue);
		viewer.setColorOfYautomatically();
		return viewer;
	}
	
	
	
	
	
}
