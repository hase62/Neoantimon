package utility;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

import sun.reflect.Reflection;

public class MatrixCorrelation {

	private MyMat M1;
	private MyMat M2;
	
	private MyMat Pearson;
	private MyMat FisherZ;
	
	private Integer itrForNullDist = 100000;
	
	private MyMat Pvalue;
	private MyMat Qvalue;
	
	private Double PvalueCutoff = null;
	private Double QvalueCutoff = null;	
	private boolean minusLogScale = false;
	
	private boolean useExpReg = false; 
	private boolean supressPcalculation  = false;
	
	
	public MatrixCorrelation(MyMat m1, MyMat m2) {
		M1 = new MyMat(m1);
		M2 = new MyMat(m2);
		if(MyFunc.isect(m1.getRowNames(), m2.getColNames()).isEmpty()){
		}
		else if(!MyFunc.isect(m1.getRowNames(), m2.getColNames()).isEmpty()){	
			M1.transpose();
		}		
		else if(!MyFunc.isect(m2.getRowNames(), m1.getColNames()).isEmpty()){	
			M2.transpose();
		}
		else if(!MyFunc.isect(m1.getRowNames(), m2.getRowNames()).isEmpty()){	
			M1.transpose();
			M2.transpose();
		}
		else {
			throw new MyException("expressionProfile and sampleLabels must have common rows or columns!");	
		}
		M1.normalizeRows();
		M2.normalizeRows();
		Pearson = new MyMat(M1.getRowNames(), M2.getRowNames());
		FisherZ = new MyMat(M1.getRowNames(), M2.getRowNames());
		Pvalue = new MyMat(M1.getRowNames(), M2.getRowNames());
		Qvalue = new MyMat(M1.getRowNames(), M2.getRowNames());
		
	}	

	public MatrixCorrelation(MyMat m1){
		M1 = new MyMat(m1);
		M1.normalizeRows();
		Pearson = new MyMat(M1.getRowNames(), M1.getRowNames());
		FisherZ = new MyMat(M1.getRowNames(), M1.getRowNames());
		Pvalue = new MyMat(M1.getRowNames(), M1.getRowNames());
		Qvalue = new MyMat(M1.getRowNames(), M1.getRowNames());
	}
	
	public void setItrForNullDist(int i){
		itrForNullDist = i;
	}
	
	public void calculate(){
		calculateCorrelation();
		if(supressPcalculation == false){
			calculatePvalue();
			calculateQvalue();
		}
	}
	
	public void useExpRegression(){
		useExpReg = true;
	}
	
	public void calculateCorrelation(){
		for(String s: Pearson.rowname){	
			for(String t: Pearson.colname){
				double c;
				if(M2 == null){
					if(s.equals(t)){
						Pearson.set(s, t, 1);
						FisherZ.set(s, t, Double.MAX_VALUE);
						continue;
					}
					if(s.compareTo(t) < 0){
						continue;
					}
					c = MyFunc.pearsonCorrelationForNormarizedList(M1.getRow(s), M1.getRow(t));	
				}else{
					c = MyFunc.pearsonCorrelationForNormarizedList(M1.getRow(s), M2.getRow(t));
				}
				double z;
				Pearson.set(s, t, c);
				if(c >= 1){
					z = Double.MAX_VALUE;
				}else if(c <= -1){
					z = -Double.MAX_VALUE;
				}else{	
					z  = 0.5*Math.log((1+c)/(1-c));
				}
				Pearson.set(s,t, c);
				FisherZ.set(s,t, z);
				if(M2 == null){
					Pearson.set(t,s, c);
					FisherZ.set(t,s, z);
				}
			}
		}
	}
	
	public void calculatePvalue(){
		for(String s: Pearson.rowname){	
		
		List <Double> nullDist = new ArrayList<Double>();
		while(nullDist.size() <itrForNullDist){
				MyMat randomM = new MyMat((M2 != null)?M2:M1);
				randomM.shuffleCols();
				for(String t: randomM.rowname){
					double c =  Math.abs(MyFunc.pearsonCorrelationForNormarizedList(M1.getRow(s), randomM.getRow(t)));	
					if(c >= 1){  
						nullDist.add(Double.MAX_VALUE);
					}else{	
						nullDist.add(0.5*Math.log((1+c)/(1-c)));
					}
					
			}
		}
		Collections.sort(nullDist);
		Collections.reverse(nullDist);
		if(useExpReg){
			MyFunc.Density D = new MyFunc.Density(nullDist); 
			double minX = MyFunc.percentile(nullDist, 0.6);
			double maxX = MyFunc.percentile(nullDist, 0.95);
			double d = (maxX - minX) / 500;
			List <Double> X = new ArrayList<Double>();
			double x;
			for(x = minX; x <= maxX; x += d ){
				X.add(x);
			}
			List <Double> Y = D.estimate(X);
			
			ExpRegression ER  = new ExpRegression();
			ER.setTrainingData(X, Y);
			List <Double> powers = new ArrayList<Double>();
			double minPower = 1.0;
			double maxPower = 3.0;
			for( double p = minPower; p <= maxPower; p += 0.1){
				powers.add(p);
			}
			try {
				ER.optimizePowerOfX(powers);
			} catch (MyException e) {
				e.printStackTrace();
			}
			
			/*for(double p = 0.9; p < 0.990001;  p += 0.01){
				double  t   = MyFunc.percentile(nullDist, p);
				double P = ER.integrate(t, upperLimit);
				System.err.println(1-p + "\t" + P + "\t" + t + "\t" + upperLimit);
			}*/
			
			double i;
			double upperLimit = MyFunc.mean(nullDist)*50;
			for(String t: Pearson.colname){	
				if(M2 == null){
					if(s.equals(t)){
						Pvalue.set(s, t, 0);
					}
					if(s.compareTo(t) < 0){
						continue;
					}
				}
				for(i=0;i<nullDist.size() && nullDist.get((int)i) > Math.abs(FisherZ.get(s, t));i++){}
				double P = i / nullDist.size();
				if(P < 0.05){
					if( Math.abs(FisherZ.get(s, t)) >= upperLimit){
						P = Double.MIN_VALUE;
					}else{
						P = ER.integrate(Math.abs(FisherZ.get(s, t)), upperLimit);
					}
				}
				Pvalue.set(s,t, P);
				if(M2 == null){
					Pvalue.set(t,s, P);
				}
			}
			
		}else{
			double i;
				for(String t: Pearson.colname){
					if(M2 == null){
						if(s.equals(t)){
							Pvalue.set(s, t, 0);
						}
						if(s.compareTo(t) < 0){
							continue;
						}
					}
					for(i=0;i<nullDist.size() && nullDist.get((int)i) > Math.abs(FisherZ.get(s, t));i++){}
					if(i==0){
						i=1;
					}
					double P = i / nullDist.size();
					Pvalue.set(s,t, P);
					if(M2 == null){
						Pvalue.set(t,s, P);
					}
				}
			
		}
		}
	}
		
		
	public void setPvalueCutoff(double d){
		PvalueCutoff = d;	
	}
	public void setQvalueCutoff(double d){
		QvalueCutoff = d;	
	}
	public void outputInMinusLogScale(){
		minusLogScale = true;
	}
	
	
	public double getPvalue(String ID1, String ID2){
		return Pvalue.get(ID1, ID2);	
	}
	
	public double getQvalue(String ID1, String ID2){
		return Qvalue.get(ID1, ID2);	
	}
	
	public Map<String, Double> getPvalueMap(){
		return  asMap(Pvalue);
	}
	public Map<String, Double> getQvalueMap(){
		return asMap( Qvalue);
	}
	public Map<String, Double> getPearsonMap(){
		return  asMap(Pearson);
	}
	public Map<String, Double> getFisherZMap(){
		return asMap(FisherZ);
	}
	
	private  Map <String, Double> asMap(MyMat M){
		Map <String, Double> m = new HashMap<String, Double>();
		if(M2 != null){
		for(String s: M.rowname){	
			for(String t: M.colname){
				m.put(s + "\t" + t,M.get(s,t));	
			}
		}
		}else{
			List <String> tmp = M.getRowNames();
			for(String s: tmp){
				for(String t: tmp){
					if(s.compareTo(t)> 0){
						m.put(s + "\t" + t,M.get(s,t));
					}		
				}
			}
		}
		return m;		
	}
	
	
	public void calculateQvalue(){
		Map <String, Double> PvalueMap = Pvalue.asMap();
		
		Map <String, Double> QvalueMap = MyFunc.calculateStoreyQvalue(PvalueMap);
		for(Map.Entry<String, Double> e: QvalueMap.entrySet()){
			List <String> tmp = MyFunc.split("\t", e.getKey());
			Qvalue.set(tmp.get(0), tmp.get(1), e.getValue());
		
			if(M2 == null ){
				Qvalue.set(tmp.get(1), tmp.get(0), e.getValue());
			}	
		}
		if(M2 == null ){
			List <String> tmp = new ArrayList<String>(M1.getRowNames());
			for(int i = 0, n = tmp.size(); i <n ;i++){
				Qvalue.set(tmp.get(i), tmp.get(i), Double.MIN_VALUE);	
			}				
		}
	}
	
	
	public String toString(){
		StringBuffer S = new StringBuffer("\t\tP value\tQ value\tPearson\tFisherZ\n");
		Map <String, Double> PvalueMap = getPvalueMap();
		Map <String, Double> QvalueMap = getQvalueMap();
		Map <String, Double> PearsonMap = getPearsonMap();
		Map <String, Double> FisherZMap = getFisherZMap();
		List<String> keys =  MyFunc.sortKeysByAscendingOrderOfValues(PvalueMap);
		for(String s: keys){
			List <String>  tmp = new ArrayList<String>();
			tmp.add(s);
			double p; 
			if(minusLogScale){
				p = (PvalueMap.get(s)==0)?Double.MAX_VALUE: -Math.log10(PvalueMap.get(s));
			}else{
				p = PvalueMap.get(s);
			}
			if(PvalueCutoff != null && (minusLogScale?(p < PvalueCutoff):(p > PvalueCutoff))){
				continue;
			}
			double q; 
			if(minusLogScale){
				q = (QvalueMap.get(s)==0)?Double.MAX_VALUE: -Math.log10(QvalueMap.get(s));
			}else{
				q = QvalueMap.get(s);
			}
			if(QvalueCutoff != null &&  (minusLogScale?(q < QvalueCutoff):(q > QvalueCutoff))){
				continue;
			}
			tmp.add(Double.toString(p));
			tmp.add(Double.toString(q));
			tmp.add(Double.toString(PearsonMap.get(s)));
			tmp.add(Double.toString(FisherZMap.get(s)));
			S.append(MyFunc.join("\t", tmp) + "\n");
		}
		return S.toString();
	}
	
	public MyMat getPvalueMatrix(){
		if(minusLogScale){
			MyMat tmp = new MyMat(Pvalue.getRowNames(), Pvalue.getColNames());
			for(String s: Pvalue.getRowNames()){
				for(String t: Pvalue.getColNames()){
					tmp.set(s, t,(Pvalue.get(s, t)==0)?Double.MAX_VALUE:-Math.log10(Pvalue.get(s, t)));
				}
			}
			return tmp;
		}else{
			return Pvalue;
		}
	}
	
	public MyMat getQvalueMatrix(){
		if(minusLogScale){
			MyMat tmp = new MyMat(Qvalue.getRowNames(), Qvalue.getColNames());
			for(String s: Qvalue.getRowNames()){
				for(String t: Qvalue.getColNames()){
					tmp.set(s, t,(Qvalue.get(s, t)==0)?Double.MAX_VALUE:-Math.log10(Qvalue.get(s, t)));
				}
			}
			return tmp;			
		}else{
			return Qvalue;
		}
	}
	
	public MyMat getPearsonMatrix(){
		return Pearson;
	}
	
	public MyMat getFisherZMatrix(){
		return FisherZ;
	}
	
	
	public static void main(String[] args) throws Exception {
		Options options = new Options();
		options.addOption("l", "minuslog", false, "output in minus log scale");
		options.addOption("p", "pcutoff", true, "overlap p-value cutoff");
		options.addOption("q", "pcutoff", true, "overlap q-value cutoff");
		options.addOption("P", "pmat", false, "get p-value matrix");
		options.addOption("Q", "qmat", false, "get q-value matrix");
		options.addOption("c", "cmat", false, "get correlation score matrix");
		options.addOption("i", "itrnull", true, "# of iteration for null dist genetation");
		options.addOption("r", "regress", false, "use regression for p-value calculation");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		
		
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] gmt_file gmt_file", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 2 && argList.size() != 1){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] gmt_file gmt_file", options);
			return;
		}
		
		MatrixCorrelation MC;
		if(argList.size() != 1){
			MC = new MatrixCorrelation(new MyMat(argList.get(0)), new MyMat(argList.get(1)));
		}else{
			MC = new MatrixCorrelation(new MyMat(argList.get(0)));
		}
		if(commandLine.hasOption("l")){
			MC.outputInMinusLogScale();
		}
		if(commandLine.hasOption("i")){
			MC.setItrForNullDist(Integer.valueOf(commandLine.getOptionValue("i")));
		}
		if(commandLine.hasOption("p")){
			MC.setPvalueCutoff(Double.valueOf(commandLine.getOptionValue("p")));
		}
		if(commandLine.hasOption("q")){
			MC.setQvalueCutoff(Double.valueOf(commandLine.getOptionValue("q")));
		}
		if(commandLine.hasOption("r")){
			MC.useExpRegression();
		}
		if(commandLine.hasOption("c")){
			MC.supressPcalculation = true;
			MC.calculate();
			System.out.print(MC.getPearsonMatrix());
		}else{
			MC.calculate();
		if(commandLine.hasOption("P")){
			System.out.print(MC.getPvalueMatrix());
		}else{
			if(commandLine.hasOption("Q")){
				System.out.print(MC.getQvalueMatrix());
			}else{
				System.out.print(MC);
				}
			}
		}				
	
	}
	
	
	
}
