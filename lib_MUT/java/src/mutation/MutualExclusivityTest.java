package mutation;

import java.util.*;
import utility.*;
import org.apache.commons.cli.*;
import sun.reflect.Reflection;

public class MutualExclusivityTest {
	
	MyMat M;
	
	double stat;
	double pval;
	int perm = 1000;
	
	public MutualExclusivityTest(MyMat M){
		this.M = M;
	}
	
	private double calculateStatistic(MyMat M){
		double s = 0;
		for (int j = 0; j < M.colSize(); j++){
			List <Double> tmp = M.getCol(j);
			double t = 0;
			for(Double d: tmp){
				if(d>0.01){
					t++;
				}
			}
			if(t > 1){
				s++;
			}
		}
		return s;
	}
	
	private void calculateStatistic(){
		stat = calculateStatistic(M);
	}
	
	private void caluculatePvalue(){
		List <Double> statNull = new ArrayList<Double>();
		while(statNull.size() < perm){
			double tmp = calculateStatistic(permutateMatrix(M));
			System.err.println(tmp);
			statNull.add(tmp);
		}
		pval = 0;
		for(Double d: statNull){
			if(d <= stat){
				pval++;
			}
		}
		if(pval==0){
			pval=1;
		}
		pval  /=  statNull.size();
	}
	
	public void perform(){
		calculateStatistic();
		caluculatePvalue();
	}
	
	private MyMat permutateMatrix(MyMat M){
		for(int i = 0; i < M.rowSize(); i++){
			List <Double> tmp = M.getRow(i);
			tmp = MyFunc.sample(tmp, tmp.size());
			for(int j = 0; j < M.colSize(); j++){
				M.set(i,j,tmp.get(j));
			}
		}
		return M;
	}
	
	private void print(){
		System.out.println(stat + "\t" + pval);
	}
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("p", "perm", true, "# of permutation");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] tabfile ", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 1){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] tabfile", options);
			return;
		}
			
		MutualExclusivityTest MET = new MutualExclusivityTest(new MyMat(argList.get(0)));
		if(commandLine.hasOption("p")){
			MET.perm = Integer.valueOf(commandLine.getOptionValue("p"));
		}
		MET.perform();
		MET.print();
	}
	
}
