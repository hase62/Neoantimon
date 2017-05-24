package bct;

import java.util.*;
import utility.*;
import org.apache.commons.cli.*;
import sun.reflect.Reflection;

public class BCTsampler {
	BCT  T;
	int c;
	int i;
	int j;
	int sumDiff;

	public static class MCMCtype{
		private final String name;
		private MCMCtype(String name){this.name = name;}
		public MCMCtype(MCMCtype b){this.name = b.name;}
		@Override
		public String toString(){return name;}
		public boolean equals(MCMCtype b){return name.equals(b.toString());};
		public static final MCMCtype  rowSumFixed = new MCMCtype("rowSumFixed");
		public static final MCMCtype columnSumFixed = new MCMCtype("columnSumFixed");
		public static final MCMCtype bothSumFixed = new MCMCtype("bothSumFixed");
	}
	
	
	interface MCMC{
		void next();
	}
	
	MCMCtype type;
	MCMC mcmc;
	
	public BCTsampler(MyMat M){
		T = new BCT(M);
		c = 0;
		sumDiff = 0;
		fixBothSum2();
	}
	
	public BCTsampler(BCT M){
		T = new BCT(M);
		c = 0;
		sumDiff = 0;
		fixBothSum2();
	}
	
	private int randCol(){
		return (int)(Math.random()*T.colSize());
	}
	private int randRow(){
		return (int)(Math.random()*T.colSize());
	}
	
	public void fixColSum(){
		type = MCMCtype.columnSumFixed;
		mcmc = new MCMC(){
			public void next(){
				for(int j = 0; j < T.colSize(); j++){
					int i = randRow();
					if(T.is1(i,j)){
						T.set0(i, j);
						while(true){
							int i2 = randRow();
							if(i != i2){
								T.set1(i2, j);
								break;
							}
						}
					}
				}
				c++;
			}
		};
	}
	
	public void fixRowSum(){
		type = MCMCtype.rowSumFixed;
		mcmc = new MCMC(){
			public void next(){
				for(int i = 0; i < T.rowSize(); i++){
					int j = randCol();
					if(T.is1(i,j)){
						T.set0(i, j);
						while(true){
							int j2 = randCol();
							if(j != j2){
								T.set1(i, j2);
								break;
							}
						}
					}
				}
				c++;
			}
		};
	}
	
	public void fixBothSum(){
		type = MCMCtype.bothSumFixed;
		mcmc = new MCMC(){
			public void next(){
				if(sumDiff == 0){
					i = randRow();
					j = randCol();
					if(T.is1(i, j)){
						T.set0(i, j);
						sumDiff = -1;
					}
				}else if(sumDiff == -1){
					int k = randRow();
					int l = randCol();
					if(i==k & j==l & T.is0(i, j)){
						T.set1(i, j);
						sumDiff = 0;
					}else if(T.is1(k, l)){
						if(Math.random() > 0.5){
							if(T.is0(k,j)){
								T.set1(k,j);
								T.set0(k,l);
								j=l;
							}
						}else{
							if(T.is0(i,l)){
								T.set1(i,l);
								T.set0(k,l);
								i=k;
							}
						}
					}
				}
				c++;
			}
		};
	}
	
	
	public void fixBothSum2(){
		type = MCMCtype.bothSumFixed;
		mcmc = new MCMC(){
			public void next(){
				if(sumDiff == 0){
					if(Math.random() > 0.5){
						i = randRow();
						j = randCol();
						if(T.is1(i, j)){
							T.set0(i, j);
							sumDiff = -1;
						}
					}else{
						i = randRow();
						j = randCol();
						if(T.is0(i, j)){
							T.set1(i, j);
							sumDiff = 1;
						}
					}
				}else if(sumDiff == -1){
					int k = randRow();
					int l = randCol();
					if(i==k & j==l & T.is0(i, j)){
						T.set1(i, j);
						sumDiff = 0;
					}else if(T.get(k, l)){
						if(Math.random() > 0.5){
							if(T.is0(k,j)){
								T.set1(k,j);
								T.set0(k,l);
								j=l;
							}
						}else{
							if(T.is0(i,l)){
								T.set1(i,l);
								T.set0(k,l);
								i=k;
							}
						}
					}
				}else if(sumDiff == 1){
					int k = randRow();
					int l = randCol();
					if(i==k & j==l & T.is1(i, j)){
							T.set0(i, j);
							sumDiff = 0;
					}else if(!T.get(k, l)){
						if(Math.random() > 0.5){
							if(T.is1(k,j)){
								T.set0(k,j);
								T.set1(k,l);
								j=l;
							}
						}else{
							if(T.is1(i,l)){
								T.set0(i,l);
								T.set1(k,l);
								i=k;
							}
						}
					}
				}
				c++;
			}
		};
	}
			
	
	public BCT getCurrentBCT(){
		return T;
	}
	
	public String getStatus(){
		return  sumDiff + " " + i + " " + j;
	}
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine = null;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] tabFile", options);
			System.exit(1);
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 1){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] tabFile", options);
			System.exit(1);
		}
		
		BCTsampler B = new BCTsampler(new MyMat(argList.get(0)));
		
		System.err.println(B.getCurrentBCT());
		for(int i = 0; i < 1000 ; i++){
			B.mcmc.next();
			System.err.println(B.getStatus() + "\n");
			System.err.println(B.getCurrentBCT());
		}
	}

}
