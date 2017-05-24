package sim;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.*;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;


import sun.reflect.Reflection;

public class Simulator2D extends Simulator {
	// parameters for simulations
	int initialPopulationCoreSize = 10; // start with [initialPopulationCoreSize]^2 cells 
	
	//variable 
	long[][] genomes;
	int limXY;
	int maxX, minX, maxY, minY;
	boolean[][] rep;
	
	
	List <Long> getGenomesAsList(){
		List <Long> tmp = new ArrayList<Long>();
		for(int i = minX; i <= maxX; i++){
			for(int j = minY; j <= maxY; j++){
				if(!isEmpty(i,j)){
					tmp.add(get(i,j));
				}
			}
		}
		return tmp;
	}
	
	void  mutate(int x, int y){
		long genome = get(x, y);	 
		long sign  = (genome<0.1)?-1L:1L;
		Set <Integer> newMutatedGenes = null;	
		for(int i=0; i < genomeSize; i++){
			if(!isMutated(genome,i) & (R.nextDouble() < mutationRate)){
				if(newMutatedGenes==null){
					newMutatedGenes = new TreeSet <Integer>(genome2mutatedGenes(genome));
				}
				newMutatedGenes.add(i);
			}
		}
		if(newMutatedGenes==null){
			return;
		}
		set(x, y,  sign*mutatedGenes2genome(newMutatedGenes));
	}

	
	boolean isMutated(int x, int y, int g){
		return isMutated(get(x,y), g);
	}
	
	
	
	int coordinate2index(int c){
		return limXY + c; 
	}
	
	public Simulator2D(){};
	
	void initializeGenomes(){
		populationSize = 0;
		time = 0;
		limXY = (int) (10*Math.pow(maxPopulationSize, 0.5));
		genomes  = new long[2*limXY+1][2*limXY+1];
		rep  = new boolean [2*limXY+1][2*limXY+1];
		genome2fitnessMap =  new HashMap <Long, Double>();
		genome2mutatedGeneMap = new HashMap <Long, Set<Integer>>();
		fillCore();
	}
	
	void fillCore(){
		int coreSize = initialPopulationCoreSize/2;
		List <Integer> I = new ArrayList<Integer>();
		maxX = coreSize;
		minX = -coreSize;
		maxY = coreSize;
		minY = -coreSize;
		for(int i = -coreSize; i <= coreSize; i++){
			I.add(i);
		}
		for(Integer x: I){
			for(Integer y: I){
				//if(x +  y <= coreSize){
					set(x,y,1);
				//}
			}
		}
		populationSize = (int)Math.pow(2*coreSize+1, 2);
	}
	
	
	void set(int x, int y, long genome){
		genomes[coordinate2index(x)][coordinate2index(y)] = genome;
	}
	
	long get(int x, int y){
		return genomes[coordinate2index(x)][coordinate2index(y)];	
	}
	
	boolean getRep(int x, int y){
		return rep[coordinate2index(x)][coordinate2index(y)];	
	}
	
	
	void setRep(int x, int y, boolean b){
		rep[coordinate2index(x)][coordinate2index(y)] = b;	
	}
	
	public boolean isEmpty(int x, int y){
		if(get(x, y) == 0){
			return true;
		}else{
			return false;
		}
	}
	
	public void clear(int x, int y){
		set(x, y, 0);
		setRep(x, y, false);
	}
	
	
	public void createNeighbor(int x, int y, long genome){
		List <List<Integer>> neighbor  = getEmptyNeighbors(x, y);
		if(!neighbor.isEmpty()){
			List<Integer> tmp2 = neighbor.get(R.nextInt(neighbor.size()));
			int x2 = tmp2.get(0);
			int y2 = tmp2.get(1);
			set(x2,y2,genome);
			if(x2 > maxX){
				maxX = x2;
			}
			if(x2 < minX){
				minX = x2; 
			}
			if(y2 > maxY){
				maxY = y2;
			}
			if(y2 < minY){
				minY = y2; 
			}		
		}else{
			int[] v = new int[8];
			double[] V = new double[8];
			//right        
			for(int i=1; x+i<=limXY; i++){
				if(isEmpty(x+i,y)){
					v[0]=i-1;
					V[0]=1.0/(i-1);
					break;
				}
			}
			//upper right
			for(int i=1; x+i<=limXY & y+i<=limXY; i++){
				if(isEmpty(x+i,y+i)){
					v[1]=i-1;
					V[1]=1.0/(i-1);
					break;
				}
			}
			//upper
			for(int i=1; y+i<=limXY; i++){
				if(isEmpty(x,y+i) ){
					v[2]=i-1;
					V[2]=1.0/(i-1);
					break;
				}
			}
			//upper left
			for(int i=1; x-i>=-limXY & y+i<=limXY; i++){
				if(isEmpty(x-i,y+i) ){
					v[3]=i-1;
					V[3]=1.0/(i-1);
					break;
				}
			}
			//left
			for(int i=1; x-i>=-limXY; i++){
				if(isEmpty(x-i,y) ){
					v[4]=i-1;
					V[4]=1.0/(i-1);
					break;
				}
			}
			//lower left
			for(int i=1; x-i>=-limXY & y-i>=-limXY; i++){
				if(isEmpty(x-i,y-i) ){
					v[5]=i-1;
					V[5]=1.0/(i-1);
					break;
				}
			}
			//lower
			for(int i=1; y-i>=-limXY; i++){
				if(isEmpty(x,y-i)){
					v[6]=i-1;
					V[6]=1.0/(i-1);
					break;
				}
			}
			//lower right
			for(int i=1; y-i>=-limXY & x+i<=limXY; i++){
				if(isEmpty(x+i,y-i) ){
					v[7]=i-1;
					V[7]=1.0/(i-1);
					break;
				}
			}
			
			double s = 0;
			for(int i = 0; i < 8;i++){
				s += V[i];
			}
			for(int i = 0; i < 8;i++){
				V[i] /= s;
			}
			double r = R.nextDouble();
			s=0;
			int d = 0;
			int D = 0;
			for(int i = 0; i < 8;i++){
				s += V[i];
				if(r < s){
					d = v[i];
					D = i+1;
					break;
				}
			}
			
			long tmp[] = new long[d];
			boolean[] tmpR = new boolean[d];
			if(D==1){
				//right
				for(int i=1;i<=d;i++){
					tmp[i-1] = get(x+i, y);
					tmpR[i-1] = getRep(x+i, y);
				}
				set(x+1, y, genome);
				setRep(x+1, y, false);
				for(int i=1;i<=d;i++){
					set(x+i+1, y, tmp[i-1]);
					setRep(x+i+1, y, tmpR[i-1]);
				}
				if(x+d+1 > maxX){
					maxX = x+d+1;
				}
			}else if(D==2){
				//upper right
				for(int i=1;i<=d;i++){
					tmp[i-1] = get(x+i, y+i);
					tmpR[i-1] = getRep(x+i, y+i);
				}
				set(x+1, y+1, genome);
				setRep(x+1, y+1, false);
				for(int i=1;i<=d;i++){
					set(x+i+1, y+i+1,tmp[i-1]);
					setRep(x+i+1, y+i+1, tmpR[i-1]);
				}
				if(x+d+1 > maxX){
					maxX = x+d+1;
				}
				if( y+d+1 > maxY){
					maxY = y+d+1;
				}
			}else if(D==3){
				//upper
				for(int i=1;i<=d;i++){
					tmp[i-1] = get(x, y+i);
					tmpR[i-1] = getRep(x, y+i);
				}
				set(x, y+1, genome);
				setRep(x, y+1, false);
				for(int i=1;i<=d;i++){
					set(x, y+i+1,tmp[i-1]);
					setRep(x, y+i+1, tmpR[i-1]);
				}
				if( y+d+1 > maxY){
					maxY = y+d+1;
				}
			}else if(D==4){
				//upper left
				for(int i=1;i<=d;i++){
					tmp[i-1] = get(x-i, y+i);
					tmpR[i-1] = getRep(x-i, y+i);
				}
				set(x-1, y+1, genome);
				setRep(x-1, y+1, false);
				for(int i=1;i<=d;i++){
					set(x-i-1, y+i+1,tmp[i-1]);
					setRep(x-i-1, y+i+1, tmpR[i-1]);
				}
				if(x-d-1 < minX){
					minX = x-d-1;
				}
				if(y+d+1 > maxY){
					maxY = y+d+1;
				}
				
			}else  if(D==5){
				//left
				for(int i=1;i<=d;i++){
					tmp[i-1] = get(x-i, y);
					tmpR[i-1] = getRep(x-i, y);
				}
				set(x-1, y, genome);
				setRep(x-1, y, false);
				for(int i=1;i<=d;i++){
					set(x-i-1, y,tmp[i-1]);
					setRep(x-i-1, y, tmpR[i-1]);
				}
				if(x-d-1 < minX){
					minX = x-d-1;
				}
			}else if(D==6){
				//lower left
				for(int i=1;i<=d;i++){
					tmp[i-1] = get(x-i, y-i);
					tmpR[i-1] = getRep(x-i, y-i);
				}
				set(x-1, y-1, genome);
				setRep(x-1, y-1, false);
				for(int i=1;i<=d;i++){
					set(x-i-1, y-i-1,tmp[i-1]);
					setRep(x-i-1, y-i-1, tmpR[i-1]);
				}
				if(x-d-1 < minX){
					minX = x-d-1;
				}
				if(y-d-1 < minY){
					minY = y-d-1;
				}
			}else if(D==7){
				//lower 
				for(int i=1;i<=d;i++){
					tmp[i-1] = get(x, y-i);
					tmpR[i-1] = getRep(x, y-i);
				}
				set(x, y-1, genome);
				setRep(x, y-1, false);
				for(int i=1;i<=d;i++){
					set(x, y-i-1,tmp[i-1]);
					setRep(x, y-i-1, tmpR[i-1]);
				}
				if(y-d-1 < minY){
					minY = y-d-1;
				}
			}else if(D==8){
				//lower right
				for(int i=1;i<=d;i++){
					tmp[i-1] = get(x+i, y-i);
					tmpR[i-1] = getRep(x+i, y-i);
				}
				set(x+1, y-1, genome);
				setRep(x+1, y-1, false);
				for(int i=1;i<=d;i++){
					set(x+i+1, y-i-1,tmp[i-1]);
					setRep(x+i+1, y-i-1, tmpR[i-1]);
				}
				if(x+d+1 > maxX){
					maxX = x+d+1;
				}
				if(y-d-1 < minY){
					minY = y-d-1;
				}
			}
		}
		
	}
	
	
	void replicate(int x, int y){
		if(symmetricReplicationProbablity>=1){
			createNeighbor(x, y, get(x,y));
		}else{
			if(isStem(x,y)){
				if(R.nextDouble() < symmetricReplicationProbablity){
					createNeighbor(x, y, get(x,y));
				}else{
					createNeighbor(x, y, -get(x,y));
				}
			}else{
				createNeighbor(x, y, get(x,y));
			}
		}	
	}
	
	
	
	double getFitness(int x, int y){
		return genome2fitness(get(x,y));
	}	
	
	boolean isStem(int x, int y){
		return (get(x, y) > 0);
	}
	
	
	void growPopulation(){
		for(int x  = minX; x <= maxX ; x++){
			for(int y  = minY; y <= maxY ; y++){	
				if(!growCell1(x,y)){
					return;
				}
			}
		}
		if(!growCell2(0,0)){
			return;
		}
		for(int r  = 1; r <= Math.abs(maxX) | r <= Math.abs(minX) | r <= Math.abs(maxY) | r <= Math.abs(minY); r++){
			//if(true){
			if(R.nextBoolean()){
				for(int y = r; y > -r; y--){
					if(!growCell2(r,y)){
						return;
					}
				}
				for(int x = r; x > -r; x--){
					if(!growCell2(x,-r)){
						return;
					}
				}
				for(int y = -r; y < r; y++){
					if(!growCell2(-r,y)){
						return;
					}
				}
				for(int x = -r; x < r; x++){
					if(!growCell2(x,r)){
						return;
					}
				}
			}else{
				for(int x = r; x > -r; x--){
					if(!growCell2(x,r)){
						return;
					}
				}
				for(int y = r; y > -r; y--){
					if(!growCell2(-r,y)){
						return;
					}
				}
				for(int x = -r; x < r; x++){
					if(!growCell2(x,-r)){
						return;
					}
				}
				for(int y = -r; y < r; y++){
					if(!growCell2(r,y)){
						return;
					}
				}		
			}
		}
	}
	
	
	
	boolean growCell1(int x, int y){
		if(isEmpty(x,y)){
			return true;
		}
		if(Math.abs(x) ==  limXY | Math.abs(y) ==  limXY){
			return false;
		}
		if((isStem(x,y) & R.nextDouble() < deathRate) | (!isStem(x,y) & R.nextDouble() < deathRateForNonStem)){
			clear(x,y);
			populationSize--;
			return true;
		}
		if(R.nextDouble() < getFitness(x,y)*growthRate){
			mutate(x,y);
			setRep(x,y,true);
		}else{
			setRep(x,y,false);
		}
		if(populationSize >= maxPopulationSize){
			return false;
		}else{
			return true;
		}
	}
	
	boolean growCell2(int x, int y){
		if(isEmpty(x,y)){
			return true;
		}
		if(Math.abs(x) ==  limXY | Math.abs(y) ==  limXY){
			return false;
		}
		if(getRep(x,y)){
			setRep(x,y,false);
			replicate(x,y);
			populationSize++;
		}
		if(populationSize >= maxPopulationSize){
			return false;
		}else{
			return true;
		}
	}
	
	
	 List <List<Integer>> getEmptyNeighbors(int x, int y){
		 List <List<Integer>> tmp = new ArrayList <List<Integer>>();
		 for(int i = x-1; i<=x+1 ; i++){
			 for(int j = y-1; j<=y+1 ; j++){
				 if(isEmpty(i,j)){
					 List<Integer> tmp2 = new ArrayList<Integer>();
					 tmp2.add(i);
					 tmp2.add(j);
					 tmp.add(tmp2);
				 }
			 }
		 }
		 return tmp;
	 }
	 
	public void simulate(){
		 for(time = 1; time <= maxTime; time++){
			 growPopulation();
			 if(time/10 == (double)time/10){
				 System.err.println("time=" + time +"\t" + "populationSize=" + populationSize) ;
			}
			 if(Math.abs(maxX) ==  limXY | Math.abs(minX) ==  limXY | Math.abs(maxY) ==  limXY | Math.abs(minY) ==  limXY){
				 System.err.println("warn: reached space limit!");
				 break;
			}
			 if(populationSize>= maxPopulationSize){
				 break;
			}
		}
		System.err.println(minX + " " + maxX + " " + minY + " " + maxY + " " + limXY + " " + count());
	}
		
		
	public String toString(){
		StringBuffer S = new StringBuffer("");
		for(int y = maxY; y>=minY; y--){
			for(int x = minX; x<maxX; x++){
			 S.append(get(x,y) + "\t");	
			}
			S.append(get(maxX,y) + "\n");	
		}
		return S.toString();
	}
	
	
	public String getParmetersAsString(){
		StringBuffer S = new StringBuffer("");
		S.append("mutationRate" + "\t" +  mutationRate + "\n");
		S.append("growthRate" + "\t" + growthRate + "\n");
		S.append("deathRate" + "\t" + deathRate + "\n");
		S.append("deathRateForNonStem" + "\t" + deathRateForNonStem + "\n");
		S.append("symmetricReplicationProbablity" + "\t" + symmetricReplicationProbablity + "\n");
		S.append("initialPopulationCoreSize" + "\t" + initialPopulationCoreSize + "\n");
		S.append("maxTime" + "\t" + maxTime + "\n");
		S.append("maxPopulationSize" + "\t" + maxPopulationSize + "\n");
		S.append("genomeSize" + "\t" + genomeSize + "\n");
		S.append("driverSize" + "\t" + driverSize + "\n");
		S.append("fitnessIncrease" + "\t" + fitnessIncrease + "\n");
		S.append("limXY" + "\t" + limXY + "\n");
		S.append("maxX" + "\t" + maxX + "\n");
		S.append("minX" + "\t" + minX + "\n");
		S.append("maxY" + "\t" + maxY + "\n");
		S.append("minY" + "\t" + minY + "\n");
		S.append("populationSize" + "\t" + populationSize + "\n");
		S.append("time" + "\t" + time + "\n");
		 return S.toString();
	}
	
	/*public String toString(){
		StringBuffer S = new StringBuffer("");
		int max = maxX;
		if(max < -minX){
			max = -minX;
		}
		if(max < maxY){
			max = maxY;
		}
		if(max < -minY){
			max = -minY;
		}
		
		for(int y = max; y>=-max; y--){
			for(int x = -max; x<max; x++){
			 S.append(get(x,y) + "\t");	
			}
			S.append(get(max,y) + "\n");	
		}
		return S.toString();
	}*/
	
	double count(){
		int tmp = 0;
		for(int x  = minX; x <= maxX ; x++){
			for(int y  = minY; y <= maxY ; y++){	
				if(!isEmpty(x,y)){
					tmp++;
				}
			}
		}
		return tmp;
	}
	
	
	public static void main(String [] args) throws Exception{	
		Options options = new Options();
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		Simulator2D S = new Simulator2D();
		options.addOption("m", "mut", true,  "mutaionRate (" + S.mutationRate  + ")" );
		options.addOption("D", "death", true, "deathRate (" + S.deathRate  + ")");
		options.addOption("N", "deathns", true, "deathRate for non stem cell (" + S.deathRateForNonStem  + ")");
		options.addOption("P", "maxpop", true,  "maxPopulationSize (" + S.maxPopulationSize  + ")");
		options.addOption("G", "gen", true,  "genomeSize (" + S.genomeSize  + ")");
		options.addOption("g", "grow", true, "growthRate (" + S.growthRate  + ")");
		options.addOption("d", "drv", true,  "driverSize (" + S.driverSize  + ")");
		options.addOption("f", "fit", true,  "fitnessIncrease (" + S.fitnessIncrease  + ")");
		options.addOption("p", "inipop", true,  "initialPopulationCoreSize (" + S.initialPopulationCoreSize   + ")");
		options.addOption("T", "maxtime", true,  "maxTime (" + S.maxTime  + ")");
		options.addOption("S", "sym", true,  "symmetric replication probablity (" +  S.symmetricReplicationProbablity + ")");
		
		options.addOption("o", "pprm", true,  "print parameters");
		
		options.addOption("s", "stat", true,  "print statistics");
		
		
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] ", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 0){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] ", options);
			return;
		}
		

		if(commandLine.hasOption("m")){
			S.mutationRate = Double.valueOf(commandLine.getOptionValue("m"));
		}
		if(commandLine.hasOption("g")){
			S.growthRate = Double.valueOf(commandLine.getOptionValue("g"));
		}
		if(commandLine.hasOption("D")){
			S.deathRate = Double.valueOf(commandLine.getOptionValue("D"));
		}
		if(commandLine.hasOption("N")){
			S.deathRateForNonStem = Double.valueOf(commandLine.getOptionValue("N"));
		}
		if(commandLine.hasOption("p")){
			S.initialPopulationCoreSize  = Integer.valueOf(commandLine.getOptionValue("p"));
		}
		if(commandLine.hasOption("P")){
			S.maxPopulationSize = Integer.valueOf(commandLine.getOptionValue("P"));
		}
		if(commandLine.hasOption("G")){
			S.genomeSize = Integer.valueOf(commandLine.getOptionValue("G"));
		}
		if(commandLine.hasOption("d")){
			S.driverSize = Integer.valueOf(commandLine.getOptionValue("d"));
		}
		if(commandLine.hasOption("f")){
			S.fitnessIncrease = Double.valueOf(commandLine.getOptionValue("f"));
		}
		if(commandLine.hasOption("T")){
			S.maxTime = Integer.valueOf(commandLine.getOptionValue("T"));
		}
		if(commandLine.hasOption("S")){
			S.symmetricReplicationProbablity  = Double.valueOf(commandLine.getOptionValue("S"));
		}
		S.initializeGenomes();
		S.simulate();
		System.out.println(S);
		
		if(commandLine.hasOption("o")){
			PrintWriter os = new PrintWriter(new FileWriter(commandLine.getOptionValue("o")));
			os.print(S.getParmetersAsString());
			os.flush();
			os.close();
		}
		if(commandLine.hasOption("s")){
			PrintWriter os = new PrintWriter(new FileWriter(commandLine.getOptionValue("s")));
			StatisticsCalculator2D SC = new StatisticsCalculator2D(S);
			SC.calculateStatistics();
			os.print(SC);
			os.flush();
			os.close();
		}
		
	}

}
