package sim;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

import sun.reflect.Reflection;

public class BigSimulator2DwithEssensialGenes {
	
	int maxBinary = 63;// 	number of digits of max long  is 63 in  the binary system 
	// parameters for simulations
	double mutationRate = 0.001;
	double growthRate = 0.0001;
	double deathRate =  0.0000001;
	double deathRateForNonStem =  0.01;
	double symmetricReplicationProbablity = 0.1;
	int maxTime = 5000000;
	int maxPopulationSize = 100000;
	int genomeSize = 300;
	int  driverSize = 10; //  driverSize  +  essensialGeneSize  <= maxBinary
	int essensialGeneSize = 0; 
	
	double fitnessIncrease  = 5;
	int initialPopulationSize = 10; 
	
	//variable
	int G; // the number of long integer used for coding genomes 
	long[][][] genomes;
	boolean [][][] state; // occupied differentiated  replicate 
	int limXY;
	int maxX, minX, maxY, minY;
	Map <Long, Double> genome2fitnessMap;
	int populationSize;
	int time;
	Map <Long, Set<Integer>> genome2mutatedGeneMap;	
	
	//for simulation
	Random R = new Random();
	
	public BigSimulator2DwithEssensialGenes(){};
	
	void initializeGenomes() throws Exception{
		if(driverSize  +  essensialGeneSize  > maxBinary |  driverSize  +  essensialGeneSize  > genomeSize){
			throw new Exception();
		}
		populationSize = 0;
		time = 0;
		limXY = (int) (2*Math.pow(maxPopulationSize, 0.5));
		G = genomeSize/maxBinary + 1; 
		genomes  = new long[2*limXY+1][2*limXY+1][G];
		state= new boolean [2*limXY+1][2*limXY+1][3];
		genome2fitnessMap =  new HashMap <Long, Double>();
		genome2mutatedGeneMap = new HashMap <Long, Set<Integer>>();
		fillCore();
	}
	
	static String decimal2binary(long d){
		return Long.toBinaryString(d);
	}
	
	static long binary2decimalLong(String numberString){
        long convertedNumber = 0L; 
        for (int i = 0; i < numberString.length(); i++){
        	if (Integer.parseInt(Character.toString(numberString.charAt(i))) == 1){
        		convertedNumber = convertedNumber + (Long.parseLong(Character.toString(numberString.charAt(i))) * (long)Math.pow(2, (numberString.length() - i - 1)));
        	}
        }
        return convertedNumber;	
	}
	
	long mutatedGenes2genomeLong(Set <Integer> mutatedGenes){
		String tmp = "";	
		for(int i = 0; i <  maxBinary; i++){
			if(mutatedGenes.contains(i)){
				tmp = "1" + tmp;
			}else{
				tmp = "0" + tmp;
			}
		}
		long tmp2 =  binary2decimalLong(tmp);
		genome2mutatedGeneMap.put(tmp2, mutatedGenes);
		return tmp2;
	}
	
	Set <Integer>  genomeLong2mutatedGenes(long genome){
		if(genome2mutatedGeneMap.containsKey(genome)){
			return genome2mutatedGeneMap.get(genome);
		}else{
			String tmp = decimal2binary(genome);
			Set <Integer> mutatedGenes = new TreeSet<Integer>();
			int l = tmp.length();
			for(int i = 0; i < l; i++){
				if(tmp.substring(l-i-1,l-i).equals("1")){
					mutatedGenes.add(i);
				}
			}
			genome2mutatedGeneMap.put(genome, mutatedGenes);
			return mutatedGenes;	
		}
	}
	
	
	Set <Integer>  genome2mutatedGenes(long genome[]){
		Set <Integer> mutatedGenes = new TreeSet<Integer>();
		for(int j = 0; j < G; j++){
			Set <Integer> tmp =  genomeLong2mutatedGenes(genome[j]);
			for(int i = 0;(i < maxBinary & j*maxBinary + i < genomeSize); i++){
				if(tmp.contains(i)){
					mutatedGenes.add(j*maxBinary+i);
				}
			}
		}
		return mutatedGenes;				
	}
	
	int coordinate2index(int c){
		return limXY + c; 
	}
	
	long[] get(int x, int y){
		return genomes[coordinate2index(x)][coordinate2index(y)];	
	}
	
	long[] getDeeply(int x, int y){
		long[] tmp = new long[G];
		for(int i = 0; i < G ; i++){
			tmp[i] = genomes[coordinate2index(x)][coordinate2index(y)][i];
		}
		return tmp;
	}

	void set(int x, int y, long[] genome){
		long[] tmp = new long[G];
		for(int i = 0; i < G ; i++){
			tmp[i] = genome[i];
		}
		genomes[coordinate2index(x)][coordinate2index(y)] = tmp;
	}
	
	
	void set(int x, int y, long[] genome, boolean[] state){
		set(x, y, genome);
		setState(x, y, state);
	}
	
	void clear(int x, int y){
		genomes[coordinate2index(x)][coordinate2index(y)] = new long[G];
		state[coordinate2index(x)][coordinate2index(y)][0] = false; //occupied
		state[coordinate2index(x)][coordinate2index(y)][1] = false; //differentiated 
		state[coordinate2index(x)][coordinate2index(y)][2] = false; //replicate
	}
	
	

	boolean[] getState(int x, int y){
		return state[coordinate2index(x)][coordinate2index(y)];	
	}
	
	
	boolean[] getStateDeeply(int x, int y){
		boolean[] tmp = new boolean[3];
		for(int i = 0; i < 3 ; i++){
			tmp[i] = state[coordinate2index(x)][coordinate2index(y)][i];
		}
		return tmp;	
	}
	
	void setState(int x, int y, boolean[] s){
		boolean[] tmp = new boolean[3];
		for(int i = 0; i < 3 ; i++){
			tmp[i] = s[i];
		}
		state[coordinate2index(x)][coordinate2index(y)] = tmp;
	}
	
	boolean isEmpty(int x, int y){
		return !state[coordinate2index(x)][coordinate2index(y)][0];
	}
	
	boolean getRep(int x, int y){
		return state[coordinate2index(x)][coordinate2index(y)][2];	
	}
	
	void setRep(int x, int y, boolean b){
		state[coordinate2index(x)][coordinate2index(y)][2] = b;	
	}
	
	void setDif(int x, int y, boolean b){
		state[coordinate2index(x)][coordinate2index(y)][1] = b;
	}
	
	boolean isStem(int x, int y){
		return !state[coordinate2index(x)][coordinate2index(y)][1];
	}
	
	boolean isMutated(int x, int y, int g){
		return isMutated(get(x,y), g);
	}
	
	boolean isMutated(long[] genome, int g){
		int n = g/maxBinary;   //  g = n*maxBinary + m
		int m = g - (n*maxBinary);
		return genomeLong2mutatedGenes(genome[n]).contains(m);
	}
	
	double genome2fitness(long[] genome){
		if(genome2fitnessMap.containsKey(genome[0])){
			return genome2fitnessMap.get(genome[0]);
		}else{
			double fitness = 1;
			for(int i=0; i< driverSize; i++){
				if(isMutated(genome,i)){
					fitness *= fitnessIncrease;
				}
			}
			for(int i=driverSize; i< driverSize + essensialGeneSize; i++){
				if(isMutated(genome,i)){
					fitness = -1;
				}
			}
			genome2fitnessMap.put(genome[0], fitness);
			return fitness;
		}
	}
	
	
	void fillCore(){
		setNormalCell(0,0);
		populationSize++;
		if(populationSize==initialPopulationSize){
			return;
		}
		for(int r  = 1; r <limXY; r++){
				for(int y = r; y > -r; y--){
					setNormalCell(r,y);
					setMaxMinXY(r, y);
					populationSize++;
					if(populationSize==initialPopulationSize){
						return;
					}
				}
				for(int x = r; x > -r; x--){
					setNormalCell(x,-r);
					setMaxMinXY(x, -r);
					populationSize++;
					if(populationSize==initialPopulationSize){
						return;
					}
				}
				for(int y = -r; y < r; y++){
					setNormalCell(-r,y);
					setMaxMinXY(-r, y);
					populationSize++;
					if(populationSize==initialPopulationSize){
						return;
					}
				}
				for(int x = -r; x < r; x++){
					setNormalCell(x,r);
					setMaxMinXY(x,r);
					populationSize++;
					if(populationSize==initialPopulationSize){
						return;
					}
				}		
		}
	}
	
	void setNormalCell(int x, int y){
		long[] tmp = new long[G];
		boolean[] tmp2 = new boolean[3];
		tmp2[0] = true;
		set(x, y,tmp, tmp2);
	}
	
	void setMaxMinXY(int x, int y){
		if(x > maxX){
			maxX = x;
		}
		if(y > maxY){
			maxY = y;
		}
		if(x < minX){
			minX = x;
		}
		if(y < minY){
			minY = y;
		}
	}
	

	void  mutate(int x, int y){
		long genome[] = get(x, y);	 
		for(int j = 0; j < G; j++){
			Set <Integer> newMutatedGenes = null;	
			for(int i=0; (i < maxBinary & j*maxBinary + i < genomeSize); i++){
				if(!isMutated(genome,j*maxBinary + i) & (R.nextDouble() < mutationRate)){
					if(newMutatedGenes==null){
						newMutatedGenes = new TreeSet <Integer>(genomeLong2mutatedGenes(genome[j]));
					}
					newMutatedGenes.add(i);
				}
			}
		
			if(newMutatedGenes!=null){
				genome[j] = mutatedGenes2genomeLong(newMutatedGenes);
			}
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
	
	public void createNeighbor(int x, int y){
		setRep(x,y,false);
		long[] genome = getDeeply(x,y); 
		boolean[] state = getStateDeeply(x, y);
		List <List<Integer>> neighbor  = getEmptyNeighbors(x, y);
		if(!neighbor.isEmpty()){
			List<Integer> tmp2 = neighbor.get(R.nextInt(neighbor.size()));
			int x2 = tmp2.get(0);
			int y2 = tmp2.get(1);
			set(x2,y2,genome,state);
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
			
			long[][] tmp = new long[d][];
			boolean[][] tmpR = new boolean[d][];
			if(D==1){
				//right
				for(int i=1;i<=d;i++){
					tmp[i-1] = getDeeply(x+i, y);
					tmpR[i-1] = getStateDeeply(x+i, y);
				}
				set(x+1, y, genome, state);
				for(int i=1;i<=d;i++){
					set(x+i+1, y, tmp[i-1], tmpR[i-1]);
				}
				if(x+d+1 > maxX){
					maxX = x+d+1;
				}
			}else if(D==2){
				//upper right
				for(int i=1;i<=d;i++){
					tmp[i-1] = getDeeply(x+i, y+i);
					tmpR[i-1] = getStateDeeply(x+i, y+i);
				}
				set(x+1, y+1, genome, state);
				for(int i=1;i<=d;i++){
					set(x+i+1, y+i+1,tmp[i-1],tmpR[i-1]);
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
					tmp[i-1] = getDeeply(x, y+i);
					tmpR[i-1] = getStateDeeply(x, y+i);
				}
				set(x, y+1, genome, state);
				for(int i=1;i<=d;i++){
					set(x, y+i+1,tmp[i-1], tmpR[i-1]);
				}
				if( y+d+1 > maxY){
					maxY = y+d+1;
				}
			}else if(D==4){
				//upper left
				for(int i=1;i<=d;i++){
					tmp[i-1] = getDeeply(x-i, y+i);
					tmpR[i-1] = getStateDeeply(x-i, y+i);
				}
				set(x-1, y+1, genome, state);
				for(int i=1;i<=d;i++){
					set(x-i-1, y+i+1,tmp[i-1],tmpR[i-1]);
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
					tmp[i-1] = getDeeply(x-i, y);
					tmpR[i-1] = getStateDeeply(x-i, y);
				}
				set(x-1, y, genome, state);
				for(int i=1;i<=d;i++){
					set(x-i-1, y, tmp[i-1],  tmpR[i-1]);
				}
				if(x-d-1 < minX){
					minX = x-d-1;
				}
			}else if(D==6){
				//lower left
				for(int i=1;i<=d;i++){
					tmp[i-1] = getDeeply(x-i, y-i);
					tmpR[i-1] = getStateDeeply(x-i, y-i);
				}
				set(x-1, y-1, genome, state);
				for(int i=1;i<=d;i++){
					set(x-i-1, y-i-1,tmp[i-1], tmpR[i-1]);
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
					tmp[i-1] = getDeeply(x, y-i);
					tmpR[i-1] = getStateDeeply(x, y-i);
				}
				set(x, y-1, genome, state);
				for(int i=1;i<=d;i++){
					set(x, y-i-1,tmp[i-1], tmpR[i-1]);
				}
				if(y-d-1 < minY){
					minY = y-d-1;
				}
			}else if(D==8){
				//lower right
				for(int i=1;i<=d;i++){
					tmp[i-1] = getDeeply(x+i, y-i);
					tmpR[i-1] = getStateDeeply(x+i, y-i);
				}
				set(x+1, y-1, genome, state);
				for(int i=1;i<=d;i++){
					set(x+i+1, y-i-1,tmp[i-1],  tmpR[i-1]);
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
			if(getFitness(x,y) < 0){
				clear(x,y);
				populationSize--;
				return true;
			}
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
	
	double getFitness(int x, int y){
		return genome2fitness(get(x,y));
	}	
	
	
	void replicate(int x, int y){	
		if(symmetricReplicationProbablity<1 & isStem(x,y)){
			//stem
			if(R.nextDouble() > symmetricReplicationProbablity){
				//asymmetric replication
				if(R.nextBoolean()){
					setDif(x, y, true);
					createNeighbor(x, y);
					setDif(x, y, false);
				}else{
					createNeighbor(x, y);
					setDif(x, y, true);
				}
			}else{
				//symmetric replication
				createNeighbor(x, y);
			}
		}else{
			//differentiated
			createNeighbor(x, y);
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
			replicate(x,y);
			populationSize++;
		}
		if(populationSize >= maxPopulationSize){
			return false;
		}else{
			return true;
		}
	}
	
	public void simulateWhilePrinting(String f, int I) throws IOException{
		int i = 0;
		printResults(f + "." + i);
		int poplationSizeCutoff = populationSize*I;
		 for(time = 1; time <= maxTime; time++){
			 growPopulation();
			 if(time/10 == (double)time/10){
				 System.err.println("time=" + time +"\t" + "populationSize=" + populationSize) ;
			}
			 if(Math.abs(maxX) ==  limXY | Math.abs(minX) ==  limXY | Math.abs(maxY) ==  limXY | Math.abs(minY) ==  limXY){
				 System.err.println("warn: reached space limit!");
				 break;
			}
			if(populationSize >= poplationSizeCutoff){
				i++;
				printResults(f + "." + i);
				poplationSizeCutoff *= I;
			}
			 if(populationSize == 0 | populationSize>= maxPopulationSize){
				 break;
			}
		}
		System.err.println(minX + " " + maxX + " " + minY + " " + maxY + " " + limXY + " " + count());
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
			 if(populationSize == 0 | populationSize>= maxPopulationSize){
				 break;
			}
		}
		System.err.println(minX + " " + maxX + " " + minY + " " + maxY + " " + limXY + " " + count());
	}
		
	
	public void printResults(String f) throws IOException{
		PrintWriter os = new PrintWriter(new FileWriter(f + ".tm"));
		PrintWriter os2 = new PrintWriter(new FileWriter(f + ".mut"));
		PrintWriter os3 = new PrintWriter(new FileWriter(f + ".prm"));
		Map <String, Integer>  tmp = new HashMap <String, Integer>();	
		int i = 0;
		
		/*int M = maxX;
		if(-minX > M){
			M = -minX;
		}
		if(maxY > M){
			M = maxY;
		}
		if(-minY > M){
			M = -minY;
		}
		maxX = M;
		
		for(int y = M; y>= -M; y--){
			for(int x = -M; x<=M; x++){
		*/
		
		for(int y = maxY; y>=minY; y--){
			for(int x = minX; x<=maxX; x++){
				if(isEmpty(x,y)){
					 os.print(0);			
				}else{
					Set <Integer> tmp2 = genome2mutatedGenes(get(x,y));
					if(!tmp.containsKey(tmp2.toString())){
						i++;
						tmp.put(tmp2.toString(), i);
						os2.print(i);
						//os2.print("\t" + tmp2.toString());
						for(int g = 0; g < genomeSize; g++){
							os2.print("\t" + (tmp2.contains(g)?1:0));
						}
						os2.print("\n");
					}
					os.print((isStem(x,y)?"":"-") + tmp.get(tmp2.toString()));
				}
				
				if(x ==  maxX){
					os.print("\n");
				}else{
					os.print("\t");
				}
			}
		}
		os3.print(getParmetersAsString());
		os.flush();
		os2.flush();
		os3.flush();
		os.close();
		os2.close();
		os3.close();
	}
	
	public String getParmetersAsString(){
		StringBuffer S = new StringBuffer("");
		S.append("mutationRate" + "\t" +  mutationRate + "\n");
		S.append("growthRate" + "\t" + growthRate + "\n");
		S.append("deathRate" + "\t" + deathRate + "\n");
		S.append("deathRateForNonStem" + "\t" + deathRateForNonStem + "\n");
		S.append("symmetricReplicationProbablity" + "\t" + symmetricReplicationProbablity + "\n");
		S.append("initialPopulationSize" + "\t" + initialPopulationSize + "\n");
		S.append("maxTime" + "\t" + maxTime + "\n");
		S.append("maxPopulationSize" + "\t" + maxPopulationSize + "\n");
		S.append("genomeSize" + "\t" + genomeSize + "\n");
		S.append("driverSize" + "\t" + driverSize + "\n");
		S.append("essensialGeneSize" + "\t" + essensialGeneSize + "\n");
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
		BigSimulator2DwithEssensialGenes S = new BigSimulator2DwithEssensialGenes();
		String outfile = "out";
		options.addOption("m", "mut", true,  "mutaionRate (" + S.mutationRate  + ")" );
		options.addOption("D", "death", true, "deathRate (" + S.deathRate  + ")");
		options.addOption("N", "deathns", true, "deathRate for non stem cell (" + S.deathRateForNonStem  + ")");
		options.addOption("P", "maxpop", true,  "maxPopulationSize (" + S.maxPopulationSize  + ")");
		options.addOption("G", "gen", true,  "genomeSize (" + S.genomeSize  + ")");
		options.addOption("g", "grow", true, "growthRate (" + S.growthRate  + ")");
		options.addOption("d", "drv", true,  "driverSize (" + S.driverSize  + ")");
		options.addOption("e", "ess", true,  "essensialGeneSize (" + S.essensialGeneSize  + ")");
		options.addOption("f", "fit", true,  "fitnessIncrease (" + S.fitnessIncrease  + ")");
		options.addOption("p", "inipop", true,  "initialPopulationSize (" + S.initialPopulationSize   + ")");
		options.addOption("T", "maxtime", true,  "maxTime (" + S.maxTime  + ")");
		options.addOption("S", "sym", true,  "symmetric replication probablity (" +  S.symmetricReplicationProbablity + ")");
		
		options.addOption("o", "outfile", true,  "out filename");
		options.addOption("s", "snap", true,  "take snap shot");
		
	
		
		
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
			S.initialPopulationSize  = Integer.valueOf(commandLine.getOptionValue("p"));
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
		if(commandLine.hasOption("e")){
			S.essensialGeneSize = Integer.valueOf(commandLine.getOptionValue("e"));
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
		if(commandLine.hasOption("o")){
			outfile = commandLine.getOptionValue("o");
		}		
		S.initializeGenomes();
		if(!commandLine.hasOption("s")){
			S.simulate();
		}else{
			S.simulateWhilePrinting(outfile, Integer.valueOf(commandLine.getOptionValue("s")));
		}
		S.printResults(outfile);	
	}
	
	

}
