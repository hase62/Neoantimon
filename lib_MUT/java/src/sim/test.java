package sim;

public class test {
	public static void main(String [] args) throws Exception{
		
		/*for(int i = 0; i<1000; i++){
			char[] c = Character.toChars(i);
			for(int j = 0; j<c.length; j++){
			System.out.println(i + " " + c[j]);
			}
			}
		Long L = Long.MAX_VALUE;
		
		String S = Long.toBinaryString(L);
		System.out.println(S.length());
		*/
		
		
		
		BigSimulator2D S = new BigSimulator2D();
		S.mutationRate = 0.01;
		S.maxPopulationSize = 10000;
		S.symmetricReplicationProbablity = 1;
		S.initialPopulationSize = 10;
		S.maxTime = 1000000000;
		S.genomeSize = 10;
		//S.maxTime = 10;
		S.initializeGenomes();
		System.out.println(S.genome2mutatedGenes(S.get(0,0))+ " " + S.getFitness(0, 0));
		for(int i = 0; i < 100; i++){
			S.mutate(0, 0);
			System.out.println(S.genome2mutatedGenes(S.get(0,0)) + " " + S.getFitness(0, 0));
		}
		//S.simulate();
		
	}

}
