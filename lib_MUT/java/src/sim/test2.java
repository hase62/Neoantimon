package sim;

public class test2 {
	public static void main(String [] args) throws Exception{
		BigSimulator2DwithSelection S = new  BigSimulator2DwithSelection();
		S.maxPopulationSize = 10;
		S.genomeSize=5;
		S.driverSize=5;
		S.mutationRate = 0.1;
		S.initializeGenomes();
		S.setSelectiveRegion1();
		int x = -1;
		int y = 1;
		System.out.println(S.selection2selectiveDriverGenes(S.getSelection(x,y)));
		System.out.println(S.genome2mutatedGenes(S.get(x,y))+ " " + S.getFitness(x, y));
		for(int i = 0; i < 20; i++){
			S.mutate(x, y);
			System.out.println(S.genome2mutatedGenes(S.get(x,y)) + " " + S.getFitness(x, y));
		}
	}
}
