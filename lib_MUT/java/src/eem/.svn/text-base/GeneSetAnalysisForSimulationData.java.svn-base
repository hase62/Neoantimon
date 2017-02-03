package eem;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import eem.*;

import utility.*;

public class GeneSetAnalysisForSimulationData{
	
	interface geneSetAnalysis{
		Map <String, Double> getResults();
	}
	
	
	geneSetAnalysis GSA;
	SimulationDataGenerator SDG;
	
	MyMat Exp;
	Map <String, List<String>> geneSets;
	Map <String, Double> scores;
	
	public  GeneSetAnalysisForSimulationData(){	}
	
	public void generateSimulationData(int i){
		double moduleGeneSubsetRate[] = {0.2, 0.1, 0.05};
		double signalStrength[] = {0.4, 0.2, 0.1};
	
		double moduleGeneSubsetRate2[] = {0.2, 0.1, 0.05};
		double signalStrength2[] = {2, 1};
		double moduleSampleSubsetRate[] = {0.3, 0.1};
		
		double moduleGeneSubsetRate3[] = {0.4, 0.2, 0.1};
		double signalStrength3[] = {0.4, 0.2, 0.1};
		
		
		SDG = new SimulationDataGenerator();
		SimulationDataGeneratorUsingBicluster S = new SimulationDataGeneratorUsingBicluster();
		SimulationDataGeneratorUsingCompositeModule S2 = new SimulationDataGeneratorUsingCompositeModule();
		switch(i){
			//coherence
			case 1:	
				SDG.setSignalStrength(signalStrength[0]);
				SDG.setModuleGeneSubsetRate(moduleGeneSubsetRate[0]);
				break;
			case 2:
				SDG.setSignalStrength(signalStrength[0]);
				SDG.setModuleGeneSubsetRate(moduleGeneSubsetRate[1]);
				break;
			case 3:
				SDG.setSignalStrength(signalStrength[0]);
				SDG.setModuleGeneSubsetRate(moduleGeneSubsetRate[2]);
				break;
			case 4:
				SDG.setSignalStrength(signalStrength[1]);
				SDG.setModuleGeneSubsetRate(moduleGeneSubsetRate[0]);
				break;
			case 5:
				SDG.setSignalStrength(signalStrength[1]);
				SDG.setModuleGeneSubsetRate(moduleGeneSubsetRate[1]);
				break;
			case 6:
				SDG.setSignalStrength(signalStrength[1]);
				SDG.setModuleGeneSubsetRate(moduleGeneSubsetRate[2]);
				break;
			case 7:
				SDG.setSignalStrength(signalStrength[2]);
				SDG.setModuleGeneSubsetRate(moduleGeneSubsetRate[0]);
				break;
			case 8:
				SDG.setSignalStrength(signalStrength[2]);
				SDG.setModuleGeneSubsetRate(moduleGeneSubsetRate[1]);
				break;
			case 9:
				SDG.setSignalStrength(signalStrength[2]);
				SDG.setModuleGeneSubsetRate(moduleGeneSubsetRate[2]);
				break;	
			
			//bicluster
			case 10:
				S.setModuleGeneSubsetRate(moduleGeneSubsetRate2[0]);
				S.setSignalStrength(signalStrength2[0]);
				S.setModuleSampleSubsetRate(moduleSampleSubsetRate[0]);
				SDG = S;
				break;	
			case 11:
				S.setModuleGeneSubsetRate(moduleGeneSubsetRate2[0]);
				S.setSignalStrength(signalStrength2[0]);
				S.setModuleSampleSubsetRate(moduleSampleSubsetRate[1]);
				SDG = S;
				break;	
			case 12:
				S.setModuleGeneSubsetRate(moduleGeneSubsetRate2[0]);
				S.setSignalStrength(signalStrength2[1]);
				S.setModuleSampleSubsetRate(moduleSampleSubsetRate[0]);
				SDG = S;
				break;	
			case 13:
				S.setModuleGeneSubsetRate(moduleGeneSubsetRate2[0]);
				S.setSignalStrength(signalStrength2[1]);
				S.setModuleSampleSubsetRate(moduleSampleSubsetRate[1]);
				SDG = S;
				break;		
			case 14:
				S.setModuleGeneSubsetRate(moduleGeneSubsetRate2[1]);
				S.setSignalStrength(signalStrength2[0]);
				S.setModuleSampleSubsetRate(moduleSampleSubsetRate[0]);
				SDG = S;
				break;	
			case 15:
				S.setModuleGeneSubsetRate(moduleGeneSubsetRate2[1]);
				S.setSignalStrength(signalStrength2[0]);
				S.setModuleSampleSubsetRate(moduleSampleSubsetRate[1]);
				SDG = S;
				break;	
			case 16:
				S.setModuleGeneSubsetRate(moduleGeneSubsetRate2[1]);
				S.setSignalStrength(signalStrength2[1]);
				S.setModuleSampleSubsetRate(moduleSampleSubsetRate[0]);
				SDG = S;
				break;	
			case 17:
				S.setModuleGeneSubsetRate(moduleGeneSubsetRate2[1]);
				S.setSignalStrength(signalStrength2[1]);
				S.setModuleSampleSubsetRate(moduleSampleSubsetRate[1]);
				SDG = S;
				break;
			case 18:
				S.setModuleGeneSubsetRate(moduleGeneSubsetRate2[2]);
				S.setSignalStrength(signalStrength2[0]);
				S.setModuleSampleSubsetRate(moduleSampleSubsetRate[0]);
				SDG = S;
				break;	
			case 19:
				S.setModuleGeneSubsetRate(moduleGeneSubsetRate2[2]);
				S.setSignalStrength(signalStrength2[0]);
				S.setModuleSampleSubsetRate(moduleSampleSubsetRate[1]);
				SDG = S;
				break;		
			case 20:
				S.setModuleGeneSubsetRate(moduleGeneSubsetRate2[2]);
				S.setSignalStrength(signalStrength2[1]);
				S.setModuleSampleSubsetRate(moduleSampleSubsetRate[0]);
				SDG = S;
				break;		
			case 21:
				S.setModuleGeneSubsetRate(moduleGeneSubsetRate2[2]);
				S.setSignalStrength(signalStrength2[1]);
				S.setModuleSampleSubsetRate(moduleSampleSubsetRate[1]);
				SDG = S;
				break;		
				
				//composite module
			case 22:	
				S2.setSignalStrength(signalStrength3[0]);
				S2.setModuleGeneSubsetRate(moduleGeneSubsetRate3[0]);
				SDG = S2;
				break;
			case 23:
				S2.setSignalStrength(signalStrength3[0]);
				S2.setModuleGeneSubsetRate(moduleGeneSubsetRate3[1]);
				SDG = S2;
				break;
			case 24:
				S2.setSignalStrength(signalStrength3[0]);
				S2.setModuleGeneSubsetRate(moduleGeneSubsetRate3[2]);
				SDG = S2;
				break;
			case 25:
				S2.setSignalStrength(signalStrength3[1]);
				S2.setModuleGeneSubsetRate(moduleGeneSubsetRate3[0]);
				SDG = S2;
				break;
			case 26:
				S2.setSignalStrength(signalStrength3[1]);
				S2.setModuleGeneSubsetRate(moduleGeneSubsetRate3[1]);
				SDG = S2;
				break;
			case 27:
				S2.setSignalStrength(signalStrength3[1]);
				S2.setModuleGeneSubsetRate(moduleGeneSubsetRate3[2]);
				SDG = S2;
				break;
			case 28:
				S2.setSignalStrength(signalStrength3[2]);
				S2.setModuleGeneSubsetRate(moduleGeneSubsetRate3[0]);
				SDG = S2;
				break;
			case 29:
				S2.setSignalStrength(signalStrength3[2]);
				S2.setModuleGeneSubsetRate(moduleGeneSubsetRate3[1]);
				SDG = S2;
				break;
			case 30:
				S2.setSignalStrength(signalStrength3[2]);
				S2.setModuleGeneSubsetRate(moduleGeneSubsetRate3[2]);
				SDG = S2;
				break;	
				
				
		}
		SDG.simulate();
		Exp = SDG.getExpression();
		geneSets = SDG.getGeneSet();
	}
	public void performGeneSetAnalysis(int i){
		switch(i){
			case 1:
				GSA = new geneSetAnalysis(){
					public Map <String, Double> getResults(){
						OneSampleOverExpression OSOE = new OneSampleOverExpression(Exp);
						OSOE.setGeneSets(geneSets);
						OSOE.perform();
						return OSOE.getPvalues();
					}
				};
				break;
			case 2:
				GSA = new geneSetAnalysis(){
						public Map <String, Double> getResults(){
							GeneSetExpressionCoherenceAnalysis GSECA= new GeneSetExpressionCoherenceAnalysis(Exp);
							GSECA.setGeneSets(geneSets);
							GSECA.setItrForPvalueCalculation(1000);
							GSECA.supressItrAdjustment();
							GSECA.perform();
							return GSECA.getPvalues();
						}
				};
				break;
			
			case 3:
				GSA = new geneSetAnalysis(){
						public Map <String, Double> getResults(){
							ParameterOptimizingCoherenceBasedEEMsearch CE = new ParameterOptimizingCoherenceBasedEEMsearch(Exp);
							CE.setGeneSets(geneSets);
							CE.recycleNullDistribution();
							CE.suppressPvalue1Cutoff();
							CE.perform();
							return CE.getPvalues();
					}
				};
				break;
			case 4:
				GSA = new geneSetAnalysis(){
						public Map <String, Double> getResults(){
							ParameterOptimizingBiclusterBasedEEMsearch BE = new ParameterOptimizingBiclusterBasedEEMsearch(Exp);
							BE.setGeneSets(geneSets);
							BE.recycleNullDistribution();
							BE.suppressPvalue1Cutoff();
							BE.perform();
							return BE.getPvalues();
					}
				};
				break;
			case 5:
				GSA = new geneSetAnalysis(){
						public Map <String, Double> getResults(){
							CoherenceBasedEEMsearch CE = new CoherenceBasedEEMsearch(Exp);
							CE.setGeneSets(geneSets);
							CE.recycleNullDistribution();
							CE.calculateCor();
							CE.suppressPvalue1Cutoff();
							CE.setRelativeRadius(0.05);
							CE.perform();
							return CE.getPvalues();
					}
				};
				break;
			case 6:
				GSA = new geneSetAnalysis(){
						public Map <String, Double> getResults(){
							BiclusterBasedEEMsearch BE = new BiclusterBasedEEMsearch(Exp);
							BE.setGeneSets(geneSets);
							BE.recycleNullDistribution();
							BE.suppressPvalue1Cutoff();
							BE.perform();
							return BE.getPvalues();
					}
				};
				break;
			
		}
		scores = GSA.getResults();
	}
	
	public List<Double> getScoresOfPositiveGeneSet(){
		List <Double> tmp = new ArrayList<Double>();
		for(String s: scores.keySet()){
			if(s.startsWith("positive")){
				tmp.add(scores.get(s));
			}
		}
		return tmp;	
	}
	
	public List<Double> getScoresOfNegativeGeneSet(){
		List <Double> tmp = new ArrayList<Double>();
		for(String s: scores.keySet()){
			if(s.startsWith("negative")){
				tmp.add(scores.get(s));
			}
		}
		return tmp;	
	}
	
	
	public void printResults2OutFile(String outFile) throws IOException{
		PrintWriter os = new PrintWriter(new FileWriter(outFile));
		os.println(MyFunc.join("\t", getScoresOfPositiveGeneSet())); 
		os.println(MyFunc.join("\t", getScoresOfNegativeGeneSet())); 
		os.flush();
		os.close();
	}
	
	
	public static void main(String[] args) throws Exception {
		int simulationDataType = Integer.valueOf(args[0]);
		int geneSetAnalysisType = Integer.valueOf(args[1]);
		String outFile = args[2];
		GeneSetAnalysisForSimulationData G = new GeneSetAnalysisForSimulationData();
		G.generateSimulationData(simulationDataType);
		G.performGeneSetAnalysis(geneSetAnalysisType);
		G.printResults2OutFile(outFile);
	}
}
