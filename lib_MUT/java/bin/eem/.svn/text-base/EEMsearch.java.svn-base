package eem;

import java.io.IOException;
import java.util.List;
import java.util.Map;

public interface EEMsearch {
	void setPvalue1Cutoff(double d);		
	void suppressPvalue1Cutoff();
	void perform();
	ExpressionModuleSet getExpressionModuleSet();
	void printResults(String outfile) throws IOException;
	void printLog(String outfile) throws IOException;	
	Map <String, Double> getPvalues();
	void setGeneSets(Map<String, List<String>> geneSets);
	void setMaxGeneSetSize(int i);
	void setMinGeneSetSize(int i);
	void setItrForPvalue2Calculation(int i);
	void recycleNullDistribution();
	String getLog();
}
