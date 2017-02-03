package eem;


import utility.*;
import utility.HierarchicalClustering.DistFuncType;

import java.util.*;

public class test {
	
	
	public static void main(String[] args) throws Exception {
		MyMat M = new MyMat(args[0]);
		Map<String, List<String>> G = MyFunc.readGeneSetFromGmtFile(args[1]);
		M.normalizeRows();
		M = M.getSubMatByRow(G.get(args[2]));
		//M = M.getCovMatForRow();
	
		ClusteredMyMat M2  = new ClusteredMyMat(M);
		
		
		
		
		//ClusteredMyMat M2  = new ClusteredMyMat(M.getSubMatByRow(MyFunc.readStringList(args[1])));
		//M2.normalizeRows();
		//M2.supressRowClustering();
		M2.performClustering();
		//M2.setRowDistFunc(DistFuncType.ABSOLUTE_CORRELATION);
		
		ClusteredMyMatViewer  MV = new ClusteredMyMatViewer(M2);
		MV.scaleColorByRow();
		MV.setOutFile(args[2] + ".pdf");
		MV.useAutoStop();
		MV.init();
		
		Heatmap H = new Heatmap(MV);
		H.setVisible(true);
	}

	
}
