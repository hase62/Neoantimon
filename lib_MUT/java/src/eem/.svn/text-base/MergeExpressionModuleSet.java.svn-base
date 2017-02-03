package eem;

public class MergeExpressionModuleSet {
	public static void main(String [] args) throws Exception{
		ExpressionModuleSet ems = new ExpressionModuleSet();
		for(int i = 0; i < args.length -1; i++){
			ems.addAll(ExpressionModuleSet.ReadFromFile(args[i]));	
		}
		ems.writeToFile(args[args.length-1]);
	}
}
