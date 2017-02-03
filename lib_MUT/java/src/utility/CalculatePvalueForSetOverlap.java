package utility;

import sun.reflect.Reflection;

public class CalculatePvalueForSetOverlap {
	public static void main(String [] args) throws Exception{
		if(args.length != 4){
			System.err.println(Reflection.getCallerClass( 1 ).getName() + " bg_size set1_size set2_size isect_size");
			return;
		}
		System.out.println(MyFunc.calculatePvalueForSetOverlap(Integer.valueOf(args[0]),Integer.valueOf(args[1]),Integer.valueOf(args[2]),Integer.valueOf(args[3])));
		
	}



}
