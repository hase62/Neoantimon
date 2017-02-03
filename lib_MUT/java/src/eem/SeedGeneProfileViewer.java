package eem;

import sun.reflect.Reflection;
import utility.*;

import java.applet.Applet;
import java.awt.BorderLayout;
import javax.swing.JFrame;


public class SeedGeneProfileViewer   extends JFrame  {

	
	
	private static final long serialVersionUID = 2913512603921292733L;



	public SeedGeneProfileViewer (Applet C) throws Exception{
		setLayout(new BorderLayout());
		add(C,BorderLayout.CENTER);
		C.init();
		pack();
		setLocation(100, 100);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		
	}
	

	public static void main(String[] args) throws Exception {
		//ExpressionModuleSet ems = ExpressionModuleSet.getFromTextFile(args[0], args[1]);
		if(args.length < 2){
			System.err.println(Reflection.getCallerClass( 1 ).getName() + " table_name id");
			return;
		}
		ExpressionModuleSet ems = ExpressionModuleSet.getFromMySQL(args[0]);
		ExpressionModule em = ems.get(args[1]);
		em.addBiclusterInformation2SampleAnnotation();
		ClusteredMyMatViewer C = new ClusteredMyMatViewer(em.getClusteredSeedGeneProfiles());
		C.scaleColorByRow();
		SeedGeneProfileViewer V =  new SeedGeneProfileViewer(C);
		V.setVisible(true);
		
	}
	
	
	
}
