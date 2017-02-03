package eem;

import utility.*;

import java.applet.Applet;
import java.awt.BorderLayout;
import javax.swing.JFrame;


public class SeedGeneProfileVeiwer   extends JFrame  {

	
	
	private static final long serialVersionUID = 2913512603921292733L;



	public SeedGeneProfileVeiwer (Applet C) throws Exception{
		setLayout(new BorderLayout());
		add(C,BorderLayout.CENTER);
		C.init();
		pack();
		setLocation(100, 100);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		
	}
	

	public static void main(String[] args) throws Exception {
		//ExpressionModuleSet ems = ExpressionModuleSet.getFromTextFile(args[0], args[1]);
		ExpressionModuleSet ems = ExpressionModuleSet.getFromMySQL(args[0]);
		ExpressionModule em = ems.get(args[1]);
		em.addBiclusterInformation2SampleAnnotation();
		ClusteredMyMatViewer C = new ClusteredMyMatViewer(em.getClusteredSeedGeneProfiles());
		C.scaleColorByRow();
		SeedGeneProfileVeiwer V =  new SeedGeneProfileVeiwer(C);
		V.setVisible(true);
		
	}
	
	
	
}
