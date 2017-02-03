package eem;

import java.applet.Applet;
import java.awt.BorderLayout;
import java.text.DecimalFormat;
import java.util.List;

import javax.swing.JFrame;
import utility.*;

public class test3   extends JFrame  {
	public test3(Applet C) throws Exception{
			setLayout(new BorderLayout());
			add(C,BorderLayout.CENTER);
			C.init();
			pack();
			setLocation(100, 100);
			setVisible(true);
	}
	public static void main(String[] args) throws Exception {
		
		StringMat A = new StringMat("test_annot.tab");
		KaplanMeierAnalysis KMA = new KaplanMeierAnalysis(A.getColMap("grade"));
		KMA.setStatus(A.getColMap("status"));
		KMA.setTime(A.getColMap("time"));
		KMA.perform();
		System.err.println("ok");
		new test3(KMA.getSurvivalCorvePlotViewer());
		
	}

}