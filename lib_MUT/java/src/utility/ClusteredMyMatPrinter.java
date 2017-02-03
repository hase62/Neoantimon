package utility;

public class ClusteredMyMatPrinter extends  ClusteredMyMatViewer{

	private String outFile;	
	public ClusteredMyMatPrinter(ClusteredMyMat M, String outFile) {
		super(M);
		this.outFile = outFile;
	}
	
	public void setup(){
		size(Math.round(Width), Math.round(Height), PDF, "out.pdf");
		textFont(font);
		noLoop();
	}
	
	public void print(){
		init();
		setup();
		draw();
	}

}
