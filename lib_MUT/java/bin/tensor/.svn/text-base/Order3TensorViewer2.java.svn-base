package tensor;
import processing.core.*;
import peasy.*;
import utility.*;
import controlP5.*;

public class Order3TensorViewer2 extends PApplet{
	
	
	private static final long serialVersionUID = 1L;


	protected Order3Tensor T;  


	ControlP5 controlP5;
	PeasyCam cam;

	protected float xSize;
	protected int xCells;
	protected float xCellSize;

	protected float ySize;
	protected int yCells;
	protected float yCellSize;

	protected float zSize;
	protected int zCells;
	protected float zCellSize;

	protected int color1 = MyColor.getRGBhex("BLUE");
	protected int color2 = MyColor.getRGBhex("WHITE");
	protected int color3 = MyColor.getRGBhex("RED");
	protected int[][][] colors;
	// camera overlay
	protected PMatrix3D currCameraMatrix;
	protected PGraphics3D g3; 


	protected int boxFill;
	protected int countX = 0;
	protected int countY = 0;
	protected int countZ = 0;


	protected int myColor = color(255,0,0);
	protected Slider xS,yS,zS;

	public Order3TensorViewer2(Order3Tensor T){
		this.T = T;
		setSize(300,300,300);
		setColor();	
	}	
	
	void setColor() {
		float max = (float)MyFunc.max(T.asList());
		float min = (float)MyFunc.min(T.asList());
		colors = new int[xCells][yCells][zCells];
		for (int x=0;x<xCells;x++) {
			for (int y=0;y<yCells;y++) {
				for (int z=0;z<zCells;z++) {
						float p = (float) ((T.get(x, y, z)-min)/(max-min));
						int color = PApplet.lerpColor(color2,color3, p, RGB);
						//int color = MyColor.lerpColor(color1,color2,color3, p, RGB);
						colors[x][y][z] = color;
		      }
		    }
		  }
		}

	
	
	
	void setSize(int xSize, int ySize, int zSize){
		xCells = T.getDimOfOrder1();
		yCells = T.getDimOfOrder2(); 
		zCells = T.getDimOfOrder3();
		this.xSize  = xSize;	
		this.ySize  = ySize;
		this.zSize  = zSize;
		xCellSize = xSize/(float)xCells;
		yCellSize = ySize/(float)yCells;
		zCellSize = zSize/zCells;
		println(xCells + " " + yCells + " " + zCells + " " + xCellSize + " " + yCellSize + " " + zCellSize);
		
		
	}	
	public void setup() {

	  size(800, 800, OPENGL);
	  //hint(ENABLE_DEPTH_SORT);
	  g3 = (PGraphics3D)g;
	  
	  controlP5 = new ControlP5(this);

		
	  lights();
	  noStroke();
	  stroke(255);
	  
	 // float fov = (float) (4*PI/3.0); 
	  //float cameraZ = (float) ((height/2.0) / tan((float) (PI * fov / 360.0))); 
	  //perspective(fov,(int)(width/height), (float)(cameraZ/2.0), (float)(cameraZ*2.0)); 
	  
	  // orthographic view
	  ortho(-width/2, width/2, -height/2, height/2, -1000, 1000); 

	  // PeasyCam Setup
	  cam = new PeasyCam(this, 750);
	  cam.setMinimumDistance(-500);
	  cam.setMaximumDistance(1500);

	  xS = (Slider)controlP5.addSlider("xSlider",0,xCells,0,10,10,200,20);
	  xS.setColorActive(color(0));
	  xS.setColorBackground( color(200,100) );
	  xS.setColorForeground(color(150,100));
	  xS.setColorLabel(color(0));
	  xS.setColorValue(color(127));
	  xS.setLabel("X");
	  xS.setLabelVisible(true);
	  xS.setNumberOfTickMarks(xCells+1);

	  
	  yS = (Slider)controlP5.addSlider("ySlider",0,yCells,0,10,30,200,20);
	  yS.setColorActive(color(0));
	  yS.setColorBackground( color(200,100) );
	  yS.setColorForeground(color(150,100));
	  yS.setColorLabel(color(0));
	  yS.setColorValue(color(127));
	  yS.setLabel("Y");
	  yS.setLabelVisible(true);
	  yS.setNumberOfTickMarks(yCells+1);
	  
	  zS = (Slider)controlP5.addSlider("zSlider",0,zCells,0,10,50,200,20);
	  zS.setColorActive(color(0));
	  zS.setColorBackground( color(200,100) );
	  zS.setColorForeground(color(150,100));
	  zS.setColorLabel(color(0));
	  zS.setColorValue(color(127));
	  zS.setLabel("Z");
	  zS.setLabelVisible(true);
	  zS.setNumberOfTickMarks(zCells+1);
	  
	  controlP5.setAutoDraw(false);
	  
	}

	

	public void draw() {
	 
	  background(255);
	  //noLoop();
	  
	  // Center and spin grid
	 /* translate(width/2, height/2, -200);
	   rotateY((float) (frameCount * 0.01));
	  rotateX((float) (frameCount * 0.01));*/
	   
	

	  for (int i=0; i<(xCells-xS.value()); i++) {
	    for (int j=0; j<(yCells-yS.value()); j++) {
	      for (int k=0; k<(zCells-zS.value()); k++) {
	        pushMatrix();


	        translate( (-xSize/2)+i*xCellSize, (-ySize/2)+j*yCellSize, (-zSize/2)+k*zCellSize );
	        
	       
	        boxFill = colors[i][j][k];

	        if (alpha(boxFill) > 0) {
	          fill(boxFill);
	          box(xCellSize, yCellSize, zCellSize);
	        }
	        popMatrix();
	      }
	    }
	  }
	  
	  
	  cam.beginHUD();
	  gui();
	  cam.endHUD();
	  
	  cam.setMouseControlled(true);
	  if(xS.isInside() || yS.isInside() || zS.isInside() ) {
	    cam.setMouseControlled(false);
	  } 
	  
	}

	void gui() {
	  currCameraMatrix = new PMatrix3D(g3.camera);
	  camera();
	  controlP5.draw();
	  g3.camera = currCameraMatrix;
	}

	


	
}
