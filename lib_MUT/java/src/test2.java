
import java.applet.Applet;
import java.awt.Graphics;

/*
<applet code="tensor/test2.class" width="150" height="150">
</applet>
*/

public class test2 extends Applet{
  public void paint(Graphics g){
    int xPoints[] = {10, 90, 75, 20};
    int yPoints[] = {20, 25, 100, 90};

    g.drawPolygon(xPoints, yPoints, 4);
  }
}