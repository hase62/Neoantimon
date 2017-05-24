package utility;

public class Integrator {

	  private final float DAMPING = 0.3f;
	  private final float ATTRACTION = 0.6f;

	  private  float value;
	  private float vel;
	  private float accel;
	  private float force;
	  private float mass = 1;

	  private float damping = DAMPING;
	  private float attraction = ATTRACTION;
	  private boolean targeting;
	  private float target;


	  public Integrator() { }


	  public  Integrator(float value) {
	    this.value = value;
	  }


	  public Integrator(float value, float damping, float attraction) {
	    this.value = value;
	    this.damping = damping;
	    this.attraction = attraction;
	  }


	  public void set(float v) {
	    value = v;
	  }


	  public boolean  update() {
	    if (targeting) {
	      force += attraction * (target - value);      
	    }

	    accel = force / mass;
	    vel = (vel + accel) * damping;
	    value += vel;

	    force = 0;
	    if(Math.abs(target - value)  < target/10000){
	    	return false;
	    }else{
	    	return true;
	    }
	    
	  }


	  public void target(float t) {
	    targeting = true;
	    target = t;
	  }


	  public void noTarget() {
	    targeting = false;
	  }
	  
	  public float value(){
		  return value;
	  }
	  
	}