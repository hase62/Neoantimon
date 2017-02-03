package utility;

import java.io.Serializable;
import java.util.*;

import eem.EEM;

public class StopWatch implements Serializable{
	private Date startDate;
	private Date stopDate;
	private int miliSecond;	
	private int second;
	private int minute;
	private int hour;
	private double elapsedTime;
	public StopWatch(){
		startDate = new Date();
		stopDate = null;
	}
	public void stop(){
		stopDate = new Date();
		elapsedTime = stopDate.getTime() - startDate.getTime();
		second = (int) elapsedTime/1000;
		miliSecond = (int)elapsedTime - second * 1000;
		if(second < 60){
			minute = 0;
			hour = 0;
			return;
		}
		minute = second/60;
		second = second%60;
		if(minute < 60){
			hour = 0;
			return;
		}
		hour = minute/60;
		minute = minute%60;
		return;
	}
	public String  toString(){
		if(stopDate == null){
			return null;
		}
		if(hour != 0){
			return hour + "hr " + minute + "min " + second + "sec";
		}
		if(minute !=0 ){
			return  minute + "min " + second + "sec";
		}
		return second + "." + miliSecond  + "sec";
	}
	public Date getStartDate(){ return startDate;}
	public Date getStopDate(){ return stopDate;}
	public double getElapsedTime(){ return elapsedTime; }
}
