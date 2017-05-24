package snp;

import java.io.*;
import java.util.*;

public class Segment {

	private int start;
	private int end;
	private int chr;
	private double value;
	private int probeCount;
	
	
	public Segment(int chr, int start, int end, int probeCount, double value){
		this.chr = chr;
		this.start = start;
		this.end = end;
		this.probeCount = probeCount;
		this.value = value;
	}
	
	public Segment(){};
	
	public Segment(Segment S){
		this.chr = S.chr;
		this.start = S.start;
		this.end = S.end;
		this.probeCount = S.probeCount;
		this.value = S.value;
	}
	
	public double value(){
		return value;
	}
	
	public int chr(){
		return chr;
	}
	
	public int start(){
		return start;
	}
	
	public void start(int i){
		start = i;
	}
	
	public int end(){
		return end;
	}
	
	public void end(int i){
		end = i;
	}
	
	public int probeCount(){
		return probeCount;
	}
	
	public void probeCount(int i){
		probeCount = i;
	}
	
	public void value(double d){
		value = d;
	}
	
	public String toString(){
		return chr + "\t" + start + "\t" + end + "\t" + probeCount + "\t" + value;
	}
	
	

}
