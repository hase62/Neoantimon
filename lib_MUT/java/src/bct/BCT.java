package bct;

import java.util.*;
import java.io.*;
import utility.*;

public class BCT {
	int n;
	int m;
	boolean[][] T; 

	public BCT(int n, int m){
		this.n = n;
		this.m = m;
		
		T = new boolean[n][];
		for(int i=0;i<n;i++){
			T[i] = new boolean[m];
		}
		
		for(int i = 0; i<n; i++){
			for(int j = 0; j<m; j++){
					T[i][j] = false;
			}
		}
	}
	
	
	public BCT(MyMat M){
		n = M.rowSize();
		m = M.colSize();
		
		T = new boolean[n][];
		for(int i=0;i<n;i++){
			T[i] = new boolean[m];
		}
		
		for(int i = 0; i<n; i++){
			for(int j = 0; j<m; j++){
				if(M.get(i, j)!=0){
					T[i][j] = true;
				}else{
					T[i][j] = false;
				}
			}
		}
	}
	
	public BCT(BCT M){
		n = M.rowSize();
		m = M.colSize();
		
		T = new boolean[n][];
		for(int i=0;i<n;i++){
			T[i] = new boolean[m];
		}
		
		for(int i = 0; i<n; i++){
			for(int j = 0; j<m; j++){
				if(M.is1(i, j)){
					T[i][j] = true;
				}else{
					T[i][j] = false;
				}
			}
		}
	}
	
	public boolean get(int i, int j){
		return T[i][j];
	}
	
	public void set1(int i, int j){
		T[i][j] = true ;
	}
	
	public void set0(int i, int j){
		T[i][j] = false ;
	}
	
	public boolean is1(int i, int j){
		return T[i][j];
	}
	
	public boolean is0(int i, int j){
		return !T[i][j];
	}
	
	public int rowSize(){
		return n;
	}
	
	public int colSize(){
		return m;
	}
	
	public BCT getSubTable(List<Integer> I, List<Integer> J){
		BCT T = new BCT(I.size(), J.size());
		for(Integer i: I){
			for(Integer j: J){
				if(is1(i,j)){
					T.set1(i,j);
				}
			}
			
		}
		return T;
	}

	public int  sum(){
		int s = 0;
		for(int i = 0; i < n; i++){
			for(int j = 0; j<m; j++){
				if(is1(i, j)){
					s++;
				}
			}
		}
		return s;
	}
	
	public int  rowSum(int i){
		int s = 0;
		for(int j = 0; j<m; j++){
			if(is1(i, j)){
				s++;
			}
		}
		return s;
	}
	
	public int  colSum(int j){
		int s = 0;
		for(int i = 0; i<n; i++){
			if(is1(i, j)){
				s++;
			}
		}
		return s;
	}
	
	
	public String toString(){
		StringBuffer S = new StringBuffer();
		for(int i = 0; i<n; i++){
			for(int j = 0; j<m; j++){
				if(get(i, j)){
					S.append(1);
				}else{
					S.append(0);
				}
			}
			S.append("\n");
		}
		return S.toString();
	}

}



