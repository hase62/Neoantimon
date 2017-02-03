package eem;

import java.sql.*;
import java.util.*;

import utility.*;


public class test4 {
	public static void main(String[] args) {
	    try {
	    	
	    MyMat E = Expression.getFromMySQL("test");	
	    System.out.print(ExpressionModuleSet.getFromTextFile("eem.out", E));
	      
	     //System.out.print(ExpressionModuleSet.getFromMySQL("test__test__eem"));
	     /* Connection con = MyFunc.getMySQLConnection("trans_fac");

	      Statement stmt = con.createStatement();
	      
	      
	      String sql = "SELECT * FROM PWM";
	      // クエリーを実行して結果セットを取得
	      ResultSet rs = stmt.executeQuery(sql);
	      // 検索された行数分ループ
	      rs.next();
	    
	        String id = rs.getString("id");
	        String mat = rs.getString("matrix");
	        System.out.println(id + "\t" + mat);
	      //}
	      
	      stmt.close();
	      con.close();*/
	    } catch (Exception e) {
	      e.printStackTrace();
	    }
	  }

	
	
	

}
