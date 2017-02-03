package utility;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.*;

public class VariableDist implements Serializable{
	private static final long serialVersionUID = -2840718245427643884L;
	private Map < String, Map <String, Double>> M;
	private Set<String> name;
	public List<String> getNames(){
		return new ArrayList<String>(name);
	}
	public boolean containsName(String name){
		return name.contains(name);
	}
	public double get(String s, String t){
		if(s.compareTo(t)> 0){
			return M.get(s).get(t);
		}
		if(s.compareTo(t) < 0){
			return M.get(t).get(s);
		}
		return 0;
	}
	
	public VariableDist(Dist D){
		name = new HashSet<String>(D.getNames());
		M = new HashMap < String, Map <String, Double>>();
		for(String s: name){
			Map <String, Double> tmp = new HashMap<String, Double>();
			for(String t: name){
				if(s.compareTo(t)> 0){
					tmp.put(t, D.get(s,t));
				}
			}
			M.put(s, tmp);
		}
	}
	public VariableDist(List <String> name){
		this.name = new HashSet<String>(name);
		M = new HashMap < String, Map <String, Double>>();
		for(String s: name){
			Map <String, Double> tmp = new HashMap<String, Double>();
			for(String t: name){
				if(s.compareTo(t)> 0){
					tmp.put(t, null);
				}
			}
			M.put(s, tmp);
		}
	}
	
	public VariableDist getSubVariableDist(List <String> names){
		VariableDist D = new VariableDist(names);
		for(String s: names){
			for(String t: names){
				if(s.compareTo(t) > 0){
					D.set(s,t,get(s,t));
				}
			}
		}
		return D;
	}
	
	
	public void set(String s, String t, double d){
		if(s.compareTo(t)> 0){
			M.get(s).put(t, d);
		}
		if(s.compareTo(t) < 0){
			M.get(t).put(s, d);
		}
	}
	
	public void remove(String s){
		for(String t : name){
			if(s.compareTo(t) < 0){
				M.get(t).remove(s);
			}
		}
		M.remove(s);
		name.remove(s);
	}
	public void add(String s){
		Map <String, Double> tmp = new HashMap<String, Double>(); 
		for(String t :name){
			if(s.compareTo(t) < 0){
				M.get(t).put(s,null);
			}
			if(s.compareTo(t) > 0){
				tmp.put(t, null);
			}
		}
		M.put(s, tmp);
		name.add(s);
	}
	
	public void add(String s, Map <String, Double> m){
		add(s);
		for(Map.Entry<String, Double> e: m.entrySet()){
			if(name.contains(e.getKey())){
				set(s, e.getKey(), e.getValue());
			}
		}
	}
	public void print(String outfile) throws IOException {
		PrintWriter os = new PrintWriter(new FileWriter(outfile));
		for(String s: name){
			for(String t: name){
				if(s.compareTo(t)> 0){
				os.println(s + "\t" + t + "\t" + get(s,t));
				}
			}
		}
		os.close();
	}
	
	public void print() throws IOException {
		PrintWriter os = new PrintWriter(System.out);
		for(String s: name){
			for(String t: name){
				if(s.compareTo(t)> 0){
				os.println(s + "\t" + t + "\t" + get(s,t));
				}
			}
		}
		os.close();
	}
}
