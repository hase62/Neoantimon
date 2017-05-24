package network;

import java.util.*;

public class NullLinkGenerator {
	Link L;
	List <String> name;
	Map <String, Integer> name2degree;
	Map <Integer, List<String>> degree2name;
	
	public NullLinkGenerator(Link L){
		this.L = L;
		name = L.getNodeName();
		setDegree();
	}
	
	private void setDegree(){
		name2degree = new HashMap <String, Integer>();
		degree2name = new HashMap <Integer, List<String>>();
		for(String n1: name){
			int i=0;
			for(String n2: name){
				if(L.get(n1, n2)){
					i++;
				}
			}
			name2degree.put(n1, i);
			if(!degree2name.containsKey(i)){
				List <String> tmp = new ArrayList<String>();
				tmp.add(n1);
				degree2name.put(i, tmp);
			}else{
				degree2name.get(i).add(n1);
			}	
		}
	}
	
	public Link getRondomNetwork(){
		Link rL = new Link(L);	
		List <String> rname =  new  ArrayList<String>(name);
		for(Integer d: degree2name.keySet()){
			List<String> tmp = new ArrayList<String>(degree2name.get(d));
			Collections.shuffle(tmp);
			for(int i = 0; i < tmp.size(); i++){
				rname.set(L.getName2index().get(degree2name.get(d).get(i)),tmp.get(i));			
			}
		}
		rL.setNodeName(rname);
		return rL;
	}
	
	public Link getRondomNetworks(List<String> target){
		Link rL = new Link(L);	
		List <String> rname =  new  ArrayList<String>(name);
		Random rnd = new Random();
		for(String s: target){
			int d = name2degree.get(s);
			List<String> tmp = new ArrayList<String>(degree2name.get(d));
			int j = L.getName2index().get(s);
			int  i = rnd.nextInt(tmp.size());
			rname.set(j, rname.get(i));
			rname.set(i, s);
		}
		rL.setNodeName(rname);
		return rL;
	}
	

}
