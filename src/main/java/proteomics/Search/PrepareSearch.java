package proteomics.Search;

import proteomics.Index.*;
import java.util.*;
import theoSeq.*;

public class PrepareSearch {

	private MassTool mass_tool_obj = null;
	private Map<String, String> parameter_map = null;
	private Map<String, String> index_chain_map = null;
	private Map<String, ChainEntry> chain_entry_map = null;
	private Map<String, Double> fix_mod_map = new HashMap<>();

	/////////////////////////////public methods///////////////////////////////////////////////////////////////////
	public PrepareSearch(Map<String, String> parameter_map) throws Exception {
		this.parameter_map = parameter_map;

		// Read fix modification
		fix_mod_map.put("G", Double.valueOf(parameter_map.get("G")));
		fix_mod_map.put("A", Double.valueOf(parameter_map.get("A")));
		fix_mod_map.put("S", Double.valueOf(parameter_map.get("S")));
		fix_mod_map.put("P", Double.valueOf(parameter_map.get("P")));
		fix_mod_map.put("V", Double.valueOf(parameter_map.get("V")));
		fix_mod_map.put("T", Double.valueOf(parameter_map.get("T")));
		fix_mod_map.put("C", Double.valueOf(parameter_map.get("C")));
		fix_mod_map.put("I", Double.valueOf(parameter_map.get("I")));
		fix_mod_map.put("L", Double.valueOf(parameter_map.get("L")));
		fix_mod_map.put("N", Double.valueOf(parameter_map.get("N")));
		fix_mod_map.put("D", Double.valueOf(parameter_map.get("D")));
		fix_mod_map.put("Q", Double.valueOf(parameter_map.get("Q")));
		fix_mod_map.put("K", Double.valueOf(parameter_map.get("K")));
		fix_mod_map.put("E", Double.valueOf(parameter_map.get("E")));
		fix_mod_map.put("M", Double.valueOf(parameter_map.get("M")));
		fix_mod_map.put("H", Double.valueOf(parameter_map.get("H")));
		fix_mod_map.put("F", Double.valueOf(parameter_map.get("F")));
		fix_mod_map.put("R", Double.valueOf(parameter_map.get("R")));
		fix_mod_map.put("Y", Double.valueOf(parameter_map.get("Y")));
		fix_mod_map.put("W", Double.valueOf(parameter_map.get("W")));
		fix_mod_map.put("Z", Double.valueOf(parameter_map.get("Z")));

		BuildIndex build_index_obj = new BuildIndex(parameter_map, fix_mod_map);
		mass_tool_obj = build_index_obj.returnMassToolObj();

		// Checking whether building a SQL database.
		System.out.println("Building database index...");
		build_index_obj.buildCLDb(parameter_map);

		BuildMap build_map_obj = new BuildMap(parameter_map, fix_mod_map);
		index_chain_map = build_map_obj.returnIndexChainMap();
		chain_entry_map = build_map_obj.returnChainEntryMap();
	}

	public MassTool returnMassTool() {
		return mass_tool_obj;
	}

	public Map<String, String> returnParameterMap() {
		return parameter_map;
	}

	public Map<String, String> returnIndexChainMap() {
		return index_chain_map;
	}

	public Map<String, ChainEntry> returnChainEntryMap() {
		return chain_entry_map;
	}

	public Map<String, Double> returnFixModMap() {
		return fix_mod_map;
	}
}
