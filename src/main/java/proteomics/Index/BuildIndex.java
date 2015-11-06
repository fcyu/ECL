package proteomics.Index;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import theoSeq.*;

public class BuildIndex {

	private float hatom_mass = 0;
	private float oatom_mass = 0;
	private float max_precursor_mass = 0;
	private int min_chain_length = 0;
	private int max_chain_length = 0;
	private MassTool mass_tool_obj = null;
	private Map<String, String> pro_seq_map = null;
	private Map<String, Float> seq_mass_map = new HashMap<>();
	private Map<String, Set<String>> seq_pro_map = new HashMap<>();
	private Set<String> for_check_duplitate = new HashSet<>();
	private Map<String, Integer> seq_term_type_map = new HashMap<>();
	private Map<String, Float> decoy_seq_mass_map = new HashMap<>();
	private Map<String, String> decoy_seq_pro_map = new HashMap<>();
	private Map<String, Integer> decoy_seq_term_type_map = new HashMap<>();
	private Map<String, Float> fix_mod_map = new HashMap<>();
	private float nterm_mass = 0;

	/////////////////////////////////public methods//////////////////////////////////////////////////////////////////
	public BuildIndex(Map<String, String> parameter_map) throws Exception {
		// initialize parameters
		max_precursor_mass = Float.valueOf(parameter_map.get("max_precursor_mass"));
		min_chain_length = Integer.valueOf(parameter_map.get("min_chain_length"));
		max_chain_length = Integer.valueOf(parameter_map.get("max_chain_length"));
		String db_path = parameter_map.get("db");
		int missed_cleavage = Integer.valueOf(parameter_map.get("missed_cleavage"));
		int nterm_linkable = Integer.valueOf(parameter_map.get("nterm_linkable"));

		// Read fix modification
		fix_mod_map.put("G", Float.valueOf(parameter_map.get("G")));
		fix_mod_map.put("A", Float.valueOf(parameter_map.get("A")));
		fix_mod_map.put("S", Float.valueOf(parameter_map.get("S")));
		fix_mod_map.put("P", Float.valueOf(parameter_map.get("P")));
		fix_mod_map.put("V", Float.valueOf(parameter_map.get("V")));
		fix_mod_map.put("T", Float.valueOf(parameter_map.get("T")));
		fix_mod_map.put("C", Float.valueOf(parameter_map.get("C")));
		fix_mod_map.put("I", Float.valueOf(parameter_map.get("I")));
		fix_mod_map.put("L", Float.valueOf(parameter_map.get("L")));
		fix_mod_map.put("N", Float.valueOf(parameter_map.get("N")));
		fix_mod_map.put("D", Float.valueOf(parameter_map.get("D")));
		fix_mod_map.put("Q", Float.valueOf(parameter_map.get("Q")));
		fix_mod_map.put("K", Float.valueOf(parameter_map.get("K")));
		fix_mod_map.put("E", Float.valueOf(parameter_map.get("E")));
		fix_mod_map.put("M", Float.valueOf(parameter_map.get("M")));
		fix_mod_map.put("H", Float.valueOf(parameter_map.get("H")));
		fix_mod_map.put("F", Float.valueOf(parameter_map.get("F")));
		fix_mod_map.put("R", Float.valueOf(parameter_map.get("R")));
		fix_mod_map.put("Y", Float.valueOf(parameter_map.get("Y")));
		fix_mod_map.put("W", Float.valueOf(parameter_map.get("W")));
		fix_mod_map.put("Z", Float.valueOf(parameter_map.get("Z")));

		// read protein database
		DbTool db_tool_obj = new DbTool(db_path);
		pro_seq_map = db_tool_obj.returnSeqMap();

		// define a new MassTool object
		mass_tool_obj = new MassTool(missed_cleavage, fix_mod_map, nterm_linkable);
		Map<String, Float> mass_table = mass_tool_obj.returnMassTable();
		hatom_mass = mass_table.get("Hatom");
		oatom_mass = mass_table.get("Oatom");
		nterm_mass = mass_table.get("Z");

		buildPepChainMap();
		buildDecoyPepChainMap();
	}

	/////////////////////////////////////public methods////////////////////////////////////////////////////////////////////
	public MassTool returnMassTool() {
		return mass_tool_obj;
	}

	public Map<String, String> returnProSeqMap() {
		return pro_seq_map;
	}

	public Map<String, Float> returnSeqMassMap() {
		return seq_mass_map;
	}

	public Map<String, Set<String>> returnSeqProMap() {
		return seq_pro_map;
	}

	public Map<String, Integer> returnSeqTermTypeMap() {
		return seq_term_type_map;
	}

	public Map<String, Float> returnDecoySeqMassMap() {
		return decoy_seq_mass_map;
	}

	public Map<String, String> returnDecoySeqProMap() {
		return decoy_seq_pro_map;
	}

	public Map<String, Integer> returnDecoySeqTermTypeMap() {
		return decoy_seq_term_type_map;
	}

	public Map<String, Float> returnFixModMap() {
		return fix_mod_map;
	}

	//////////////////////////////////////////private methods////////////////////////////////////////////////////////
	private void buildPepChainMap() {
		Set<String> pro_id_set = pro_seq_map.keySet();
		for (String pro_id : pro_id_set) {
			String pro_seq = pro_seq_map.get(pro_id);
			Set<String> seq_set = mass_tool_obj.buildChainSet(pro_seq);
			for (String seq : seq_set) {
				if ((seq.length() < min_chain_length) || (seq.length() > max_chain_length)) {
					continue;
				}

				float mass_temp = mass_tool_obj.calResidueMass(seq) + nterm_mass + 2 * hatom_mass + oatom_mass;
				if (mass_temp <= max_precursor_mass) {
					seq_mass_map.put(seq, mass_temp);

					// Add the sequence to the check set for decoy duplicate check
					String template_seq = seq.replace("L", "I"); // "L" and "I" have the same mass.
					template_seq = template_seq.replace("K", "Q"); // "K" and "Q" have the close mass.
					for_check_duplitate.add(template_seq);

					if (pro_seq.startsWith(seq)) {
						seq_term_type_map.put(seq, 1);
					} else if (pro_seq.endsWith(seq)) {
						seq_term_type_map.put(seq, 2);
					} else {
						seq_term_type_map.put(seq, 0);
					}

					if (seq_pro_map.containsKey(seq)) {
						seq_pro_map.get(seq).add(pro_id);
					} else {
						Set<String> pro_list = new HashSet<>();
						pro_list.add(pro_id);
						seq_pro_map.put(seq, pro_list);
					}
				}
			}
		}
	}

	private void buildDecoyPepChainMap() throws Exception {
		Set<String> seq_set = seq_pro_map.keySet();
		for (String original_seq : seq_set) {
			String decoy_seq;
			if (original_seq.endsWith("K") || original_seq.endsWith("R")) {
				String temp = reverseSeq(original_seq.substring(0, original_seq.length() - 1));
				decoy_seq = temp + original_seq.charAt(original_seq.length() - 1);
			} else {
				decoy_seq = reverseSeq(original_seq);
			}

			// Check duplicate
			String new_decoy_seq = decoy_seq.replace("L", "I");
			new_decoy_seq = new_decoy_seq.replace("K", "Q");
			if (for_check_duplitate.contains(new_decoy_seq)) {
				// the decoy sequence is the same as the target sequence
				continue;
			}

			float decoy_mass_temp = seq_mass_map.get(original_seq);
			decoy_seq_mass_map.put(decoy_seq, decoy_mass_temp);
			String pro_id = seq_pro_map.get(original_seq).iterator().next();
			String decoy_pro_id = "DECOY_" + pro_id;
			String pro_seq = pro_seq_map.get(pro_id);
			if (pro_seq.startsWith(original_seq)) {
				decoy_seq_term_type_map.put(decoy_seq, 1);
			} else if (pro_seq.endsWith(original_seq)) {
				decoy_seq_term_type_map.put(decoy_seq, 2);
			} else {
				decoy_seq_term_type_map.put(decoy_seq, 0);
			}

			decoy_seq_pro_map.put(decoy_seq, decoy_pro_id);
		}
	}

	private String reverseSeq(String sequence) {
		Pattern fix_pattern = Pattern.compile("[K]");
		String decoy_str = "";
		Matcher fix_matcher = fix_pattern.matcher(sequence);
		int sequence_length = sequence.length();
		int idx_1 = 0;
		int idx_2;
		while (idx_1 < sequence_length) {
			String fix_aa;
			if (fix_matcher.find()) {
				idx_2 = fix_matcher.start();
				fix_aa = sequence.substring(idx_2, idx_2 + 1);
			} else {
				idx_2 = sequence_length;
				fix_aa = "";
			}
			String part = sequence.substring(idx_1, idx_2);

			// Reverse part sequence
			decoy_str += new StringBuilder(part).reverse().toString() + fix_aa;
			idx_1 = idx_2 + 1;
		}

		return decoy_str;
	}
}
