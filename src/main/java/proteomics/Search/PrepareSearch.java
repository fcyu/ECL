package proteomics.Search;

import proteomics.Index.*;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import theoSeq.*;

public class PrepareSearch {

	private Map<String, ChainEntry> chain_entry_map = new HashMap<>();
	private BuildIndex build_index_obj = null;

	/////////////////////////////public methods///////////////////////////////////////////////////////////////////
	public PrepareSearch(Map<String, String> parameter_map) throws Exception {
		int nterm_linkable = Integer.valueOf(parameter_map.get("nterm_linkable"));
		int[] temp = string2Array(parameter_map.get("common_ion_charge"));
		int min_charge = temp[0];
		temp = string2Array(parameter_map.get("xlink_ion_charge"));
		Arrays.sort(temp);
		int max_charge = temp[temp.length - 1];
		Pattern link_site_pattern = Pattern.compile("K(?!$)");

		build_index_obj = new BuildIndex(parameter_map);
		MassTool mass_tool_obj = build_index_obj.returnMassTool();

		Map<String, Float> seq_mass_map = build_index_obj.returnSeqMassMap();
		Map<String, Set<String>> seq_pro_map = build_index_obj.returnSeqProMap();
		Map<String, Integer> seq_term_type_map = build_index_obj.returnSeqTermTypeMap();
		Map<String, Float> decoy_seq_mass_map = build_index_obj.returnDecoySeqMassMap();
		Map<String, String> decoy_seq_pro_map = build_index_obj.returnDecoySeqProMap();
		Map<String, Integer> decoy_seq_term_type_map = build_index_obj.returnDecoySeqTermTypeMap();

		// target
		for (String chain_seq : seq_mass_map.keySet()) {
			float chain_mass = seq_mass_map.get(chain_seq);
			Set<String> protein_id_set = seq_pro_map.get(chain_seq);
			int term_type = seq_term_type_map.get(chain_seq);
			float[][] chain_ion_array = mass_tool_obj.buildChainIonArray(chain_seq, min_charge, max_charge);

			// Generate a chain link site map for further use.
			Matcher link_site_matcher = link_site_pattern.matcher(chain_seq);
			List<Integer> link_site_list = new LinkedList<>();
			while (link_site_matcher.find()) {
				int link_site = link_site_matcher.start();
				link_site_list.add(link_site);
			}

			if ((term_type == 1) && !link_site_list.contains(0) && (nterm_linkable == 1)) {
				link_site_list.add(0);
			}

			Iterator<String> iterator_temp = protein_id_set.iterator();
			String protein_id = iterator_temp.next();
			while(iterator_temp.hasNext()) {
				protein_id += "&" + iterator_temp.next();
			}
			ChainEntry chain_entry = new ChainEntry(chain_mass, protein_id, "1", link_site_list, chain_ion_array);
			chain_entry_map.put(chain_seq, chain_entry);
		}

		// decoy
		for (String decoy_chain_seq : decoy_seq_mass_map.keySet()) {
			float decoy_chain_mass = decoy_seq_mass_map.get(decoy_chain_seq);
			String decoy_protein_id = decoy_seq_pro_map.get(decoy_chain_seq);
			int decoy_term_type = decoy_seq_term_type_map.get(decoy_chain_seq);
			float[][] decoy_chain_ion_array = mass_tool_obj.buildChainIonArray(decoy_chain_seq, min_charge, max_charge);

			// Generate a chain link site map for further use.
			Matcher link_site_matcher = link_site_pattern.matcher(decoy_chain_seq);
			List<Integer> decoy_link_site_list = new LinkedList<>();
			while (link_site_matcher.find()) {
				int link_site = link_site_matcher.start();
				decoy_link_site_list.add(link_site);
			}

			if ((decoy_term_type == 1) && !decoy_link_site_list.contains(0) && (nterm_linkable == 1)) {
				decoy_link_site_list.add(0);
			}

			ChainEntry chain_entry = new ChainEntry(decoy_chain_mass, decoy_protein_id, "0", decoy_link_site_list, decoy_chain_ion_array);
			chain_entry_map.put(decoy_chain_seq, chain_entry);
		}
	}

	public Map<String, ChainEntry> returnChainEntryMap() {
		return chain_entry_map;
	}

	public BuildIndex returnBuildIndex() {
		return build_index_obj;
	}

	private int[] string2Array(String str) {
		String[] temp = str.split(",");
		int[] output = new int[temp.length];
		for (int i = 0; i < output.length; ++i) {
			output[i] = Integer.valueOf(temp[i]);
		}

		return output;
	}
}
