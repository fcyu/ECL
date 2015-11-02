package proteomics.Index;

import java.util.*;
import java.sql.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import theoSeq.*;
import proteomics.Validation.*;

public class BuildIndex {

	private static final String TARGET = "1";
	private static final String DECOY = "0";

	private double hatom_mass = 0;
	private double oatom_mass = 0;
	private float max_precursor_mass = 0;
	private int min_chain_length = 0;
	private int max_chain_length = 0;
	private String db_path = "";
	private MassTool mass_tool_obj = null;
	private Map<String, String> pro_seq_map = null;
	final private Map<String, Float> seq_mass_map = new HashMap<>();
	final private Map<String, Set<String>> seq_pro_map = new HashMap<>();
	final private Set<String> for_check_duplitate = new HashSet<>();
	final private Map<String, Integer> seq_term_type_map = new HashMap<>();
	final private Map<String, Float> decoy_seq_mass_map = new HashMap<>();
	final private Map<String, Set<String>> decoy_seq_pro_map = new HashMap<>();
	final private Map<String, Integer> decoy_seq_term_type_map = new HashMap<>();
	private double nterm_mass = 0;

	/////////////////////////////////public methods//////////////////////////////////////////////////////////////////
	public BuildIndex(Map<String, String> parameter_map, Map<String, Double> fix_mod_map) throws Exception {
		// initialize parameters
		max_precursor_mass = Float.valueOf(parameter_map.get("max_precursor_mass"));
		min_chain_length = Integer.valueOf(parameter_map.get("min_chain_length"));
		max_chain_length = Integer.valueOf(parameter_map.get("max_chain_length"));
		String decoy_db = parameter_map.get("decoy_db");
		db_path = parameter_map.get("db");
		int missed_cleavage = Integer.valueOf(parameter_map.get("missed_cleavage"));
		String cl_aa = parameter_map.get("cl_aa");
		int nterm_linkable = Integer.valueOf(parameter_map.get("nterm_linkable"));

		// read protein database
		DbTool db_tool_obj = new DbTool(db_path);
		pro_seq_map = db_tool_obj.returnSeqMap();

		// define a new MassTool object
		mass_tool_obj = new MassTool(missed_cleavage, fix_mod_map, nterm_linkable);
		Map<String, Double> mass_table = mass_tool_obj.returnMassTable();
		hatom_mass = mass_table.get("Hatom");
		oatom_mass = mass_table.get("Oatom");
		nterm_mass = mass_table.get("Z");
	}

	/////////////////////////////////////public methods////////////////////////////////////////////////////////////////////
	public void buildCLDb(Map<String, String> parameter_map) throws Exception {
		buildPepChainMap();
		buildDecoyPepChainMap();

		// Create index SQL database
		Connection index_conn = null;
		Statement index_statement = null;
		try {
			Class.forName("org.sqlite.JDBC").newInstance();
			index_conn = DriverManager.getConnection("jdbc:sqlite:" + ":memory:");
			index_statement = index_conn.createStatement();
			index_statement.executeUpdate("create table peptide_chain (chain_index integer primary key asc, mass real not null, peptide_sequence text not null, protein_id text not null, type text not null, term_type integer not null);");
			index_statement.executeUpdate("create index mass on peptide_chain (mass);");
			index_statement.executeUpdate("create index seq on peptide_chain (peptide_sequence);");
		} catch (SQLException ex) {
			System.err.println("SQLException: " + ex.getMessage());
			System.exit(1);
		}

		// Build peptide_chain table
		// It contains target and decoy sequence
		PreparedStatement index_prepared_statement = index_conn.prepareStatement("insert into peptide_chain (mass, peptide_sequence, protein_id, type, term_type) values (?, ?, ?, ?, ?);");
		try {
			index_conn.setAutoCommit(false);
			Set<String> seq_set = seq_pro_map.keySet();
			for (String seq : seq_set) {
				if ((seq.length() < min_chain_length) || (seq.length() > max_chain_length)) {
					continue;
				}

				index_prepared_statement.setFloat(1, seq_mass_map.get(seq));
				index_prepared_statement.setString(2, seq);
				Set<String> pro_id_set = seq_pro_map.get(seq);
				String temp = pro_id_set.toString().replace(", ", "&");
				String temp_2 = temp.substring(1, temp.length() - 1);
				index_prepared_statement.setString(3, temp_2);
				index_prepared_statement.setString(4, TARGET);
				index_prepared_statement.setInt(5, seq_term_type_map.get(seq));
				index_prepared_statement.executeUpdate();
			}
			index_conn.commit();
		} catch (SQLException ex) {
			System.err.println("SQLException: " + ex.getMessage());
			System.exit(1);
		} finally {
			index_prepared_statement.close();
			index_conn.setAutoCommit(true);
		}

		index_prepared_statement = index_conn.prepareStatement("insert into peptide_chain (mass, peptide_sequence, protein_id, type, term_type) values (?, ?, ?, ?, ?);");
		try {
			index_conn.setAutoCommit(false);
			Set<String> seq_set = decoy_seq_pro_map.keySet();
			for (String seq : seq_set) {
				if ((seq.length() < min_chain_length) || (seq.length() > max_chain_length)) {
					continue;
				}

				index_prepared_statement.setFloat(1, decoy_seq_mass_map.get(seq));
				index_prepared_statement.setString(2, seq);
				Set<String> pro_id_set = decoy_seq_pro_map.get(seq);
				String temp = pro_id_set.toString().replace(", ", "&");
				String temp_2 = temp.substring(1, temp.length() - 1);
				index_prepared_statement.setString(3, temp_2);
				index_prepared_statement.setString(4, DECOY);
				index_prepared_statement.setInt(5, decoy_seq_term_type_map.get(seq));
				index_prepared_statement.executeUpdate();
			}
			index_conn.commit();
		} catch (SQLException ex) {
			System.err.println("SQLException: " + ex.getMessage());
			System.exit(1);
		} finally {
			index_prepared_statement.close();
			index_conn.setAutoCommit(true);
		}

		// Backup SQL database to disk
		try {
			index_statement.executeUpdate("backup to " + db_path + ".db");
		} catch (SQLException ex) {
			System.err.println("SQLException: " + ex.getMessage());
			System.exit(1);
		}

		index_statement.close();
		index_conn.close();
	}

	public MassTool returnMassToolObj() {
		return mass_tool_obj;
	}

	//////////////////////////////////////////private methods////////////////////////////////////////////////////////
	private void buildPepChainMap() {
		Set<String> pro_id_set = pro_seq_map.keySet();
		for (String pro_id : pro_id_set) {
			String pro_seq = pro_seq_map.get(pro_id);
			Set<String> seq_set = mass_tool_obj.buildChainSet(pro_seq);
			for (String seq : seq_set) {
				float mass_temp = (float) mass_tool_obj.calResidueMass(seq) + (float) nterm_mass + 2 * (float) hatom_mass + (float) oatom_mass; // calMass just calculate the residue mass, so we should add a H2O
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
						Set<String> pro_list = seq_pro_map.get(seq);
						pro_list.add(pro_id);
						seq_pro_map.put(seq, pro_list);
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
			String decoy_seq = "";
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

			Set<String> decoy_pro_set = new HashSet<>();
			decoy_pro_set.add(decoy_pro_id);
			decoy_seq_pro_map.put(decoy_seq, decoy_pro_set);
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
