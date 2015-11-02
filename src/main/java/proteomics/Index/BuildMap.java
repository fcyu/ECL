package proteomics.Index;

import theoSeq.MassTool;

import java.sql.*;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class BuildMap {

    private Map<String, ChainEntry> chain_entry_map = new HashMap<>();
    private Map<String, String> index_chain_map = new HashMap<>();

    public BuildMap(Map<String, String> parameter_map, Map<String, Double> fix_mod_map) throws Exception {
        int nterm_linkable = Integer.valueOf(parameter_map.get("nterm_linkable"));
        MassTool mass_tool_obj = new MassTool(Integer.valueOf(parameter_map.get("missed_cleavage")), fix_mod_map, nterm_linkable);
        int[] temp = string2Array(parameter_map.get("common_ion_charge"));
        int min_charge = temp[0];
        temp = string2Array(parameter_map.get("xlink_ion_charge"));
        int max_charge = temp[temp.length - 1];

        Pattern link_site_pattern = Pattern.compile(parameter_map.get("cl_aa") + "(?!$)");

        // Create a index SQL connection and read the database from disk.
        Connection index_conn = null;
        Statement index_statement = null;
        ResultSet index_rs = null;
        try {
            Class.forName("org.sqlite.JDBC").newInstance();
            index_conn = DriverManager.getConnection("jdbc:sqlite::memory:");
            index_statement = index_conn.createStatement();
            index_statement.executeUpdate("restore from " + parameter_map.get("db") + ".db");
        } catch (SQLException ex) {
            System.err.println(ex.getMessage());
            System.exit(1);
        }

        // Extract chain sequence and mass from peptide_chain table for further use
        index_rs = index_statement.executeQuery("select chain_index, mass, peptide_sequence, protein_id, type, term_type from peptide_chain;");
        while (index_rs.next()) {
            String chain_index = String.valueOf(index_rs.getInt(1)) + "-"; //TODO: Left for vairable modification
            float chain_mass = index_rs.getFloat(2);
            String chain_seq = index_rs.getString(3);
            String protein_id = index_rs.getString(4);
            String chain_type = index_rs.getString(5);
            int term_type = index_rs.getInt(6);
            double[][] chain_ion_array = mass_tool_obj.buildChainIonArray(chain_seq, min_charge, max_charge);

            index_chain_map.put(chain_index, chain_seq);

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

            ChainEntry chain_entry = new ChainEntry(chain_mass, chain_index, protein_id, chain_type, link_site_list, chain_ion_array);
            chain_entry_map.put(chain_seq, chain_entry);
        }

        index_statement.close();
        index_conn.close();
    }

    public Map<String, ChainEntry> returnChainEntryMap() {
        return chain_entry_map;
    }

    public Map<String, String> returnIndexChainMap() {
        return index_chain_map;
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
