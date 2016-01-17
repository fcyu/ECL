package theoSeq;

import java.util.*;
import java.util.regex.*;

public class MassTool {

    private static final float DOUBLE_ZERO = 1e-6f;
    private static final float PROTON_MASS = 1.00727646688f;

    private final Map<String, Float> mass_table = new HashMap<>();
    private final int missed_cleavage;
    private final int nterm_linkable;

    public MassTool(final int missed_cleavage, Map<String, Float> fix_mod_map, int nterm_linkable) {
        this.missed_cleavage = missed_cleavage;
        this.nterm_linkable = nterm_linkable;
        mass_table.put("G", 57.021464f + fix_mod_map.get("G"));
        mass_table.put("A", 71.037114f + fix_mod_map.get("A"));
        mass_table.put("S", 87.032028f + fix_mod_map.get("S"));
        mass_table.put("P", 97.052764f + fix_mod_map.get("P"));
        mass_table.put("V", 99.068414f + fix_mod_map.get("V"));
        mass_table.put("T", 101.047678f + fix_mod_map.get("I"));
        mass_table.put("C", 103.009184f + fix_mod_map.get("C"));
        mass_table.put("I", 113.084064f + fix_mod_map.get("I"));
        mass_table.put("L", 113.084064f + fix_mod_map.get("L"));
        mass_table.put("N", 114.042927f + fix_mod_map.get("N"));
        mass_table.put("D", 115.026943f + fix_mod_map.get("D"));
        mass_table.put("Q", 128.058578f + fix_mod_map.get("Q"));
        mass_table.put("K", 128.094963f + fix_mod_map.get("K"));
        mass_table.put("E", 129.042593f + fix_mod_map.get("E"));
        mass_table.put("M", 131.040485f + fix_mod_map.get("M"));
        mass_table.put("H", 137.058912f + fix_mod_map.get("H"));
        mass_table.put("F", 147.068414f + fix_mod_map.get("F"));
        mass_table.put("R", 156.101111f + fix_mod_map.get("R"));
        mass_table.put("Y", 163.063329f + fix_mod_map.get("Y"));
        mass_table.put("W", 186.079313f + fix_mod_map.get("W"));
        mass_table.put("X", 118.805717f);
        mass_table.put("Z", fix_mod_map.get("Z"));
        mass_table.put("Hatom", 1.007825032f);
        mass_table.put("Natom", 14.00307401f);
        mass_table.put("Oatom", 15.99491462f);
        mass_table.put("Patom", 30.97376151f);
        mass_table.put("Satom", 31.97207069f);
    }

    public float calResidueMass(String seq) { // Caution: nterm mod is not included!!!!!
        float total_mass = 0;
        int length = seq.length();
        for (int idx = 0; idx < length; ++idx) {
            total_mass += mass_table.get(seq.substring(idx, idx + 1));
        }

        return total_mass;
    }

    public Set<String> buildChainSet(String pro_seq) {
        Map<Integer, List<int[]>> digest_range_map = digestTrypsin(pro_seq);
        Set<String> chain_seq_set = new HashSet<>();

        for (int i = 0; i <= missed_cleavage; ++i) {
            for (int[] digest_range_1 : digest_range_map.get(i)) {
                String sub_string = pro_seq.substring(digest_range_1[0], digest_range_1[1]);
                if (sub_string.substring(0, sub_string.length() - 1).contains("K")) {
                    // If there is a K in middle, this peptide is a chain.
                    chain_seq_set.add(sub_string);
                }
            }
            if (nterm_linkable == 1) {
                if (digest_range_map.get(i).size() > 0) {
                    int[] digest_range = digest_range_map.get(i).get(0);
                    String sub_string = pro_seq.substring(digest_range[0], digest_range[1]);
                    chain_seq_set.add(sub_string);
                }
            }
        }
        return chain_seq_set;
    }

    public float[][] buildChainIonArray(String pep_chain, int min_charge, int max_charge) {
        // [NOTE] The b/y-ions charge 0
        int charge_num = max_charge - min_charge + 1;
        float[][] chain_ion_array = new float[2 * charge_num][pep_chain.length()];
        float b_ion_mass = mass_table.get("Z");
        float y_ion_mass = calResidueMass(pep_chain) + mass_table.get("Z") + 2 * mass_table.get("Hatom") + mass_table.get("Oatom");

        for (int charge = min_charge; charge <= max_charge; ++charge) {
            float b_ion_mass_charge = (b_ion_mass + charge * PROTON_MASS) / charge;
            float y_ion_mass_charge = (y_ion_mass + charge * PROTON_MASS) / charge;
            int charge_idx = charge - min_charge;
            int charge_idx_2 = 2 * charge_idx;
            int charge_idx_2_1 = charge_idx_2 + 1;

            for (int idx = 0; idx < pep_chain.length(); ++idx) {
                // y-ion
                chain_ion_array[charge_idx_2_1][idx] = y_ion_mass_charge;

                String aa = pep_chain.substring(idx, idx + 1);

                // b-ion
                b_ion_mass_charge += mass_table.get(aa) / charge;
                chain_ion_array[charge_idx_2][idx] = b_ion_mass_charge;

                // Calculate next y-ion
                if (idx == 0) {
                    y_ion_mass_charge -= (mass_table.get(aa) / charge + mass_table.get("Z") / charge);
                } else {
                    y_ion_mass_charge -= mass_table.get(aa) / charge;
                }
            }
        }

        return chain_ion_array;
    }

    public Map<String, Float> returnMassTable() {
        return mass_table;
    }

    public float[] buildVector(float[][] ion_matrix, int precursor_charge, int min_charge) {
        int col_num = ion_matrix[0].length;

        int max_row = Math.min(ion_matrix.length / 2, precursor_charge - min_charge) * 2;
        float[] temp = new float[max_row * col_num];
        for (int i = 0; i < max_row; ++i) {
            System.arraycopy(ion_matrix[i], 0, temp, i * col_num, col_num);
        }
        Arrays.sort(temp);

        int left_border = 0;
        int right_border = temp.length;
        boolean has_duplicate = false;
        for (int i = 0; i < temp.length - 1; ++i) {
            if (Math.abs(temp[i] - temp[i+1]) < DOUBLE_ZERO) {
                temp[i] = 0;
                ++left_border;
                has_duplicate = true;
            }
        }

        if (has_duplicate) {
            Arrays.sort(temp);
            return Arrays.copyOfRange(temp, left_border, right_border);
        } else {
            return Arrays.copyOfRange(temp, left_border, right_border);
        }
    }

    public float[][] buildPseudoCLIonArray(float[][] seq_ion, int link_site, int[] common_ion_charge, int[] xlink_ion_charge, float additional_mass) {
        int charge_num = xlink_ion_charge[xlink_ion_charge.length - 1] - common_ion_charge[0] + 1;
        int col_num = seq_ion[0].length;
        float[][] cl_ion_array = new float[2 * charge_num][col_num];

        // Common ion only.
        for (int charge = common_ion_charge[0]; charge < xlink_ion_charge[0]; ++charge) {
            int charge_idx = charge - common_ion_charge[0];
            int charge_idx_2 = 2 * charge_idx;
            int charge_idx_2_1 = charge_idx_2 + 1;
            System.arraycopy(seq_ion[charge_idx_2], 0, cl_ion_array[charge_idx_2], 0, link_site);
            System.arraycopy(seq_ion[charge_idx_2_1], link_site + 1, cl_ion_array[charge_idx_2_1], link_site + 1, col_num - link_site - 1);
        }

        // Common and xlink ion
        for (int charge = xlink_ion_charge[0]; charge <= common_ion_charge[common_ion_charge.length - 1]; ++charge) {
            int charge_idx = charge - common_ion_charge[0]; // Caution!
            int charge_idx_2 = 2 * charge_idx;
            int charge_idx_2_1 = charge_idx_2 + 1;
            System.arraycopy(seq_ion[charge_idx_2], 0, cl_ion_array[charge_idx_2], 0, link_site);
            System.arraycopy(seq_ion[charge_idx_2_1], link_site + 1, cl_ion_array[charge_idx_2_1], link_site + 1, col_num - link_site - 1);
            float addition_mz = additional_mass / charge;
            for (int idx = 0; idx < col_num; ++idx) {
                if (idx < link_site) {
                    cl_ion_array[charge_idx_2_1][idx] = seq_ion[charge_idx_2_1][idx] + addition_mz;
                } else if (idx == link_site) {
                    cl_ion_array[charge_idx_2][idx] = seq_ion[charge_idx_2][idx] + addition_mz;
                    cl_ion_array[charge_idx_2_1][idx] = seq_ion[charge_idx_2_1][idx] + addition_mz;
                } else {
                    cl_ion_array[charge_idx_2][idx] = seq_ion[charge_idx_2][idx] + addition_mz;
                }
            }
        }

        // Xlink ion only
        for (int charge = common_ion_charge[common_ion_charge.length - 1] + 1; charge <= xlink_ion_charge[xlink_ion_charge.length - 1]; ++charge) {
            int charge_idx = charge - common_ion_charge[0]; // Caution!
            int charge_idx_2 = 2 * charge_idx;
            int charge_idx_2_1 = charge_idx_2 + 1;

            // Alpha chain
            float addition_mz = additional_mass / charge;
            for (int idx = 0; idx < col_num; ++idx) {
                if (idx < link_site) {
                    cl_ion_array[charge_idx_2_1][idx] = seq_ion[charge_idx_2_1][idx] + addition_mz;
                } else if (idx == link_site) {
                    cl_ion_array[charge_idx_2][idx] = seq_ion[charge_idx_2][idx] + addition_mz;
                    cl_ion_array[charge_idx_2_1][idx] = seq_ion[charge_idx_2_1][idx] + addition_mz;
                } else {
                    cl_ion_array[charge_idx_2][idx] = seq_ion[charge_idx_2][idx] + addition_mz;
                }
            }
        }

        return cl_ion_array;
    }

    private Map<Integer, List<int[]>> digestTrypsin(String pro_seq) {
        // Cut a protein
        List<Integer> cut_point_list = new LinkedList<>();
        int length = pro_seq.length();
        Pattern rep_1 = Pattern.compile("(K|R)(?=[^P])");
        int idx_start = 0;
        Matcher match_obj = rep_1.matcher(pro_seq);
        cut_point_list.add(0);
        while (idx_start < length) {
            if (match_obj.find()) {
                int cut_point = match_obj.end();
                cut_point_list.add(cut_point);
                idx_start = cut_point;
            } else {
                cut_point_list.add(length);
                break;
            }
        }

        Collections.sort(cut_point_list);

        // Deal with missed cleavage
        Map<Integer, List<int[]>> digest_range_map = new HashMap<>();
        for (int time = 0; time <= missed_cleavage; ++time) {
            List<int[]> temp = new LinkedList<>();
            int left_point;
            int right_point;
            for (int i = 0; i + 1 + time < cut_point_list.size(); ++i) {
                left_point = cut_point_list.get(i);
                right_point = cut_point_list.get(i + 1 + time);
                temp.add(new int[]{left_point, right_point});
            }
            digest_range_map.put(time, temp);
        }

        return digest_range_map;
    }
}
