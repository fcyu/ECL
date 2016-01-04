package proteomics.Search;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import uk.ac.ebi.pride.tools.jmzreader.model.*;
import uk.ac.ebi.pride.tools.mzxml_parser.*;
import proteomics.Index.BuildIndex;
import proteomics.Index.ChainEntry;
import proteomics.Math.CalScore;
import proteomics.Spectrum.PreSpectrum;
import theoSeq.MassTool;

import java.io.File;
import java.util.*;

public class Search {

    private static final float PROTON_MASS = 1.00727646688f;
    private static final Logger logger = LoggerFactory.getLogger(Search.class);

    private int[] ms1_charge = null;
    private int[] common_ion_charge = null;
    private int[] xlink_ion_charge = null;
    private float ms1_tolerance = 0;
    private int ms1_tolerance_unit = 0;
    private float ms2_tolerance = 0;
    private float linker_mass = 0;
    private float min_precursor_mass = 0;
    private float max_precursor_mass = 0;
    private int min_peak_num = 0;
    private Map<String, ChainEntry> chain_entry_map = null;
    private Map<String, Float> fix_mod_map = null;
    private MassTool mass_tool_obj = null;
    private Map<Integer, SpectrumEntry> num_spectrum_map = new HashMap<>();
    private float single_dot_product_t = 1e-6f;
    private TreeMap<Integer, List<SpectrumEntry>> mass1000_spectrum_map = new TreeMap<>();
    private TreeMap<Integer, Set<String>> mass1000_chain_map = new TreeMap<>();
    private int linker_mass1000 = 0;

    /////////////////////////////////////////public methods////////////////////////////////////////////////////////////
    public Search(PrepareSearch ps, Map<String, String> parameter_map) throws Exception {
        chain_entry_map = ps.returnChainEntryMap();
        BuildIndex build_index_obj = ps.returnBuildIndex();
        fix_mod_map = build_index_obj.returnFixModMap();
        mass_tool_obj = build_index_obj.returnMassTool();
        ms1_charge = string2Array(parameter_map.get("ms1_charge"));
        common_ion_charge = string2Array(parameter_map.get("common_ion_charge"));
        xlink_ion_charge = string2Array(parameter_map.get("xlink_ion_charge"));
        ms1_tolerance_unit = Integer.valueOf(parameter_map.get("ms1_tolerance_unit"));
        ms1_tolerance = Float.valueOf(parameter_map.get("ms1_tolerance")) * 1.1f; // consider rounding error.
        ms2_tolerance = Float.valueOf(parameter_map.get("ms2_tolerance"));
        linker_mass = Float.valueOf(parameter_map.get("cl_mass")) - 2 * Float.valueOf(parameter_map.get("K"));
        linker_mass1000 = round1000(linker_mass);
        min_precursor_mass = Float.valueOf(parameter_map.get("min_precursor_mass"));
        max_precursor_mass = Float.valueOf(parameter_map.get("max_precursor_mass"));
        min_peak_num = Integer.valueOf(parameter_map.get("min_peak_num"));
    }

    public List<FinalResultEntry> doSearch(String msxml_path) throws Exception {
        // Read mzxml
		logger.info("Reading and processing spectra...");
        MzXMLFile spectra_parser = null;
        try {
            File spectra_file = new File(msxml_path);
            if ((!spectra_file.exists() || (spectra_file.isDirectory()))) {
                logger.error("The spectra file not found.");
            }
            spectra_parser = new MzXMLFile(spectra_file);
        } catch (MzXMLParsingException ex) {
            logger.error(ex.getMessage());
            System.exit(1);
        }

        // Get current time
        Iterator<Spectrum> spectrum_iterator = spectra_parser.getSpectrumIterator();
        while (spectrum_iterator.hasNext()) {
            Spectrum spectrum = spectrum_iterator.next();

            if (spectrum.getMsLevel() != 2) {
                continue;
            }

            if (spectrum.getPrecursorCharge() == null) {
                continue;
            }
            int precursor_charge = spectrum.getPrecursorCharge();

            if ((precursor_charge < ms1_charge[0]) || (precursor_charge > ms1_charge[ms1_charge.length - 1])) {
                continue;
            }

            float precursor_mz = spectrum.getPrecursorMZ().floatValue();
            float precursor_mass = precursor_mz * precursor_charge - precursor_charge * PROTON_MASS;
            if ((precursor_mass > max_precursor_mass + 1) || (precursor_mass < min_precursor_mass - 1)) {
                continue;
            }

            Map<Double, Double> raw_mz_intensity_map = spectrum.getPeakList();
            if (raw_mz_intensity_map.size() < min_peak_num) {
                continue;
            }

            float[][] mz_intensity_array = PreSpectrum.preProcessSpec(raw_mz_intensity_map, precursor_mass, precursor_charge, ms2_tolerance, precursor_mass);
            if (mz_intensity_array[0].length <= min_peak_num) {
                continue;
            }

            int scan_num = Integer.valueOf(spectrum.getId());
            SpectrumEntry spectrum_entry = new SpectrumEntry(scan_num, spectrum.getPrecursorIntensity().floatValue(), precursor_mz, precursor_mass, precursor_charge, mz_intensity_array);

            if (mass1000_spectrum_map.containsKey(round1000(precursor_mass))) {
                mass1000_spectrum_map.get(round1000(precursor_mass)).add(spectrum_entry);
            } else {
                List<SpectrumEntry> spectrum_list = new LinkedList<>();
                spectrum_list.add(spectrum_entry);
                mass1000_spectrum_map.put(round1000(precursor_mass), spectrum_list);
            }

            num_spectrum_map.put(scan_num, spectrum_entry);
        }

        // Generate a mass chain map.
        for (String chain : chain_entry_map.keySet()) {
            ChainEntry chain_entry = chain_entry_map.get(chain);
            int mass1000 = round1000(chain_entry.chain_mass);
            if (mass1000_chain_map.containsKey(mass1000)) {
                mass1000_chain_map.get(mass1000).add(chain);
            } else {
                Set<String> temp = new HashSet<>();
                temp.add(chain);
                mass1000_chain_map.put(mass1000, temp);
            }
        }

        Map<Integer, ResultEntry> result_map = new HashMap<>();
        logger.info("Searching cross-linked peptides...");

        Integer[] mass1000_array = mass1000_chain_map.keySet().toArray(new Integer[mass1000_chain_map.size()]);
        int max_spectrum_mass1000 = mass1000_spectrum_map.lastKey();
        int mass1000_1_max = (max_spectrum_mass1000 - linker_mass1000) / 2;
        int mass1000_1_max_idx = 0;
        for (int i = 0; i < mass1000_array.length; ++i) {
            if (mass1000_array[i] >= mass1000_1_max) {
                mass1000_1_max_idx = i;
                break;
            }
        }

        // There is no spectrum meet the requirement.
        if (mass1000_1_max_idx == 0) {
            return new LinkedList<>();
        }

        int last_progress = -1;
        for (int i = 0; i <= mass1000_1_max_idx; ++i) {
            int progress = i * 20 / mass1000_1_max_idx;
            if (progress != last_progress) {
                logger.info("Searching progress: {}%", progress * 5);
                last_progress = progress;
            }

            int mass1000_1 = mass1000_array[i];

            int start_mass = 2 * mass1000_1 + linker_mass1000;
            if (ms1_tolerance_unit == 0) {
                start_mass -= round1000(ms1_tolerance);
            } else if (ms1_tolerance_unit == 1) {
                start_mass -= (int) (start_mass * ms1_tolerance / 1e6f);
            }
            if (start_mass >= max_spectrum_mass1000) {
                break;
            }

            NavigableMap<Integer, List<SpectrumEntry>> sub_spectra_map = mass1000_spectrum_map.subMap(start_mass, true, max_spectrum_mass1000, true);
            if (sub_spectra_map.size() == 0) {
                continue;
            }

            TreeMap<Integer, List<SemiPSM>> exp_mass1000_2_semiPSM_map = new TreeMap<>();
            for (String chain_seq_1 : mass1000_chain_map.get(mass1000_1)) {
                ChainEntry chain_entry_1 = this.chain_entry_map.get(chain_seq_1);
                List<Integer> link_site_list_1 = chain_entry_1.linksite_list;
                float[][] seq_ion_1 = chain_entry_1.chain_ion_array;
                linearSearch(sub_spectra_map, chain_seq_1, mass1000_1, seq_ion_1, link_site_list_1, exp_mass1000_2_semiPSM_map);
            }

            if (exp_mass1000_2_semiPSM_map.size() == 0) {
                continue;
            }

            for (int exp_mass1000_2 : exp_mass1000_2_semiPSM_map.keySet()) {
                int exp_mass1000 = mass1000_1 + exp_mass1000_2 + linker_mass1000;
                int mass1000_2_left = 0;
                int mass1000_2_right = 0;
                if (ms1_tolerance_unit == 0) {
                    mass1000_2_left = exp_mass1000_2 - round1000(ms1_tolerance);
                    mass1000_2_right = exp_mass1000_2 + round1000(ms1_tolerance);
                } else if (ms1_tolerance_unit == 1) {
                    mass1000_2_left = (int) (exp_mass1000 / (1 + ms1_tolerance / 1e6)) - mass1000_1 - linker_mass1000;
                    mass1000_2_right = (int) (exp_mass1000 / (1 - ms1_tolerance / 1e6)) - mass1000_1 - linker_mass1000;
                }
                NavigableMap<Integer, Set<String>> mass1000_2_map = mass1000_chain_map.subMap(mass1000_2_left, true, mass1000_2_right, true);
                List<SemiPSM> semi_psm_list = exp_mass1000_2_semiPSM_map.get(exp_mass1000_2);
                for (int mass1000_2 : mass1000_2_map.keySet()) {
                    Set<String> chain_set_2 = mass1000_chain_map.get(mass1000_2);
                    for (String chain_seq_2 : chain_set_2) {
                        ChainEntry chain_entry_2 = chain_entry_map.get(chain_seq_2);
                        float[][] seq_ion_2 = chain_entry_2.chain_ion_array;
                        for (int link_site_2 : chain_entry_2.linksite_list) {
                            float[][] pseudo_ion_array = mass_tool_obj.buildPseudoCLIonArray(seq_ion_2, link_site_2, common_ion_charge, xlink_ion_charge, back1000(mass1000_1) + linker_mass);
                            Map<Integer, float[]> charge_theo_mz_map = new HashMap<>();
                            float[] theo_mz_2;
                            for (SemiPSM semi_psm : semi_psm_list) {
                                int scan_num = semi_psm.scan_num;
                                SpectrumEntry spectrum_entry = num_spectrum_map.get(scan_num);
                                int precursor_charge = spectrum_entry.precursor_charge;
                                if (charge_theo_mz_map.containsKey(precursor_charge)) {
                                    theo_mz_2 = charge_theo_mz_map.get(precursor_charge);
                                } else {
                                    theo_mz_2 = mass_tool_obj.buildVector(pseudo_ion_array, precursor_charge, common_ion_charge[0]);
                                    charge_theo_mz_map.put(precursor_charge, theo_mz_2);
                                }
                                CalScore cal_score_2 = new CalScore(spectrum_entry.mz_intensity_array, theo_mz_2, ms2_tolerance);
                                float dot_product_2 = cal_score_2.cal_dot_product();
                                if (dot_product_2 <= single_dot_product_t) {
                                    continue;
                                }

                                // Calculate final score
                                double score = (semi_psm.dot_product + dot_product_2) / Math.sqrt(semi_psm.ion_num + theo_mz_2.length);

                                // record result
                                if (result_map.containsKey(scan_num)) {
                                    ResultEntry last_result = result_map.get(scan_num);
                                    if (score > last_result.score) {
                                        float total_mass = back1000(mass1000_1 + mass1000_2) + linker_mass;
                                        float abs_ppm = (float) (Math.abs(spectrum_entry.precursor_mass - total_mass) * 1e6 / total_mass);
                                        ResultEntry result_entry = new ResultEntry(semi_psm.chain_seq, chain_seq_2, semi_psm.link_site, link_site_2, abs_ppm, score, last_result.score);
                                        result_map.put(scan_num, result_entry);
                                    } else if (score > last_result.second_score) {
                                        result_map.get(scan_num).second_score = score;
                                    }
                                } else {
                                    float total_mass = back1000(mass1000_1 + mass1000_2) + linker_mass;
                                    float abs_ppm = (float) (Math.abs(spectrum_entry.precursor_mass - total_mass) * 1e6 / total_mass);
                                    ResultEntry result_entry = new ResultEntry(semi_psm.chain_seq, chain_seq_2, semi_psm.link_site, link_site_2, abs_ppm, score, -1);
                                    result_map.put(scan_num, result_entry);
                                }
                            }
                        }
                    }
                }
            }
        }

        // Build a fix modification only map for further use.
        Map<String, Float> fix_mod_only = new HashMap<>();
        for (String aa : fix_mod_map.keySet()) {
            float mod_mass = fix_mod_map.get(aa);
            if (mod_mass != 0.0) {
                fix_mod_only.put(aa, mod_mass);
            }
        }

        List<FinalResultEntry> search_result = new LinkedList<>();
        for (int spectrum_num : result_map.keySet()) {
            int rank = 1;
            ResultEntry result_entry = result_map.get(spectrum_num);
            String chain_seq_1 = result_entry.chain_seq_1;
            String chain_seq_2 = result_entry.chain_seq_2;
            ChainEntry chain_entry_1 = chain_entry_map.get(chain_seq_1);
            ChainEntry chain_entry_2 = chain_entry_map.get(chain_seq_2);
            String type = chain_entry_1.chain_type + chain_entry_2.chain_type; // 11 means target, 00, 10, 01 are decoy

            // Add fix modification annotation string.
            Set<String> fix_mod_aa = fix_mod_only.keySet();
            String fix_mod_1 = "";
            String fix_mod_2 = "";
            for (String aa : fix_mod_aa) {
                if (chain_seq_1.contains(aa)) {
                    fix_mod_1 += fix_mod_only.get(aa) + "@" + aa + ";";
                }
                if (chain_seq_2.contains(aa)) {
                    fix_mod_2 += fix_mod_only.get(aa) + "@" + aa + ";";
                }
            }
            if (fix_mod_aa.contains("Z")) {
                fix_mod_1 += fix_mod_only.get("Z") + "@Z;";
                fix_mod_2 += fix_mod_only.get("Z") + "@Z;";
            }

            String mod_1 = fix_mod_1;
            String mod_2 = fix_mod_2;

            SpectrumEntry spectrum_entry = num_spectrum_map.get(spectrum_num);

            int precursor_charge = spectrum_entry.precursor_charge;

            // analyze cross-linking type: inter or intra
            String[] pro_1;
            if (chain_entry_1.chain_type.contentEquals("0")) {
                pro_1 = new String[]{chain_entry_1.pro_id.substring(6)};
            } else {
                pro_1 = chain_entry_1.pro_id.split("&");
            }

            String[] pro_2;
            if (chain_entry_2.chain_type.contentEquals("0")) {
                pro_2 = new String[]{chain_entry_2.pro_id.substring(6)};
            } else {
                pro_2 = chain_entry_2.pro_id.split("&");
            }

            String cl_type = "intra_protein";
            boolean keep = false;
            for (String temp_1 : pro_1) {
                for (String temp_2 : pro_2) {
                    if (temp_1.contentEquals(temp_2)) {
                        keep = true;
                        break;
                    }
                }
            }
            if ((!keep) || chain_seq_1.contentEquals(chain_seq_2) || chain_seq_1.contains(chain_seq_2) || chain_seq_2.contains(chain_seq_1)) {
                cl_type = "inter_protein";
            }

            double delta_score = 1;
            if (Math.abs(result_entry.second_score + 1) > 1e-6) {
                delta_score = result_entry.second_score / result_entry.score;
            }

            FinalResultEntry re = new FinalResultEntry(spectrum_num, rank, precursor_charge, spectrum_entry.precursor_mz, result_entry.abs_ppm, result_entry.score, delta_score, chain_seq_1, result_entry.link_site_1, mod_1, chain_entry_1.pro_id, chain_seq_2, result_entry.link_site_2, mod_2, chain_entry_2.pro_id, cl_type, type, -1);
            search_result.add(re);
        }

        return search_result;
    }

    ///////////////////////////////////////private methods//////////////////////////////////////////////////////////////
    private void linearSearch(NavigableMap<Integer, List<SpectrumEntry>> sub_spectra_map, String chain_seq, int mass1000_1, float[][] seq_ion, List<Integer> link_site_list_1, TreeMap<Integer, List<SemiPSM>> exp_mass1000_2_semiPSM_map) {
        for (int exp_mass1000 : sub_spectra_map.keySet()) {
            // Get the MS1 mass range.
            int mass1000_2_left = 0;
            int mass1000_2_right = 0;
            if (ms1_tolerance_unit == 0) {
                mass1000_2_left = exp_mass1000 - round1000(ms1_tolerance) - mass1000_1 - linker_mass1000;
                mass1000_2_right = exp_mass1000 + round1000(ms1_tolerance) - mass1000_1 - linker_mass1000;
            } else if (ms1_tolerance_unit == 1) {
                mass1000_2_left = (int) (exp_mass1000 / (1 + ms1_tolerance / 1e6f)) - mass1000_1 - linker_mass1000;
                mass1000_2_right = (int) (exp_mass1000 / (1 - ms1_tolerance / 1e6f)) - mass1000_1 - linker_mass1000;
            }

            NavigableMap<Integer, Set<String>> sub_map = mass1000_chain_map.subMap(mass1000_2_left, true, mass1000_2_right, true);
            if (sub_map.size() == 0) {
                continue;
            }

            for (int link_site_1 : link_site_list_1) {
                float[][] pseudo_ion_array = mass_tool_obj.buildPseudoCLIonArray(seq_ion, link_site_1, common_ion_charge, xlink_ion_charge, back1000(exp_mass1000 - mass1000_1));
                List<SpectrumEntry> spectrum_list = sub_spectra_map.get(exp_mass1000);
                Map<Integer, float[]> charge_mz_map = new HashMap<>();
                for (SpectrumEntry spectrum_entry : spectrum_list) {
                    int precursor_charge = spectrum_entry.precursor_charge;
                    float[] theo_mz;
                    if (charge_mz_map.containsKey(precursor_charge)) {
                        theo_mz = charge_mz_map.get(precursor_charge);
                    } else {
                        theo_mz = mass_tool_obj.buildVector(pseudo_ion_array, precursor_charge, common_ion_charge[0]);
                        charge_mz_map.put(precursor_charge, theo_mz);
                    }

                    // Calculate dot produce
                    CalScore cal_score = new CalScore(spectrum_entry.mz_intensity_array, theo_mz, ms2_tolerance);
                    float dot_product = cal_score.cal_dot_product();

                    // Record result
                    if (dot_product > single_dot_product_t) {
                        int exp_mass1000_2 = exp_mass1000 - mass1000_1 - linker_mass1000;
                        if (exp_mass1000_2_semiPSM_map.containsKey(exp_mass1000_2)) {
                            exp_mass1000_2_semiPSM_map.get(exp_mass1000_2).add(new SemiPSM(spectrum_entry.scan_num, dot_product, theo_mz.length, link_site_1, chain_seq));
                        } else {
                            List<SemiPSM> temp = new LinkedList<>();
                            temp.add(new SemiPSM(spectrum_entry.scan_num, dot_product, theo_mz.length, link_site_1, chain_seq));
                            exp_mass1000_2_semiPSM_map.put(exp_mass1000_2, temp);
                        }
                    }
                }
            }
        }
    }

    private int[] string2Array(String str) {
        String[] temp = str.split(",");
        int[] output = new int[temp.length];
        for (int i = 0; i < output.length; ++i) {
            output[i] = Integer.valueOf(temp[i]);
        }

        return output;
    }

    private int round1000(float a) {
        return (int) (a * 1000);
    }

    private float back1000(int a) {
        return (float) a / (float) 1000;
    }
}