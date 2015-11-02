package proteomics.Search;

import org.systemsbiology.jrap.stax.MSXMLParser;
import org.systemsbiology.jrap.stax.Scan;
import org.systemsbiology.jrap.stax.ScanHeader;
import proteomics.Index.ChainEntry;
import proteomics.LogEntry;
import proteomics.Math.CalScore;
import proteomics.Spectrum.PreSpectrum;
import theoSeq.MassTool;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class Search {

    private static final double PROTON_MASS = 1.00727646688;

    private LogEntry log_entry = null;
    private int[] ms1_charge = null;
    private int[] common_ion_charge = null;
    private int[] xlink_ion_charge = null;
    private float ms1_tolerance = 0;
    private int ms1_tolerance_unit = 0;
    private double ms2_tolerance = 0;
    private double linker_mass = 0;
    private float min_precursor_mass = 0;
    private float max_precursor_mass = 0;
    private int min_peak_num = 0;
    private Map<String, String> index_chain_map = null;
    private Map<String, ChainEntry> chain_entry_map = null;
    private Map<String, Double> fix_mod_map = null;
    private MassTool mass_tool = null;
    private String output_str = null;
    private Map<Integer, SpectrumEntry> num_spectrum_map = new HashMap<>();
    private float single_dot_product_t = 1e-6f;
    private TreeMap<Integer, List<SpectrumEntry>> mass1000_spectrum_map = new TreeMap<>();
    private TreeMap<Integer, Set<String>> mass1000_chain_map = new TreeMap<>();
    private int linker_mass1000 = 0;

    /////////////////////////////////////////public methods////////////////////////////////////////////////////////////
    public Search(PrepareSearch ps, LogEntry log_entry) throws Exception {
        this.log_entry = log_entry;

        Map<String, String> parameter_map = ps.returnParameterMap();
        index_chain_map = ps.returnIndexChainMap();
        chain_entry_map = ps.returnChainEntryMap();
        fix_mod_map = ps.returnFixModMap();
        mass_tool = ps.returnMassTool();
        ms1_charge = string2Array(parameter_map.get("ms1_charge"));
        common_ion_charge = string2Array(parameter_map.get("common_ion_charge"));
        xlink_ion_charge = string2Array(parameter_map.get("xlink_ion_charge"));
        ms1_tolerance_unit = Integer.valueOf(parameter_map.get("ms1_tolerance_unit"));
        ms1_tolerance = Float.valueOf(parameter_map.get("ms1_tolerance"));
        ms2_tolerance = Double.valueOf(parameter_map.get("ms2_tolerance"));
        String cl_aa = parameter_map.get("cl_aa");
        linker_mass = Double.valueOf(parameter_map.get("cl_mass")) - 2 * Double.valueOf(parameter_map.get(cl_aa)); // TODO: different link sites have different mod. The compensation should be different.
        linker_mass1000 = round1000((float) linker_mass);
        min_precursor_mass = Float.valueOf(parameter_map.get("min_precursor_mass"));
        max_precursor_mass = Float.valueOf(parameter_map.get("max_precursor_mass"));
        min_peak_num = Integer.valueOf(parameter_map.get("min_peak_num"));
    }

    public List<FinalResultEntry> doSearch(String msxml_path) throws Exception {
        // Read mzxml
        MSXMLParser msxml_parser = null;
        try {
            File f = new File(msxml_path);
            if ((!f.exists() || (f.isDirectory()))) {
                FileNotFoundException file_error = new FileNotFoundException("mzXML file not found.");
                throw file_error;
            }
            msxml_parser = new MSXMLParser(msxml_path);
        } catch (FileNotFoundException ex) {
            System.err.println(ex.getMessage());
            System.exit(1);
        }

        // Get current time
        double time_spectra_start = System.nanoTime();
        System.out.println("Reading and processing spectra...");
        for (int j = 1; j <= msxml_parser.getMaxScanNumber(); ++j) {
            Scan scan = msxml_parser.rap(j);
            ScanHeader header = scan.getHeader();

            if (header.getMsLevel() == 1) {
                continue;
            }

            if (header.getPrecursorCharge() == -1) {
                continue;
            }

            int precursor_charge = header.getPrecursorCharge();
            if ((precursor_charge < ms1_charge[0]) || (precursor_charge > ms1_charge[ms1_charge.length - 1])) {
                continue;
            }

            float precursor_mz = header.getPrecursorMz();
            float precursor_mass = precursor_mz * precursor_charge - precursor_charge * (float) PROTON_MASS;
            if ((precursor_mass > max_precursor_mass + 1) || (precursor_mass < min_precursor_mass - 1)) {
                continue;
            }

            double[][] raw_mz_intensity_array = scan.getMassIntensityList();
            if (raw_mz_intensity_array[0].length < min_peak_num) {
                continue;
            }

            double[][] mz_intensity_array = PreSpectrum.preProcessSpec(raw_mz_intensity_array, precursor_mass, precursor_charge, ms2_tolerance, precursor_mass);
            if (mz_intensity_array[0].length <= min_peak_num) {
                continue;
            }

            SpectrumEntry spectrum_entry = new SpectrumEntry(j, header.getPrecursorIntensity(), precursor_mz, precursor_mass, precursor_charge, mz_intensity_array);

            if (mass1000_spectrum_map.containsKey(round1000(precursor_mass))) {
                List<SpectrumEntry> spectrum_list = mass1000_spectrum_map.get(round1000(precursor_mass));
                spectrum_list.add(spectrum_entry);
                mass1000_spectrum_map.put(round1000(precursor_mass), spectrum_list);
            } else {
                List<SpectrumEntry> spectrum_list = new LinkedList<>();
                spectrum_list.add(spectrum_entry);
                mass1000_spectrum_map.put(round1000(precursor_mass), spectrum_list);
            }

            num_spectrum_map.put(j, spectrum_entry);
        }

        double time_spectra_end = System.nanoTime();
        double spectra_duration = (time_spectra_end - time_spectra_start) * 1e-9;
        System.out.println("Duration: " + (int) spectra_duration + " seconds");
        output_str += "Duration: " + (int) spectra_duration + " seconds" + "\r\n";

        // Generate a mass chain map.
        Set<String> chain_set_1 = chain_entry_map.keySet();
        for (String chain : chain_set_1) {
            ChainEntry chain_entry = chain_entry_map.get(chain);
            int mass1000 = round1000(chain_entry.chain_mass);
            if (mass1000_chain_map.containsKey(mass1000)) {
                Set<String> temp = mass1000_chain_map.get(mass1000);
                temp.add(chain);
                mass1000_chain_map.put(mass1000, temp);
            } else {
                Set<String> temp = new HashSet<>();
                temp.add(chain);
                mass1000_chain_map.put(mass1000, temp);
            }
        }

        // Get current time
        Map<Integer, ResultEntry> result_map = new HashMap<>();
        double time_search_start = System.nanoTime();
        System.out.println("Searching cross-linked peptides...");
        output_str += "Searching cross-linked peptides..." + "\r\n";

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
        int last_progress = 0;
        float step = 0;
        float total = (float) mass1000_1_max_idx;
        for (int i = 0; i <= mass1000_1_max_idx; ++i) {
            ++step;
            int progress = (int) Math.floor(step * 10 / total);
            if (progress != last_progress) {
                System.out.print(progress + "0% ");
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
            chain_set_1 = mass1000_chain_map.get(mass1000_1);
            for (String chain_1 : chain_set_1) {
                ChainEntry chain_entry_1 = this.chain_entry_map.get(chain_1);
                List<Integer> link_site_list_1 = chain_entry_1.linksite_list;
                double[][] seq_ion_1 = chain_entry_1.chain_ion_array;
                linearSearch(sub_spectra_map, mass1000_1, chain_entry_1.chain_index, seq_ion_1, link_site_list_1, exp_mass1000_2_semiPSM_map);
            }

            if (exp_mass1000_2_semiPSM_map.size() == 0) {
                continue;
            }

            Set<Integer> exp_mass1000_2_set = exp_mass1000_2_semiPSM_map.keySet();
            for (int exp_mass1000_2 : exp_mass1000_2_set) {
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
                Set<Integer> mass1000_2_set = mass1000_2_map.keySet();
                for (int mass1000_2 : mass1000_2_set) {
                    Set<String> chain_set_2 = mass1000_chain_map.get(mass1000_2);
                    for (String chain_2 : chain_set_2) {
                        ChainEntry chain_entry_2 = chain_entry_map.get(chain_2);
                        double[][] seq_ion_2 = chain_entry_2.chain_ion_array;
                        List<Integer> link_site_list_2 = chain_entry_2.linksite_list;
                        for (int link_site_2 : link_site_list_2) {
                            double[][] pseudo_ion_array = mass_tool.buildPseudoCLIonArray(seq_ion_2, link_site_2, common_ion_charge, xlink_ion_charge, back1000(mass1000_1) + (float) linker_mass);
                            Map<Integer, double[]> charge_theo_mz_map = new HashMap<>();
                            double[] theo_mz_2;
                            for (SemiPSM semi_psm : semi_psm_list) {
                                int scan_num = semi_psm.scan_num;
                                SpectrumEntry spectrum_entry = num_spectrum_map.get(scan_num);
                                int precursor_charge = spectrum_entry.precursor_charge;
                                if (charge_theo_mz_map.containsKey(precursor_charge)) {
                                    theo_mz_2 = charge_theo_mz_map.get(precursor_charge);
                                } else {
                                    theo_mz_2 = mass_tool.buildVector(pseudo_ion_array, precursor_charge, common_ion_charge[0]);
                                    charge_theo_mz_map.put(precursor_charge, theo_mz_2);
                                }
                                CalScore cal_score_2 = new CalScore(spectrum_entry.mz_intensity_array, theo_mz_2, ms2_tolerance);
                                double dot_product_2 = cal_score_2.cal_dot_product();
                                if (dot_product_2 <= single_dot_product_t) {
                                    continue;
                                }

                                // Calculate final score
                                double xcorr = (semi_psm.dot_product + dot_product_2) / Math.sqrt(semi_psm.ion_num + theo_mz_2.length);

                                // record result
                                if (result_map.containsKey(scan_num)) {
                                    ResultEntry last_result = result_map.get(scan_num);
                                    if (xcorr > last_result.xcorr) {
                                        float total_mass = back1000(mass1000_1 + mass1000_2) + (float) linker_mass;
                                        float abs_ppm = (float) (Math.abs(spectrum_entry.precursor_mass - total_mass) * 1e6 / total_mass);
                                        String cl_index = semi_psm.chain_index + "-" + chain_entry_2.chain_index + "-" + semi_psm.link_site + "-" + link_site_2;
                                        ResultEntry result_entry = new ResultEntry(abs_ppm, xcorr, last_result.xcorr, cl_index);
                                        result_map.put(scan_num, result_entry);
                                    } else if (xcorr > last_result.scond_xcorr) {
                                        ResultEntry result_entry = new ResultEntry(last_result.abs_ppm, last_result.xcorr, xcorr, last_result.cl_index);
                                        result_map.put(scan_num, result_entry);
                                    }
                                } else {
                                    float total_mass = back1000(mass1000_1 + mass1000_2) + (float) linker_mass;
                                    float abs_ppm = (float) (Math.abs(spectrum_entry.precursor_mass - total_mass) * 1e6 / total_mass);
                                    String cl_index = semi_psm.chain_index + "-" + chain_entry_2.chain_index + "-" + semi_psm.link_site + "-" + link_site_2;
                                    ResultEntry result_entry = new ResultEntry(abs_ppm, xcorr, 0, cl_index);
                                    result_map.put(scan_num, result_entry);
                                }
                            }
                        }
                    }
                }
            }
        }

        // Build a fix modification only map for further use.
        Map<String, Double> fix_mod_only = new HashMap<>();
        Set<String> aa_set = fix_mod_map.keySet();
        for (String aa : aa_set) {
            double mod_mass = fix_mod_map.get(aa);
            if (mod_mass != 0.0) {
                fix_mod_only.put(aa, mod_mass);
            }
        }

        List<FinalResultEntry> search_result = new LinkedList<>();
        Set<Integer> spectrum_num_set = result_map.keySet();
        for (int spectrum_num : spectrum_num_set) {
            int rank = 1;
            ResultEntry result_entry = result_map.get(spectrum_num);
            String cl_index = result_entry.cl_index;
            String[] temp = cl_index.split("-");
            String chain_index_1 = temp[0] + "-" + temp[1];
            String var_mod_1 = temp[1];
            String chain_index_2 = temp[2] + "-" + temp[3];
            String var_mod_2 = temp[3];
            String alpha_chain = index_chain_map.get(chain_index_1);
            String beta_chain = index_chain_map.get(chain_index_2);
            ChainEntry chain_entry_1 = chain_entry_map.get(alpha_chain);
            ChainEntry chain_entry_2 = chain_entry_map.get(beta_chain);
            String type = chain_entry_1.chain_type + chain_entry_2.chain_type; // 11 means target, 00, 10, 01 are decoy

            // Add fix modification annotation string.
            Set<String> fix_mod_aa = fix_mod_only.keySet();
            String fix_mod_1 = "";
            String fix_mod_2 = "";
            for (String aa : fix_mod_aa) {
                if (alpha_chain.contains(aa)) {
                    fix_mod_1 += fix_mod_only.get(aa) + "@" + aa + ";";
                }
                if (beta_chain.contains(aa)) {
                    fix_mod_2 += fix_mod_only.get(aa) + "@" + aa + ";";
                }
            }
            if (fix_mod_aa.contains("Z")) {
                fix_mod_1 += fix_mod_only.get("Z") + "@Z;";
                fix_mod_2 += fix_mod_only.get("Z") + "@Z;";
            }

            String mod_1 = fix_mod_1 + var_mod_1;
            String mod_2 = fix_mod_2 + var_mod_2;

            ScanHeader scan_header = msxml_parser.rap(spectrum_num).getHeader();

            int precursor_charge = scan_header.getPrecursorCharge();

            String pro_1 = chain_entry_1.pro_id;
            if (chain_entry_1.chain_type.contentEquals("0")) {
                pro_1 = chain_entry_1.pro_id.substring(6);
            }

            String pro_2 = chain_entry_2.pro_id;
            if (chain_entry_2.chain_type.contentEquals("0")) {
                pro_2 = chain_entry_2.pro_id.substring(6);
            }

            String cl_type = "intra_protein";
            if (!pro_1.contentEquals(pro_2)) {
                cl_type = "inter_protein";
            }

            double delta_xcorr = result_entry.scond_xcorr / result_entry.xcorr;

            FinalResultEntry re = new FinalResultEntry(spectrum_num, cl_index, rank, precursor_charge, scan_header.getPrecursorMz(), result_entry.abs_ppm, result_entry.xcorr, delta_xcorr, alpha_chain, mod_1, chain_entry_1.pro_id, beta_chain, mod_2, chain_entry_2.pro_id, cl_type, type, -1);
            search_result.add(re);
        }

        Calendar cal_search_end = Calendar.getInstance();
        double time_search_end = System.nanoTime();
        spectra_duration = (time_search_end - time_search_start) * 1e-9;

        System.out.println();
        System.out.println("Duration: " + (int) spectra_duration + " seconds");
        output_str += "Duration: " + (int) spectra_duration + " seconds" + "\r\n";

        log_entry.output_str = output_str;

        return search_result;
    }

    ///////////////////////////////////////private methods//////////////////////////////////////////////////////////////
    private void linearSearch(NavigableMap<Integer, List<SpectrumEntry>> sub_spectra_map, int mass1000_1, String chain_index, double[][] seq_ion, List<Integer> link_site_list_1, TreeMap<Integer, List<SemiPSM>> exp_mass1000_2_semiPSM_map) {
        Set<Integer> exp_mass1000_set = sub_spectra_map.keySet();
        for (int exp_mass1000 : exp_mass1000_set) {
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
                double[][] pseudo_ion_array = mass_tool.buildPseudoCLIonArray(seq_ion, link_site_1, common_ion_charge, xlink_ion_charge, back1000(exp_mass1000 - mass1000_1));
                List<SpectrumEntry> spectrum_list = sub_spectra_map.get(exp_mass1000);
                Map<Integer, double[]> charge_mz_map = new HashMap<>();
                for (SpectrumEntry spectrum_entry : spectrum_list) {
                    int precursor_charge = spectrum_entry.precursor_charge;
                    double[] theo_mz;
                    if (charge_mz_map.containsKey(precursor_charge)) {
                        theo_mz = charge_mz_map.get(precursor_charge);
                    } else {
                        theo_mz = mass_tool.buildVector(pseudo_ion_array, precursor_charge, common_ion_charge[0]);
                        charge_mz_map.put(precursor_charge, theo_mz);
                    }

                    // Calculate dot produce
                    CalScore cal_score = new CalScore(spectrum_entry.mz_intensity_array, theo_mz, ms2_tolerance);
                    double dot_product = cal_score.cal_dot_product();

                    // Record result
                    if (dot_product > single_dot_product_t) {
                        int exp_mass1000_2 = exp_mass1000 - mass1000_1 - linker_mass1000;
                        if (exp_mass1000_2_semiPSM_map.containsKey(exp_mass1000_2)) {
                            List<SemiPSM> temp = exp_mass1000_2_semiPSM_map.get(exp_mass1000_2);
                            temp.add(new SemiPSM(spectrum_entry.scan_num, dot_product, theo_mz.length, link_site_1, chain_index));
                            exp_mass1000_2_semiPSM_map.put(exp_mass1000_2, temp);
                        } else {
                            List<SemiPSM> temp = new LinkedList<>();
                            temp.add(new SemiPSM(spectrum_entry.scan_num, dot_product, theo_mz.length, link_site_1, chain_index));
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