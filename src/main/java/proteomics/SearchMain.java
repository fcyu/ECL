package proteomics;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Parameter.Parameter;
import proteomics.Search.FinalResultEntry;
import proteomics.Search.PrepareSearch;
import proteomics.Search.Search;
import proteomics.Validation.CalFDR;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class SearchMain {

    public static final Logger logger = LoggerFactory.getLogger(SearchMain.class);

    public static void main(String[] args) throws Exception {
        // Process inputs
        if (args.length != 2) {
            help();
        }

        // Set parameters
        String parameter_path = args[0].trim();
        String msxml_path = args[1].trim();

        // Get the parameter map
        Parameter parameter = new Parameter(parameter_path);
        Map<String, String> parameter_map = parameter.returnParameterMap();

        // Prepare search
        logger.info("Indexing database...");
        PrepareSearch ps = new PrepareSearch(parameter_map);
        Map<String, String> pro_annotate_map = ps.returnBuildIndex().getProAnnotateMap();

        // Searching...
        Search search = new Search(ps, parameter_map);
        List<FinalResultEntry> search_results = search.doSearch(msxml_path);

        if (search_results.isEmpty()) {
            logger.warn("There is no useful PSM.");
        } else {
            // save result
            logger.info("Estimating q-value...");
            List<List<FinalResultEntry>> picked_result = pickResult(search_results);
            CalFDR cal_fdr_obj = new CalFDR(picked_result.get(0));
            List<FinalResultEntry> intra_result = cal_fdr_obj.includeStats();
            Collections.sort(intra_result, Collections.<FinalResultEntry>reverseOrder());
            cal_fdr_obj = new CalFDR(picked_result.get(1));
            List<FinalResultEntry> inter_result = cal_fdr_obj.includeStats();
            Collections.sort(inter_result, Collections.<FinalResultEntry>reverseOrder());
            logger.info("Saving results...");
            saveResult(intra_result, inter_result, pro_annotate_map, msxml_path);
        }
        logger.info("Done.");
    }

    private static void saveResult(List<FinalResultEntry> intra_result, List<FinalResultEntry> inter_result, Map<String, String> pro_annotate_map, String id_file_name) throws Exception {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(id_file_name + ".intra.target.csv"))) {
            writer.write("scan_num,spectrum_precursor_mz,charge,score,delta_score,abs_ppm,peptide_1,site_1,mod_1,protein_1,protein_annotation_1,peptide_2,site_2,mod_2,protein_2,protein_annotation_2,q_value\n");
            for (FinalResultEntry re : intra_result) {
                if (re.type.contentEquals("11")) {
                    int link_site_1 = re.link_site_1 + 1;
                    int link_site_2 = re.link_site_2 + 1;
                    String pro_1 = re.pro_id_1;
                    if (re.pro_id_1.contains("&")) {
                        pro_1 = re.pro_id_1.split("&")[0];
                    }
                    String pro_2 = re.pro_id_2;
                    if (re.pro_id_2.contains("&")) {
                        pro_2 = re.pro_id_2.split("&")[0];
                    }
                    String annotate_1 = pro_annotate_map.get(pro_1).replace(",", ";");
                    String annotate_2 = pro_annotate_map.get(pro_2).replace(",", ";");
                    writer.write(re.spectrum_id + "," + re.spectrum_precursor_mz + "," + re.charge + "," + String.format("%.4f", re.score) + "," + String.format("%.2f", re.delta_score) + "," + String.format("%.2f", re.ppm) + "," + re.seq_1 + "," + link_site_1 + "," + re.mod_1 + "," + re.pro_id_1 + ',' + annotate_1 + "," + re.seq_2 + "," + link_site_2 + "," + re.mod_2 + "," + re.pro_id_2 + "," + annotate_2 + "," + String.format("%.4f", re.qvalue) + "\n");
                }
            }
        } catch (IOException ex) {
            logger.error(ex.getMessage());
            System.exit(1);
        }

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(id_file_name + ".intra.decoy.csv"))) {
            writer.write("scan_num,spectrum_precursor_mz,charge,score,delta_score,abs_ppm,peptide_1,site_1,mod_1,protein_1,protein_annotation_1,peptide_2,site_2,mod_2,protein_2,protein_annotation_2\n");
            for (FinalResultEntry re : intra_result) {
                if (re.type.contentEquals("11")) {
                    int link_site_1 = re.link_site_1 + 1;
                    int link_site_2 = re.link_site_2 + 1;
                    writer.write(re.spectrum_id + "," + re.spectrum_precursor_mz + "," + re.charge + "," + String.format("%.4f", re.score) + "," + String.format("%.2f", re.delta_score) + "," + String.format("%.2f", re.ppm) + "," + re.seq_1 + "," + link_site_1 + "," + re.mod_1 + "," + re.pro_id_1 + ",decoy," + re.seq_2 + "," + link_site_2 + "," + re.mod_2 + "," + re.pro_id_2 + ",decoy\n");
                }
            }
        } catch (IOException ex) {
            logger.error(ex.getMessage());
            System.exit(1);
        }

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(id_file_name + ".inter.target.csv"))) {
            writer.write("scan_num,spectrum_precursor_mz,charge,score,delta_score,abs_ppm,peptide_1,site_1,mod_1,protein_1,protein_annotation_1,peptide_2,site_2,mod_2,protein_2,protein_annotation_2,q_value\n");
            for (FinalResultEntry re : inter_result) {
                if (re.type.contentEquals("11")) {
                    int link_site_1 = re.link_site_1 + 1;
                    int link_site_2 = re.link_site_2 + 1;
                    String pro_1 = re.pro_id_1;
                    if (re.pro_id_1.contains("&")) {
                        pro_1 = re.pro_id_1.split("&")[0];
                    }
                    String pro_2 = re.pro_id_2;
                    if (re.pro_id_2.contains("&")) {
                        pro_2 = re.pro_id_2.split("&")[0];
                    }
                    String annotate_1 = pro_annotate_map.get(pro_1).replace(",", ";");
                    String annotate_2 = pro_annotate_map.get(pro_2).replace(",", ";");
                    writer.write(re.spectrum_id + "," + re.spectrum_precursor_mz + "," + re.charge + "," + String.format("%.4f", re.score) + "," + String.format("%.2f", re.delta_score) + "," + String.format("%.2f", re.ppm) + "," + re.seq_1 + "," + link_site_1 + "," + re.mod_1 + "," + re.pro_id_1 + ',' + annotate_1 + "," + re.seq_2 + "," + link_site_2 + "," + re.mod_2 + "," + re.pro_id_2 + "," + annotate_2 + "," + String.format("%.4f", re.qvalue) + "\n");
                }
            }
        } catch (IOException ex) {
            logger.error(ex.getMessage());
            System.exit(1);
        }

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(id_file_name + ".inter.decoy.csv"))) {
            writer.write("scan_num,spectrum_precursor_mz,charge,score,delta_score,abs_ppm,peptide_1,site_1,mod_1,protein_1,protein_annotation_1,peptide_2,site_2,mod_2,protein_2,protein_annotation_2,q_value\n");
            for (FinalResultEntry re : inter_result) {
                if (re.type.contentEquals("11")) {
                    int link_site_1 = re.link_site_1 + 1;
                    int link_site_2 = re.link_site_2 + 1;
                    writer.write(re.spectrum_id + "," + re.spectrum_precursor_mz + "," + re.charge + "," + String.format("%.4f", re.score) + "," + String.format("%.2f", re.delta_score) + "," + String.format("%.2f", re.ppm) + "," + re.seq_1 + "," + link_site_1 + "," + re.mod_1 + "," + re.pro_id_1 + ",decoy," + re.seq_2 + "," + link_site_2 + "," + re.mod_2 + "," + re.pro_id_2 + ",decoy\n");
                }
            }
        } catch (IOException ex) {
            logger.error(ex.getMessage());
            System.exit(1);
        }
    }

    private static List<List<FinalResultEntry>> pickResult(List<FinalResultEntry> search_result) {
        List<List<FinalResultEntry>> picked_result = new LinkedList<>();
        List<FinalResultEntry> inter_protein_result = new LinkedList<>();
        List<FinalResultEntry> intra_protein_result = new LinkedList<>();

        for (FinalResultEntry result_entry : search_result) {
            if (result_entry.cl_type.contentEquals("intra_protein")) {
                intra_protein_result.add(result_entry);
            } else {
                inter_protein_result.add(result_entry);
            }
        }

        picked_result.add(intra_protein_result);
        picked_result.add(inter_protein_result);

        return picked_result;
    }

    private static void help() {
        String help_str = "ECL version 20160117\r\n"
                + "A cross-linked peptides identification tool.\r\n"
                + "Author: Fengchao Yu\r\n"
                + "Email: fyuab@connect.ust.hk\r\n"
                + "ECL usage: java -Xmx25g -jar /path/to/ECL.jar <parameter_file> <data_file>\r\n"
                + "\t<parameter_file>: parameter file. Can be download along with ECL.\r\n"
                + "\t<data_file>: spectra data file (mzXML)\r\n"
                + "\texample: java -Xmx32g -jar ECL.jar parameter.def data.mzxml\r\n";
        System.out.print(help_str);
        System.exit(1);
    }
}
