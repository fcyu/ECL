package proteomics;

import proteomics.Parameter.Parameter;
import proteomics.Search.FinalResultEntry;
import proteomics.Search.PrepareSearch;
import proteomics.Search.Search;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.*;

public class SearchMain {

    public static void main(String[] args) throws Exception {
        // Get current time
        SimpleDateFormat sdf = new SimpleDateFormat("yyyy.MM.dd 'at' HH:mm:ss Z");
        Calendar cal_start = Calendar.getInstance();
        Date date_start = cal_start.getTime();
        float time_start = System.nanoTime();

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

        LogEntry log_entry = new LogEntry("");

        // Prepare search
        System.out.println("Indexing database...");
        log_entry.output_str += "Indexing database...";
        PrepareSearch ps = new PrepareSearch(parameter_map);

        // Searching...
        Search search = new Search(ps, log_entry, parameter_map);
        List<FinalResultEntry> search_results = search.doSearch(msxml_path);

        if (search_results.isEmpty()) {
            System.out.println("There is no PSM.");
        } else {
            // save result
            List<List<FinalResultEntry>> picked_result = pickResult(search_results);
            System.out.println("Saving results...");
            saveResult2CSV(picked_result.get(0), picked_result.get(1), msxml_path);
        }

        // Get end time
        Calendar cal_end = Calendar.getInstance();
        Date date_end = cal_end.getTime();
        float time_end = System.nanoTime();

        double duration = (time_end - time_start) * 1e-9;

        try (BufferedWriter target_writer = new BufferedWriter(new FileWriter(msxml_path + ".log.txt"))) {
            target_writer.write(msxml_path + " finished.\r\n"
                    + "Started on: " + sdf.format(date_start) + "\r\n"
                    + "Ended on: " + sdf.format(date_end) + "\r\n"
                    + "Duration: " + (int) duration + " second" + "\r\n"
                    + "\r\n"
                    + "Log: \r\n"
                    + log_entry.output_str + "\r\n");
        } catch (IOException ex) {
            System.err.println(ex.getMessage());
        }

        System.out.println("Done.");
    }

    //////////////////////////////////////////private methods///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    private static void saveResult2CSV(List<FinalResultEntry> search_results_intra, List<FinalResultEntry> search_results_inter, String id_file_name) throws Exception {
        try (BufferedWriter intra_writer = new BufferedWriter(new FileWriter(id_file_name + ".intra.xls"))) {
            intra_writer.write("id\tlabel\tscannr\txcorr\tdeltaXcorr\tabsppm\tpeptide\tproteinId1\n");
            for (FinalResultEntry re : search_results_intra) {
                int link_site_1 = re.link_site_1 + 1;
                int link_site_2 = re.link_site_2 + 1;
                if (re.type.contentEquals("11")) {
                    intra_writer.write(re.spectrum_id + "." + re.charge + "\t" + "1" + "\t" + re.spectrum_id + "\t" + re.xcorr + "\t" + re.delta_xcorr + "\t" + re.abs_ppm + "\t" + "-." + re.seq_1 + "-" + link_site_1 + "-" + re.seq_2 + "-" + link_site_2 + ".-" + "\t" + re.pro_id_1 + "-" + re.pro_id_2 + "\n");
                } else {
                    intra_writer.write(re.spectrum_id + "." + re.charge + "\t" + "-1" + "\t" + re.spectrum_id + "\t" + re.xcorr + "\t" + re.delta_xcorr + "\t" + re.abs_ppm + "\t" + "-." + re.seq_1 + "-" + link_site_1 + "-" + re.seq_2 + "-" + link_site_2 + ".-" + "\t" + re.pro_id_1 + "-" + re.pro_id_2 + "\n");
                }
            }
        } catch (IOException ex) {
            System.err.println("IOException: " + ex.getMessage());
            System.exit(1);
        }

        try (BufferedWriter inter_writer = new BufferedWriter(new FileWriter(id_file_name + ".inter.xls"))) {
            inter_writer.write("id\tlabel\tscannr\txcorr\tdeltaXcorr\tabsppm\tpeptide\tproteinId1\n");
            for (FinalResultEntry re : search_results_inter) {
                int link_site_1 = re.link_site_1 + 1;
                int link_site_2 = re.link_site_2 + 1;
                if (re.type.contentEquals("11")) {
                    inter_writer.write(re.spectrum_id + "." + re.charge + "\t" + "1" + "\t" + re.spectrum_id + "\t" + re.xcorr + "\t" + re.delta_xcorr + "\t" + re.abs_ppm + "\t" + "-." + re.seq_1 + "-" + link_site_1 + "-" + re.seq_2 + "-" + link_site_2 + ".-" + "\t" + re.pro_id_1 + "-" + re.pro_id_2 + "\n");
                } else {
                    inter_writer.write(re.spectrum_id + "." + re.charge + "\t" + "-1" + "\t" + re.spectrum_id + "\t" + re.xcorr + "\t" + re.delta_xcorr + "\t" + re.abs_ppm + "\t" + "-." + re.seq_1 + "-" + link_site_1 + "-" + re.seq_2 + "-" + link_site_2 + ".-" + "\t" + re.pro_id_1 + "-" + re.pro_id_2 + "\n");
                }
            }
        } catch (IOException ex) {
            System.err.println("IOException: " + ex.getMessage());
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
        String help_str = "ECL version 20151111\r\n"
                + "A cross-linked peptides identification tool.\r\n"
                + "Author: Fengchao Yu\r\n"
                + "Email: fyuab@connect.ust.hk\r\n"
                + "ECL usage: java -Xmx32g -jar /path/to/ECL.jar <parameter_file> <data_file>\r\n"
                + "\t<parameter_file>: parameter file. Can be download along with ECL.\r\n"
                + "\t<data_file>: spectra data file (mzXML)\r\n"
                + "\texample: java -Xmx32g -jar ECL.jar parameter.def data.mzxml";
        System.out.print(help_str);
        System.exit(1);
    }
}
