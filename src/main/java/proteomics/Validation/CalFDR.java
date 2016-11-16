package proteomics.Validation;

import java.util.*;
import proteomics.Search.*;

public class CalFDR {

    private static final int BIN_NUM = 1000;
    private static final float DELTA_SCORE_T = 0.95f;

    private float[] qvalue_array = null;
    private List<FinalResultEntry> results;
    private final double bin_size;

    public CalFDR(List<FinalResultEntry> results) {
        this.results = results;

        // find the max_score
        double max_score = 0;
        for (FinalResultEntry entry : results) {
            if (entry.score > max_score) {
                max_score = entry.score;
            }
        }

        bin_size = max_score / BIN_NUM;

        final int array_length = 1 + BIN_NUM;
        int[] decoy_count_vector = new int[array_length];
        int[] target_count_vector = new int[array_length];
        int[] fuse_count_vector = new int[array_length];
        float[] fdr_array = new float[array_length];
        qvalue_array = new float[array_length];

        for (FinalResultEntry re : results) {
            if (re.delta_score > DELTA_SCORE_T) {
                continue;
            }
            if (re.type.contentEquals("00")) {
                int idx = (int) Math.floor(re.score / bin_size);
                ++decoy_count_vector[idx];
            } else if (re.type.contentEquals("11")) {
                int idx = (int) Math.floor(re.score / bin_size);
                ++target_count_vector[idx];
            } else {
                int idx = (int) Math.floor(re.score / bin_size);
                ++fuse_count_vector[idx];
            }
        }

        // Calculate FDR
        for (int idx_1 = 0; idx_1 < array_length; ++idx_1) {
            int decoy_count = 0;
            int fuse_count = 0;
            int target_count = 0;
            for (int idx_2 = idx_1; idx_2 < array_length; ++idx_2) {
                decoy_count += decoy_count_vector[idx_2];
                target_count += target_count_vector[idx_2];
                fuse_count += fuse_count_vector[idx_2];
            }

            float fdr;
            if ((target_count == 0) && ((decoy_count > 0) || (fuse_count > 0))) {
                fdr = 1;
            } else if ((target_count == 0) && (decoy_count == 0) && (fuse_count == 0)) {
                fdr = 0;
            } else if (fuse_count < decoy_count) {
                fdr = (float) decoy_count / (float) target_count;
            } else {
                fdr = (float) (fuse_count - decoy_count) / (float) target_count;
            }

            fdr = Math.min(fdr, 1); // Adjust those fdrs that are larger than 1
            fdr_array[idx_1] = fdr;
        }

        // Convert FDR to qvalue
        float last_q_value = fdr_array[0];
        qvalue_array[0] = last_q_value;

        for (int idx_1 = 1; idx_1 < array_length; ++idx_1) {
            float q_value = fdr_array[idx_1];
            if (q_value >= last_q_value) {
                qvalue_array[idx_1] = last_q_value;
            } else {
                qvalue_array[idx_1] = q_value;
                last_q_value = q_value;
            }
        }
    }

    public List<FinalResultEntry> includeStats() {
        for (FinalResultEntry re : results) {
            if (re.type.contentEquals("11") && (re.delta_score <= DELTA_SCORE_T)) {
                int idx = (int) Math.floor(re.score / bin_size);
                re.qvalue = qvalue_array[idx];
            }
        }

        return results;
    }
}