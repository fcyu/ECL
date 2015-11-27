package proteomics.Validation;

import java.util.*;
import proteomics.Search.*;

public class CalFDR {

	private static final float MAX_SCORE = 1;
	private static final float PRECISION = 0.001f;
	private static final float DELTA_SCORE_T = 0.95f;

	private float[] qvalue_array = null;
    private List<FinalResultEntry> results;

	public CalFDR(List<FinalResultEntry> results, boolean unique_peptide) {
        this.results = results;
		final int array_length = 1 + (int) Math.ceil(MAX_SCORE / PRECISION);
		float[] decoy_count_vector = new float[array_length];
		float[] target_count_vector = new float[array_length];
		float[] fuse_count_vector = new float[array_length];
		float[] fdr_array = new float[array_length];
		qvalue_array = new float[array_length];

        if (unique_peptide) {
            Map<String, FinalResultEntry> peptide_result_map = new HashMap<>();
            for (FinalResultEntry result_entry : results) {
                if (result_entry.delta_score > DELTA_SCORE_T) {
                    continue;
                }
                String peptide = result_entry.seq_1 + "-" + result_entry.seq_2;
                if (peptide_result_map.containsKey(peptide)) {
                    if (peptide_result_map.get(peptide).score < result_entry.score) {
                        peptide_result_map.put(peptide, result_entry);
                    }
                } else {
                    peptide_result_map.put(peptide, result_entry);
                }
            }
            for (FinalResultEntry re : peptide_result_map.values()) {
                if (re.type.contentEquals("00")) {
                    int idx = (int) Math.floor(re.score / PRECISION);
                    ++decoy_count_vector[idx];
                } else if (re.type.contentEquals("11")) {
                    int idx = (int) Math.floor(re.score / PRECISION);
                    ++target_count_vector[idx];
                } else {
                    int idx = (int) Math.floor(re.score / PRECISION);
                    ++fuse_count_vector[idx];
                }
            }
        } else {
            for (FinalResultEntry re : results) {
                if (re.delta_score > DELTA_SCORE_T) {
                    continue;
                }
                if (re.type.contentEquals("00")) {
                    int idx = (int) Math.floor(re.score / PRECISION);
                    ++decoy_count_vector[idx];
                } else if (re.type.contentEquals("11")) {
                    int idx = (int) Math.floor(re.score / PRECISION);
                    ++target_count_vector[idx];
                } else {
                    int idx = (int) Math.floor(re.score / PRECISION);
                    ++fuse_count_vector[idx];
                }
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
			if (fuse_count < decoy_count) {
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
        List<FinalResultEntry> final_results = new LinkedList<>();
		for (FinalResultEntry re : results) {
			if (re.type.contentEquals("11")) {
                if (re.delta_score <= DELTA_SCORE_T) {
                    int idx = (int) Math.floor(re.score / PRECISION);
                    float qvalue = qvalue_array[idx];
                    final_results.add(new FinalResultEntry(re.spectrum_id, re.rank, re.charge, re.spectrum_precursor_mz, re.abs_ppm, re.score, re.delta_score, re.seq_1, re.link_site_1, re.mod_1, re.pro_id_1, re.seq_2, re.link_site_2, re.mod_2, re.pro_id_2, re.cl_type, re.type, qvalue));
                }
			}
		}

		return final_results;
	}
}