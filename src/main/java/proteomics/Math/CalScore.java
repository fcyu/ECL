package proteomics.Math;


public class CalScore {

    private float[][] alignment_matrix = null;

    public CalScore(float[][] exp_matrix, float[] theo_vector, float ms2_tolerance) {
        alignment_matrix = new float[2][exp_matrix[0].length + theo_vector.length];

        int start = 0;
        int idx = 0;
        for (float theo_mz : theo_vector) {
            if (exp_matrix[0][0] - theo_mz > ms2_tolerance) {
                continue;
            }

            int j = start;
            while (j < exp_matrix[0].length) {
                float exp_mz = exp_matrix[0][j];
                float temp = exp_mz - theo_mz;
                if ((temp >= -ms2_tolerance) && (temp <= ms2_tolerance)) {
                    alignment_matrix[0][idx] = exp_mz;
                    alignment_matrix[1][idx] = exp_matrix[1][j];
                    ++idx;
                } else if (temp > ms2_tolerance) {
                    start = j - 1;
                    break;
                }
                ++j;
            }
        }
    }

    public float cal_dot_product() {
        float dot_value = 0;
        for (int i = 0; i < alignment_matrix[0].length; ++i) {
            dot_value += alignment_matrix[1][i];
        }

        return dot_value;
    }
}