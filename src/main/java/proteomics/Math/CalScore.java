/*
 * Copyright 2015-2017 The Hong Kong University of Science and Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package proteomics.Math;


public class CalScore {

    private double[][] alignment_matrix = null;

    public CalScore(float[][] exp_matrix, float[] theo_vector, float ms2_tolerance) {
        alignment_matrix = new double[2][exp_matrix[0].length + theo_vector.length];

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

    public double cal_dot_product() {
        double dot_value = 0;
        for (int i = 0; i < alignment_matrix[0].length; ++i) {
            dot_value += alignment_matrix[1][i];
        }

        return dot_value;
    }
}