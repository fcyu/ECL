package proteomics.Math;


public class CalScore {

	private double[][] alignment_matrix = null;

	public CalScore(double[][] exp_matrix, double[] theo_vector, double xcorr_tolerance) {
		alignment_matrix = new double[2][exp_matrix[0].length + theo_vector.length];

		int start = 0;
		int idx = 0;
		for (double theo_mz : theo_vector) {
			if (exp_matrix[0][0] - theo_mz > xcorr_tolerance) {
				continue;
			}
			
			int j = start;
			while (j < exp_matrix[0].length) {
				double exp_mz = exp_matrix[0][j];
				double temp = exp_mz - theo_mz;
				if ((temp >= -xcorr_tolerance) && (temp <= xcorr_tolerance)) {
					alignment_matrix[0][idx] = exp_mz;
					alignment_matrix[1][idx] = exp_matrix[1][j];
					++idx;
				} else if (temp > xcorr_tolerance) {
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