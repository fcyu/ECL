package proteomics.Spectrum;

import java.util.*;

public class PreSpectrum {

    private static final double PROTON_MASS = 1.00727646688;
	private static final double DEFAULT_INTENSITY = 1; // DO NOT change. Otherwise, change CalScore accordingly.
    private static final Double DOUBLE_ZERO = 1e-6;

    public static double[][] preProcessSpec(double[][] peaks_array, float precursor_mass, int precursor_charge, double ms2_tolerance, double max_mz) {
        // remove precursor peak from spectrum
        TreeMap<Double, Double> temp = removePrecursorPeak(peaks_array, precursor_mass, precursor_charge, ms2_tolerance);

        // reduce noise
        TreeMap<Double, Double> denoised_pl_map = deNoise(new TreeMap<>(temp.subMap(0d, max_mz)));

        // normalize
        TreeMap<Double, Double> normzlized_pl_map = normalizeSpec(denoised_pl_map);

        return prepareXcorr(normzlized_pl_map);
    }

    private static TreeMap<Double, Double> removePrecursorPeak(double[][] peak_array, float precursor_mass, int precursor_charge, double ms2_tolerance) {
        TreeMap<Double, Double> mz_intensity_map = new TreeMap<>();

        for (int i = 0; i < peak_array[0].length; ++i) {
            for (int charge = precursor_charge; charge > 0; --charge) {
                float mz = (float) (precursor_mass + charge * PROTON_MASS) / charge;
                if ((peak_array[1][i] > DOUBLE_ZERO) && (Math.abs(peak_array[0][i] - mz) > ms2_tolerance)) {
                    mz_intensity_map.put(peak_array[0][i], peak_array[1][i]);
                }
            }
        }

        return mz_intensity_map;
    }

//    private static TreeMap<Double, Double> deNoise(SortedMap<Double, Double> pl_map) {
//        TreeMap<Double, Double> denoised_pl_map = new TreeMap<>();
//        for (double mz : pl_map.keySet()) {
//            double left_mz = Math.max(pl_map.firstKey(), mz - 25);
//            double right_mz = Math.min(pl_map.lastKey(), mz + 25);
//            SortedMap<Double, Double> sub_map = pl_map.subMap(left_mz, right_mz);
//            Double[] intensity_array = sub_map.values().toArray(new Double[sub_map.size()]);
//            if (intensity_array.length > 6) {
//                Arrays.sort(intensity_array);
//                double threshold = intensity_array[intensity_array.length - 6];
//                if (pl_map.get(mz) > threshold) {
//                    denoised_pl_map.put(mz, pl_map.get(mz));
//                }
//            } else {
//                denoised_pl_map.put(mz, pl_map.get(mz));
//            }
//        }
//
//        return denoised_pl_map;
//    }

    private static TreeMap<Double, Double> deNoise(TreeMap<Double, Double> pl_map) {
        TreeMap<Double, Double> denoised_pl_map = new TreeMap<>();
        double window_size = (pl_map.lastKey() - pl_map.firstKey()) / 10 + 1;
        for (int i = 0; i < 10; ++i) {
            double left_mz = pl_map.firstKey() + i * window_size;
            double right_mz = Math.min(left_mz + window_size, pl_map.lastKey());
            NavigableMap<Double, Double> sub_pl_map = pl_map.subMap(left_mz, true, right_mz, true);
            if (sub_pl_map.size() > 4) {
                double noise_intensity = estimateNoiseIntensity(sub_pl_map);
                for (double mz : sub_pl_map.keySet()) {
                    if (sub_pl_map.get(mz) > noise_intensity) {
                        denoised_pl_map.put(mz, sub_pl_map.get(mz));
                    }
                }
            } else {
                for (double mz : sub_pl_map.keySet()) {
                    denoised_pl_map.put(mz, sub_pl_map.get(mz));
                }
            }
        }
        return denoised_pl_map;
    }

    private static double estimateNoiseIntensity(Map<Double, Double> pl) {
        Set<Double> intensity_set = new HashSet<>(pl.values());
        Double[] unique_intensity_vector = intensity_set.toArray(new Double[intensity_set.size()]);
        Arrays.sort(unique_intensity_vector);
        double[] cum = new double[unique_intensity_vector.length];
        for (int i = 0; i < unique_intensity_vector.length; ++i) {
            for (double intensity : pl.values()) {
                if (intensity <= unique_intensity_vector[i]) {
                    ++cum[i];
                }
            }
        }
        double[][] diff = new double[2][unique_intensity_vector.length - 1];
        for (int i = 0; i < unique_intensity_vector.length - 1; ++i) {
            diff[0][i] = cum[i + 1] - cum[i];
            diff[1][i] = unique_intensity_vector[i + 1] - unique_intensity_vector[i];
        }
        double[] diff_2 = new double[unique_intensity_vector.length - 1];
        for (int i = 0; i < unique_intensity_vector.length - 1; ++i) {
            diff_2[i] = diff[0][i] / (diff[1][i] + DOUBLE_ZERO);
        }
        double max_value = 0;
        int max_idx = 0;
        for (int i = 0; i < unique_intensity_vector.length - 1; ++i) {
            if (diff_2[i] > max_value) {
                max_value = diff_2[i];
                max_idx = i;
            }
        }

        return unique_intensity_vector[max_idx]; // TODO: improve
    }

    private static TreeMap<Double, Double> normalizeSpec(TreeMap<Double, Double> pl_map) {
        // sqrt the intensity and find the highest intensity.
        TreeMap<Double, Double> sqrt_pl_map = new TreeMap<>();
        double highest_intensity = 0;
        Set<Double> mz_set = pl_map.keySet();
        for (double mz : mz_set) {
            if (pl_map.get(mz) > DOUBLE_ZERO) {
                double sqrt_intensity = Math.sqrt(pl_map.get(mz));
                if (sqrt_intensity > highest_intensity) {
                    highest_intensity = sqrt_intensity;
                }
                sqrt_pl_map.put(mz, sqrt_intensity);
            }
        }

        // normalize the highest intensity to DEFAULT_INTENSITY
        double factor = DEFAULT_INTENSITY / highest_intensity;
        mz_set = sqrt_pl_map.keySet();
        for (double mz : mz_set) {
            sqrt_pl_map.put(mz, sqrt_pl_map.get(mz) * factor);
        }

        // divide the spectrum into 10 windows and normalize each windows to DEFAULT_INTENSITY
        TreeMap<Double, Double> windowed_pl_map = new TreeMap<>();
        double window_size = (sqrt_pl_map.lastKey() - sqrt_pl_map.firstKey()) / 10 + 1;
        for (int i = 0; i < 10; ++i) {
            // find the max intensity in each window
            double left_mz = sqrt_pl_map.firstKey() + i * window_size;
            double right_mz = Math.min(left_mz + window_size, sqrt_pl_map.lastKey());
            NavigableMap<Double, Double> sub_map = sqrt_pl_map.subMap(left_mz, true, right_mz, true);
            if (!sub_map.isEmpty()) {
                Double[] intensity_array = sub_map.values().toArray(new Double[sub_map.size()]);
                double temp_1 = DEFAULT_INTENSITY / intensity_array[intensity_array.length - 1];
                double temp_2 = 0.05 * intensity_array[intensity_array.length - 1];
                for (double mz : sub_map.keySet()) {
                    if (sub_map.get(mz) > temp_2) {
                        windowed_pl_map.put(mz, sub_map.get(mz) * temp_1);
                    }
                }
            }
        }

        return windowed_pl_map;
    }

    private static double[][] prepareXcorr(TreeMap<Double, Double> mz_intensity_map) {
        double[][] output = new double[2][mz_intensity_map.size()];

        Collection<Double> intensity_list = mz_intensity_map.values();
        double temp = 0;
        for (double intensity : intensity_list) {
            temp += intensity * intensity;
        }
        temp = Math.sqrt(temp);

        Set<Double> mz_set = mz_intensity_map.keySet();
        int i = 0;
        for (double mz : mz_set) {
            output[0][i] = mz;
            output[1][i] = mz_intensity_map.get(mz) / temp;
            ++i;
        }

        return output;
    }
}
