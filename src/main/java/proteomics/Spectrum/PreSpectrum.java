package proteomics.Spectrum;

import java.util.*;

public class PreSpectrum {

    private static final float PROTON_MASS = 1.00727646688f;
	private static final float DEFAULT_INTENSITY = 1; // DO NOT change. Otherwise, change CalScore accordingly.
    private static final float FLOAT_ZERO = 1e-6f;

    public static float[][] preProcessSpec(Map<Double, Double> peaks_map, float precursor_mass, int precursor_charge, float ms2_tolerance, float max_mz) {
        // remove precursor peak from spectrum
        TreeMap<Float, Float> temp = removePrecursorPeak(peaks_map, precursor_mass, precursor_charge, ms2_tolerance);

        // reduce noise
        TreeMap<Float, Float> denoised_pl_map = deNoise(new TreeMap<>(temp.subMap(0f, max_mz)));

        // normalize
        TreeMap<Float, Float> normzlized_pl_map = normalizeSpec(denoised_pl_map);

        return prepareXcorr(normzlized_pl_map);
    }

    private static TreeMap<Float, Float> removePrecursorPeak(Map<Double, Double> peak_map, float precursor_mass, int precursor_charge, float ms2_tolerance) {
        TreeMap<Float, Float> mz_intensity_map = new TreeMap<>();

        for (double mz : peak_map.keySet()) {
            for (int charge = precursor_charge; charge > 0; --charge) {
                float temp = (precursor_mass + charge * PROTON_MASS) / charge;
                if ((peak_map.get(mz) > FLOAT_ZERO) && (Math.abs(peak_map.get(mz) - temp) > ms2_tolerance)) {
                    mz_intensity_map.put((float) mz, peak_map.get(mz).floatValue());
                }
            }
        }

        return mz_intensity_map;
    }

    private static TreeMap<Float, Float> deNoise(TreeMap<Float, Float> pl_map) {
        TreeMap<Float, Float> denoised_pl_map = new TreeMap<>();
        float window_size = (pl_map.lastKey() - pl_map.firstKey()) / 10 + 1;
        for (int i = 0; i < 10; ++i) {
            float left_mz = pl_map.firstKey() + i * window_size;
            float right_mz = Math.min(left_mz + window_size, pl_map.lastKey());
            NavigableMap<Float, Float> sub_pl_map = pl_map.subMap(left_mz, true, right_mz, true);
            if (sub_pl_map.size() > 4) {
                float noise_intensity = estimateNoiseIntensity(sub_pl_map);
                for (float mz : sub_pl_map.keySet()) {
                    if (sub_pl_map.get(mz) > noise_intensity) {
                        denoised_pl_map.put(mz, sub_pl_map.get(mz));
                    }
                }
            } else {
                for (float mz : sub_pl_map.keySet()) {
                    denoised_pl_map.put(mz, sub_pl_map.get(mz));
                }
            }
        }
        return denoised_pl_map;
    }

    private static float estimateNoiseIntensity(Map<Float, Float> pl) {
        Set<Float> intensity_set = new HashSet<>(pl.values());
        Float[] unique_intensity_vector = intensity_set.toArray(new Float[intensity_set.size()]);
        Arrays.sort(unique_intensity_vector);
        float[] cum = new float[unique_intensity_vector.length];
        for (int i = 0; i < unique_intensity_vector.length; ++i) {
            for (float intensity : pl.values()) {
                if (intensity <= unique_intensity_vector[i]) {
                    ++cum[i];
                }
            }
        }
        float[][] diff = new float[2][unique_intensity_vector.length - 1];
        for (int i = 0; i < unique_intensity_vector.length - 1; ++i) {
            diff[0][i] = cum[i + 1] - cum[i];
            diff[1][i] = unique_intensity_vector[i + 1] - unique_intensity_vector[i];
        }
        float[] diff_2 = new float[unique_intensity_vector.length - 1];
        for (int i = 0; i < unique_intensity_vector.length - 1; ++i) {
            diff_2[i] = diff[0][i] / (diff[1][i] + FLOAT_ZERO);
        }
        float max_value = 0;
        int max_idx = 0;
        for (int i = 0; i < unique_intensity_vector.length - 1; ++i) {
            if (diff_2[i] > max_value) {
                max_value = diff_2[i];
                max_idx = i;
            }
        }

        return unique_intensity_vector[max_idx]; // TODO: improve
    }

    private static TreeMap<Float, Float> normalizeSpec(TreeMap<Float, Float> pl_map) {
        // sqrt the intensity and find the highest intensity.
        TreeMap<Float, Float> sqrt_pl_map = new TreeMap<>();
        float highest_intensity = 0;
        for (float mz : pl_map.keySet()) {
            if (pl_map.get(mz) > FLOAT_ZERO) {
                float sqrt_intensity = (float) Math.sqrt(pl_map.get(mz));
                if (sqrt_intensity > highest_intensity) {
                    highest_intensity = sqrt_intensity;
                }
                sqrt_pl_map.put(mz, sqrt_intensity);
            }
        }

        // normalize the highest intensity to DEFAULT_INTENSITY
        float factor = DEFAULT_INTENSITY / highest_intensity;
        for (float mz : sqrt_pl_map.keySet()) {
            sqrt_pl_map.put(mz, sqrt_pl_map.get(mz) * factor);
        }

        // divide the spectrum into 10 windows and normalize each windows to DEFAULT_INTENSITY
        TreeMap<Float, Float> windowed_pl_map = new TreeMap<>();
        float window_size = (sqrt_pl_map.lastKey() - sqrt_pl_map.firstKey()) / 10 + 1;
        for (int i = 0; i < 10; ++i) {
            // find the max intensity in each window
            float left_mz = sqrt_pl_map.firstKey() + i * window_size;
            float right_mz = Math.min(left_mz + window_size, sqrt_pl_map.lastKey());
            NavigableMap<Float, Float> sub_map = sqrt_pl_map.subMap(left_mz, true, right_mz, true);
            if (!sub_map.isEmpty()) {
                Float[] intensity_array = sub_map.values().toArray(new Float[sub_map.size()]);
                float temp_1 = DEFAULT_INTENSITY / intensity_array[intensity_array.length - 1];
                float temp_2 = (float) 0.05 * intensity_array[intensity_array.length - 1];
                for (float mz : sub_map.keySet()) {
                    if (sub_map.get(mz) > temp_2) {
                        windowed_pl_map.put(mz, sub_map.get(mz) * temp_1);
                    }
                }
            }
        }

        return windowed_pl_map;
    }

    private static float[][] prepareXcorr(TreeMap<Float, Float> mz_intensity_map) {
        float[][] output = new float[2][mz_intensity_map.size()];

        Collection<Float> intensity_list = mz_intensity_map.values();
        float temp = 0;
        for (float intensity : intensity_list) {
            temp += intensity * intensity;
        }
        temp = (float) Math.sqrt(temp);

        int i = 0;
        for (float mz : mz_intensity_map.keySet()) {
            output[0][i] = mz;
            output[1][i] = mz_intensity_map.get(mz) / temp;
            ++i;
        }

        return output;
    }
}
