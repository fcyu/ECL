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

package proteomics.Search;

public class SpectrumEntry {
    public int scan_num = 0;
    public float precursor_intensity = 0;
    public float precursor_mz = 0;
    public float precursor_mass = 0;
    public int precursor_charge = 0;
    public float[][] mz_intensity_array = null;

    public SpectrumEntry(int scan_num, float precursor_intensity, float precursor_mz, float precursor_mass, int precursor_charge, float[][] mz_intensity_array) {
        this.scan_num = scan_num;
        this.precursor_intensity = precursor_intensity;
        this.precursor_mz = precursor_mz;
        this.precursor_mass = precursor_mass;
        this.precursor_charge = precursor_charge;
        this.mz_intensity_array = mz_intensity_array;
    }
}