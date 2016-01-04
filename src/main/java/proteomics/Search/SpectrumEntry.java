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