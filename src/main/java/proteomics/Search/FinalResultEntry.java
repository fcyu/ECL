package proteomics.Search;

public class FinalResultEntry implements Comparable<FinalResultEntry> {

    public int spectrum_id = 0;
    public int rank = 0;
    public int charge = 0;
    public float spectrum_precursor_mz = 0;
    public float abs_ppm = 0;
    public double score = 0;
    public double delta_score = 0;
    public String seq_1 = null;
    public int link_site_1 = -1;
    public String mod_1 = null;
    public String pro_id_1 = null;
    public String seq_2 = null;
    public int link_site_2 = -1;
    public String mod_2 = null;
    public String pro_id_2 = null;
    public String cl_type = null;
    public String type = null;
    public int C13_correction = 0;
    public float qvalue = 0;

    public FinalResultEntry(int spectrum_id, int rank, int charge, float spectrum_precursor_mz, float abs_ppm, double score, double delta_score, String seq_1, int link_site_1, String mod_1, String pro_id_1, String seq_2, int link_site_2, String mod_2, String pro_id_2, String cl_type, String type, int C13_correction, float qvalue) {
        this.spectrum_id = spectrum_id;
        this.rank = rank;
        this.charge = charge;
        this.spectrum_precursor_mz= spectrum_precursor_mz;
        this.abs_ppm = abs_ppm;
        this.score = score;
        this.delta_score = delta_score;
        this.seq_1 = seq_1;
        this.link_site_1 = link_site_1;
        this.mod_1 = mod_1;
        this.pro_id_1 = pro_id_1;
        this.seq_2 = seq_2;
        this.link_site_2 = link_site_2;
        this.mod_2 = mod_2;
        this.pro_id_2 = pro_id_2;
        this.cl_type = cl_type;
        this.type = type;
        this.C13_correction = C13_correction;
        this.qvalue = qvalue;
    }

    public int compareTo(FinalResultEntry other) {
        if (this.score > other.score) {
            return 1;
        } else if (this.score < other.score) {
            return -1;
        } else {
            return 0;
        }
    }
}
