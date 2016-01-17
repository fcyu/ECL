package proteomics.Search;

public class ResultEntry {

    public String chain_seq_1 = "";
    public String chain_seq_2 = "";
    public float ppm = 0;
    public double score = 0;
    public double second_score = 0;
    public int link_site_1 = 0;
    public int link_site_2 = 0;
    public int C13_correction = 0;

    public ResultEntry(String chain_seq_1, String chain_seq_2, int link_site_1, int link_site_2, float ppm, double score, double second_score, int C13_correction) {
        this.chain_seq_1 = chain_seq_1;
        this.chain_seq_2 = chain_seq_2;
        this.link_site_1 = link_site_1;
        this.link_site_2 = link_site_2;
        this.ppm = ppm;
        this.score = score;
        this.second_score = second_score;
        this.C13_correction = C13_correction;
    }
}
