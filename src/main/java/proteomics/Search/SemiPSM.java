package proteomics.Search;


public class SemiPSM {
    public int scan_num = 0;
    public float dot_product = 0;
    public int ion_num = 0;
    public int link_site = 0;
    public String chain_seq = "";
    public int C13_correction = 0;

    public SemiPSM(int scan_num, float dot_product, int ion_num, int link_site, String chain_seq, int C13_correction) {
        this.scan_num = scan_num;
        this.dot_product = dot_product;
        this.ion_num = ion_num;
        this.link_site = link_site;
        this.chain_seq = chain_seq;
        this.C13_correction = C13_correction;
    }
}

