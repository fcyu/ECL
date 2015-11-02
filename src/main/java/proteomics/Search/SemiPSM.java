package proteomics.Search;


public class SemiPSM {
    public int scan_num = 0;
    public double dot_product = 0;
    public int ion_num = 0;
    public int link_site = 0;
    public String chain_index = "";

    public SemiPSM(int scan_num, double dot_product, int ion_num, int link_site, String chain_index) {
        this.scan_num = scan_num;
        this.dot_product = dot_product;
        this.ion_num = ion_num;
        this.link_site = link_site;
        this.chain_index = chain_index;
    }
}

