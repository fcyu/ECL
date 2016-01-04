package proteomics.Index;

import java.util.*;

public class ChainEntry {

    public float chain_mass = 0;
    public String pro_id = null;
    public String chain_type = null;
    public List<Integer> linksite_list = null;
    public float[][] chain_ion_array = null;

    public ChainEntry(float chain_mass, String pro_id, String chain_type, List<Integer> linksite_list, float[][] chain_ion_array) {
        this.chain_mass = chain_mass;
        this.pro_id = pro_id;
        this.chain_type = chain_type;
        this.linksite_list = linksite_list;
        this.chain_ion_array = chain_ion_array;
    }
}