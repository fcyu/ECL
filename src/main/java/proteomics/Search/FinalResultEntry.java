package proteomics.Search;

public class FinalResultEntry {

	public int spectrum_id = 0;
	public String cl_index = null;
	public int rank = 0;
	public int charge = 0;
	public float spectrum_precursor_mz = 0;
	public float abs_ppm = 0;
	public double xcorr = 0;
	public double delta_xcorr = 0;
	public String sequence_1 = null;
	public String mod_1 = null;
	public String protein_id_1 = null;
	public String sequence_2 = null;
	public String mod_2 = null;
	public String protein_id_2 = null;
	public String cl_type = null;
	public String type = null;
	public double qvalue = 0;

	public FinalResultEntry(int spectrum_id, String cl_index, int rank, int charge, float spectrum_precursor_mz, float abs_ppm, double xcorr, double delta_xcorr, String sequence_1, String mod_1, String protein_id_1, String sequence_2, String mod_2, String protein_id_2, String cl_type, String type, double qvalue) {
		this.spectrum_id = spectrum_id;
		this.cl_index = cl_index;
		this.rank = rank;
		this.charge = charge;
		this.spectrum_precursor_mz= spectrum_precursor_mz;
		this.abs_ppm = abs_ppm;
		this.xcorr = xcorr;
		this.delta_xcorr = delta_xcorr;
		this.sequence_1 = sequence_1;
		this.mod_1 = mod_1;
		this.protein_id_1 = protein_id_1;
		this.sequence_2 = sequence_2;
		this.mod_2 = mod_2;
		this.protein_id_2 = protein_id_2;
		this.cl_type = cl_type;
		this.type = type;
		this.qvalue = qvalue;
	}
}
