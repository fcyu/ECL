package proteomics.Search;

public class ResultEntry {

	public float abs_ppm = 0;
    public double xcorr = 0;
	public String cl_index = null;
    public double scond_xcorr = 0;

    public ResultEntry(float abs_ppm, double xcorr, double second_xcorr, String cl_index) {
		this.abs_ppm = abs_ppm;
        this.xcorr = xcorr;
		this.cl_index = cl_index;
        this.scond_xcorr = second_xcorr;
    }
}
