package proteomics.Validation;

import java.util.*;
import java.util.regex.*;

public class GenerateDecoySeq {

	private String type = null;
	private Pattern fix_pattern = null;

	public GenerateDecoySeq(String type) {
		this.type = type;

		fix_pattern = Pattern.compile("[KR]");
	}

	public String generateSeq(String sequence) {
		String decoy_seq = new String();
		if (type.contentEquals("shuffle")) {
			decoy_seq = shuffleSeq(sequence);
		} else if (type.contentEquals("reverse")) {
			decoy_seq = reverseSeq(sequence);
		}

		return decoy_seq;
	}

	/**
	 * Shuffle a sequnce based on Fisher-Yates shuffle method
	 *
	 * @param sequence
	 * @return
	 */
	private String shuffleSeq(String sequence) {
		String decoy_str = new String();
		Random rnd = new Random();
		Matcher fix_matcher = fix_pattern.matcher(sequence);
		int sequence_length = sequence.length();
		int idx_1 = 0;
		int idx_2 = 0;
		while (idx_1 < sequence_length) {
			String fix_aa = new String();
			if (fix_matcher.find()) {
				idx_2 = fix_matcher.start();
				fix_aa = sequence.substring(idx_2, idx_2 + 1);
			} else {
				idx_2 = sequence_length;
				fix_aa = "";
			}
			String part = sequence.substring(idx_1, idx_2);

			// Fisher-Yates shuffle method
			char[] part_array = part.toCharArray();
			int array_num = part_array.length;
			for (int idx_3 = array_num - 1; idx_3 > 0; --idx_3) {
				int idx_4 = rnd.nextInt(idx_3 + 1);
				char temp = part_array[idx_4];
				part_array[idx_4] = part_array[idx_3];
				part_array[idx_3] = temp;
			}

			// Put shuffled string to the result sequence.
			decoy_str += new String(part_array) + fix_aa;

			idx_1 = idx_2 + 1;
		}

		return decoy_str;
	}

	private String reverseSeq(String sequence) {
		String decoy_str = new String();
		Matcher fix_matcher = fix_pattern.matcher(sequence);
		int sequence_length = sequence.length();
		int idx_1 = 0;
		int idx_2 = 0;
		while (idx_1 < sequence_length) {
			String fix_aa = new String();
			if (fix_matcher.find()) {
				idx_2 = fix_matcher.start();
				fix_aa = sequence.substring(idx_2, idx_2 + 1);
			} else {
				idx_2 = sequence_length;
				fix_aa = "";
			}
			String part = sequence.substring(idx_1, idx_2);

			// Reverse part sequence
			StringBuilder temp = new StringBuilder(part);
			StringBuilder reverse = temp.reverse();
			decoy_str += reverse.toString() + fix_aa;
			idx_1 = idx_2 + 1;
		}

		return decoy_str;
	}
}
