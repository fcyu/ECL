package proteomics.Parameter;

import java.io.*;
import java.util.*;
import java.util.regex.*;

public class Parameter {
	private String parameter_file = "./src/main/resources/parameter.def";
	private Map<String, String> parameter_map = new HashMap<>();
	
	public Parameter() throws Exception {
		Pattern comment_line_pattern = Pattern.compile("^#.*");
		Pattern comment_pattern = Pattern.compile("([^#]+)#?.*");
		
		BufferedReader parameter_reader = new BufferedReader(new FileReader(parameter_file));
		try {
			String line = "";
			while ((line = parameter_reader.readLine()) != null) {
				line = line.trim();
				Matcher comment_line_matcher = comment_line_pattern.matcher(line);
				if (!comment_line_matcher.matches()) {
					// This is not a comment line
					Matcher line_matcher = comment_pattern.matcher(line);
					String parameter_string;
					if (line_matcher.matches()) {
						parameter_string = line_matcher.group(1).trim();
						String[] parameter_array = parameter_string.split("=");
						String parameter_name = parameter_array[0].trim();
						String parameter_value = parameter_array[1].trim();
						parameter_map.put(parameter_name, parameter_value);
					}
				}
			}
		} catch (IOException ex) {
			System.err.println("IOException: " + ex.getMessage());
			System.exit(1);
		} catch (IllegalStateException ex) {
			System.err.println("IllegalStateException: " + ex.getMessage());
			System.exit(1);
		} finally {
			parameter_reader.close();
		}
	}
	
	
	public Parameter(String parameter_file) throws Exception {
		this.parameter_file = parameter_file;
		
		Pattern comment_line_pattern = Pattern.compile("^#.*");
		Pattern comment_pattern = Pattern.compile("([^#]+)#?.*");
		
		BufferedReader parameter_reader = new BufferedReader(new FileReader(parameter_file));
		try {
			String line = "";
			while ((line = parameter_reader.readLine()) != null) {
				line = line.trim();
				Matcher comment_line_matcher = comment_line_pattern.matcher(line);
				if (!comment_line_matcher.matches()) {
					// This is not a comment line
					Matcher line_matcher = comment_pattern.matcher(line);
					String parameter_string;
					if (line_matcher.matches()) {
						parameter_string = line_matcher.group(1).trim();
						String[] parameter_array = parameter_string.split("=");
						String parameter_name = parameter_array[0].trim();
						String parameter_value = parameter_array[1].trim();
						parameter_map.put(parameter_name, parameter_value);
					}
				}
			}
		} catch (IOException ex) {
			System.err.println("IOException: " + ex.getMessage());
			System.exit(1);
		} catch (IllegalStateException ex) {
			System.err.println("IllegalStateException: " + ex.getMessage());
			System.exit(1);
		} finally {
			parameter_reader.close();
		}
	}
	
	
	public Map<String, String> returnParameterMap() {
		return parameter_map;
	}
}