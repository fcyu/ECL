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

		try (BufferedReader parameter_reader = new BufferedReader(new FileReader(parameter_file))) {
			String line;
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
			parameter_reader.close();
		} catch (IOException | IllegalStateException ex) {
			ex.printStackTrace();
			System.exit(1);
		}
	}
	
	
	public Parameter(String parameter_file) throws Exception {
		this.parameter_file = parameter_file;
		
		Pattern comment_line_pattern = Pattern.compile("^#.*");
		Pattern comment_pattern = Pattern.compile("([^#]+)#?.*");

		try (BufferedReader parameter_reader = new BufferedReader(new FileReader(parameter_file))) {
			String line;
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
			parameter_reader.close();
		} catch (IOException | IllegalStateException ex) {
			ex.printStackTrace();
			System.exit(1);
		}
	}
	
	
	public Map<String, String> returnParameterMap() {
		return parameter_map;
	}
}