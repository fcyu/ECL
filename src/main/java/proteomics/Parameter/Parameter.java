/*
 * Copyright 2015-2017 The Hong Kong University of Science and Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package proteomics.Parameter;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.*;
import java.util.regex.*;

public class Parameter {

    private static final Logger logger = LoggerFactory.getLogger(Parameter.class);

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
            logger.error(ex.getMessage());
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
            logger.error(ex.getMessage());
            System.exit(1);
        }
    }


    public Map<String, String> returnParameterMap() {
        return parameter_map;
    }
}