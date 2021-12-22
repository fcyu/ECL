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

package theoSeq;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.*;
import java.util.regex.*;

public class DbTool {

    private static final Logger logger = LoggerFactory.getLogger(DbTool.class);

    private Map<String, String> pro_seq_map = new HashMap<>();
    private Map<String, String> pro_annotate_map = new HashMap<>();

    public DbTool(String db_name) throws Exception {
        String id = "";
        String annotate;
        String seq = "";

        boolean new_pro = true;

        Pattern header_pattern = Pattern.compile(">([^\\s]*)(.*)");

        try (BufferedReader db_reader = new BufferedReader(new FileReader(db_name))) {
            String line;
            while ((line = db_reader.readLine()) != null) {
                line = line.trim();
                Matcher head_matcher = header_pattern.matcher(line);
                if (head_matcher.matches()) {
                    // This line is a header
                    if (!new_pro) {
                        // This isn't the first protein
                        pro_seq_map.put(id, seq);
                    }
                    id = head_matcher.group(1);
                    annotate = head_matcher.group(2);
                    pro_annotate_map.put(id, annotate);
                    new_pro = true;
                } else if (!line.isEmpty()) {
                    // This line is a body
                    if (new_pro) {
                        seq = line;
                        new_pro = false;
                    } else {
                        seq += line;
                    }
                }
            }
            // Last protein
            pro_seq_map.put(id, seq);
        } catch (IOException | PatternSyntaxException ex) {
            logger.error(ex.getMessage());
            System.exit(1);
        }
    }

    public Map<String, String> returnSeqMap() {
        return pro_seq_map;
    }

    public Map<String, String> returnAnnotateMap() {
        return pro_annotate_map;
    }
}
