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

package proteomics.Search;

public class ResultEntry {

    public String chain_seq_1 = "";
    public String chain_seq_2 = "";
    public float ppm = 0;
    public double score = 0;
    public double second_score = 0;
    public int link_site_1 = 0;
    public int link_site_2 = 0;
    public int C13_correction = 0;

    public ResultEntry(String chain_seq_1, String chain_seq_2, int link_site_1, int link_site_2, float ppm, double score, double second_score, int C13_correction) {
        this.chain_seq_1 = chain_seq_1;
        this.chain_seq_2 = chain_seq_2;
        this.link_site_1 = link_site_1;
        this.link_site_2 = link_site_2;
        this.ppm = ppm;
        this.score = score;
        this.second_score = second_score;
        this.C13_correction = C13_correction;
    }
}
