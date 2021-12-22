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


public class SemiPSM {
    public int scan_num = 0;
    public double dot_product = 0;
    public int ion_num = 0;
    public int link_site = 0;
    public String chain_seq = "";
    public int C13_correction = 0;

    public SemiPSM(int scan_num, double dot_product, int ion_num, int link_site, String chain_seq, int C13_correction) {
        this.scan_num = scan_num;
        this.dot_product = dot_product;
        this.ion_num = ion_num;
        this.link_site = link_site;
        this.chain_seq = chain_seq;
        this.C13_correction = C13_correction;
    }
}

