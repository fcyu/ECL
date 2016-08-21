# ECL
A fast and exhaustive cross-linked peptides identification tool.

## How to use it?
Requirements: Java 1.7 or later

Usage:
```
java -Xmx32g -jar /path/to/ECL.jar <parameter_file> <data_file>
```
1. <parameter_file>: parameter file. Can be downloaded along with ECL.
2. <data_file>: spectra data file (mzXML).

example: java -Xmx25g -jar ECL.jar parameter.def data.mzxml

## Cite
Yu, Fengchao, Ning Li, and Weichuan Yu. "ECL: an exhaustive search tool for the identification of cross-linked peptides using whole database." BMC Bioinformatics 17.1 (2016): 1.

[bibtex]
```bibtex
@article{yu2016ecl,
  title={ECL: an exhaustive search tool for the identification of cross-linked peptides using whole database},
  author={Yu, Fengchao and Li, Ning and Yu, Weichuan},
  journal={BMC Bioinformatics},
  volume={17},
  number={1},
  pages={1},
  year={2016},
  publisher={BioMed Central}
}
```
