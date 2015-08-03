# pqpq_python
Rewrite PQPQ (originally developped in MATLAB by Dr. Jenny Shorshed) in python, do clustering at peptide or psm level.

Written by Yafeng Zhu at Karolinska institutet, Sweden.

Email: yafeng.zhu@ki.se

Reference:

Enhanced information output from shotgun proteomics data by protein quantification and peptide quality control (PQPQ). Forshed et al.

Mol Cell Proteomics. 2011 Oct. http://www.ncbi.nlm.nih.gov/pubmed/21734112

Application in splice variant analysis in shotgun proteomics:

SpliceVista, a tool for splice variant identification and visualization in shotgun proteomics data. Zhu et al.

Mol Cell Proteomics. 2014 Jun. http://www.ncbi.nlm.nih.gov/pubmed/24692640

INPUT FORMAT:

A tab delimited PSM or Peptide table. First line is header. 
Peptide sequences and protein accessions locate in first and second column respectively.


EXAMPLE1

python pqpq2.py --in file1 --out file2 --metric correlation --method complete --ratio True --log True --reference 1

'--metric': Available options are correlation (default),cosine, euclidean. 
For all available clustering metric, check on website http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html?highlight=distance.pdist#scipy.spatial.distance.pdist


'--method': options are: complete (default),average, single, weighted.

'--t': distance cutoff used to split clusters. Default value is set 0.4.

'--in': input filename. Tab delimited PSM or Peptide table.

'--out': output filename.

'--log': if set “True”, will do log2 transformation before clustering.

'--ratio’: if set “True”, will calculate ratios using reference sample.

'--reference': an integer to indicate which sample is used as reference to calculate ratios. 1 corresponds to first sample.


EXAMPLE2

If you have a PSM table, but want to do clustering at peptide level. You can use script grouppsm.py to merge psm data into peptide level. Note: In order to detect potential PTMs clustering, you need to differentiate peptide sequences of PTMs to those without PTMs before you group PSMs table to peptide table. It is recommmended to do quantitative clustering at peptide level.

For example:

FYT+79.966DIDELGK is a peptide with phospho mod at Threonine(T), same peptide seq without mods can be written as FYTDIDELGK. Or it can be FYtDIDELGK (modifiied amimo acid is in lower case) vs FYTDIDELGK. As long as you make the sequence different in a way, modified peptides won't be groupped into peptides without modifications.

python grouppsm.py --in file1 --out file2 --method median


'--method': options are: mean or median (default).




