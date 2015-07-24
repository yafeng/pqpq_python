# pqpq_python
Rewirte pqpq in python, do clustering at peptide level
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

python pqpq2.py --in file1 --out file2 --metric correlation --method --ratio True --log True --reference 1

'--metric': Available options are correlation (default),cosine, euclidean. 

'--method': options are: complete (default),average, single, weighted.

'--t': distance cutoff used to split clusters

'--in': input filename. Tab delimited PSM or Peptide table.

'--out': output filename.

'--log': if set “True”, will do log2 transformation before clustering.

'--ratio’: if set “True”, will calculate ratios using reference sample.

'--reference': an integer to indicate which sample is used as reference to calculate ratios. 1 corresponds to first sample.


EXAMPLE2

If you have a PSM table, but want to do clustering at peptide level. You can use script grouppsm.py to merge psm data into peptide level.

python grouppsm.py --in file1 --out file2 --method median


'--method': options are: mean or median (default).




