# NIEE
source code for paper "Detecting Tipping Points of Complex Diseases by Network Information Entropy"

## Detecting Tipping Points of Complex Diseases by Network Information Entropy

The progression of complex diseases often involves abrupt and non-linear changes, characterized by sudden shifts that trigger critical transformations. Identifying these critical states or tipping points is crucial for understanding disease progression and developing effective interventions. To address this challenge, we present a model-free method named Network Information Entropy of Edges (NIEE). Leveraging dynamic network biomarkers (DNB), sample-specific network (SSN), and information entropy theories, NIEE detects critical states or tipping points in diverse data types, including bulk, single sample expression data. Demonstrating its effectiveness, we applied NIEE to several real disease datasets, successfully detecting critical points before disease onset. Our findings underscore NIEE's potential to enhance comprehension of complex disease development.

## File description

```NIEE.py``` is the main script of NIEE.

```Preprocessing.ipynb``` is the preprocessing script, aiming to format expression data and generate the ENSP_to_symbol dictionary.

```Result_analysis.ipynb``` is the result analysis script, aiming to screen out the key local network.

## Parameters

The parameters of NIEE.py are as follows:

```-s``` StringDB dataframe input (https://www.string-db.org/cgi/download) and user can use other background network instead of StringDB.

```-sl``` StringDB score limit (default: 450, range: 0-1000).

```-dfexpr``` gene expression dataframe input (row for genes, column for samples).

```-dfanno``` annotation dataframe input.

```-dic``` dictionary for ENSP to symbol (generated from ```Preprocessing.ipynb```).

```-o``` output name.

## Supplementary information

Folders ```KEGG```,```MIC``` and ```Weight``` contain test scripts that support method comparison calculations in Supplementary Note S6-S8.

Folder ```KEGG``` contains the KEGG pathway background network calculation scripts.

Folder ```MIC``` contains the Maximal Information Coefficient calculation scripts.

Folder ```Weight``` contains the weight of protein sequence length calculation scripts.


## Data avalability

The datasets supporting the conclusions of this article can be downloaded from [NCBI GEO](http://www.ncbi.nlm.nih.gov/geo/) and [TCGA](http://cancergenome.nih.gov) repositories.

And the source code is accessible on [GitHub](https://github.com/lllvcs/NIEE).

## Contact

If you have any question, please contact us at [Chengshang LYU](mailto:22395524+lllvcs@users.noreply.github.com) and [Xiaoping Liu](mailto:xpliu@ucas.ac.cn).

Or just feel free to open an issue.
