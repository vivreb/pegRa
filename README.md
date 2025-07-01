# pegRa: Peptide Group Analysis

pegRa is a pipeline to evaluate proteomics data on a peptide level. It is based on the assumption that if a protein changes in abundance, all the peptides from this protein follow the same trend. If a peptide does not follow the same trend, it likely reports on a post-translational modification (PTM), a splicing event, or in case a limited proteolysis (LiP) step was done, a structural change.

During data preprocessing, peptides are filtered and imputed. For each protein, peptides are then assigned into two groups based on fold change. We assume that the larger group reports on the protein abundance change. Each peptide in the smaller group is then tested against the larger group, to see if it is an outlier. Peptides that are significantly different from the large group are considered structural changes.

## Inputs and outputs

The input should be the output from a mass spectrometry search engine. We use Spectronaut, but files from other search engines can in principle be adapted to fit the input file format. An example file can be found on PRIDE Project PXD035183 (https://www.ebi.ac.uk/pride/archive/projects/PXD035183), file "directDIA_Report_rapamycin_vs_DMSO_1.tsv". 

The output is a table listing all the significantly changing proteins, as well as a PDF file showing the plots of the peptide profiles. Optionally, the table listing all peptides can also be exported.

## How to use

First save pegRa.R, example.R, and the example file listed above into the same folder. Run example.R line by line. 

