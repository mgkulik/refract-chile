# REFRACT Chile MGKulik
Repository to deposit relevant code generated during my trip to Santiago through the REFRACT program between Jan/April 2023.

## Main function


The following files are available:

**get_peptide_score.py:** Main function to calculate the score for the provided peptide with 6 standard amino acids.
In position 0 only C, F, W and Y are accepted and they are not accepted in the other 5 positions.
At the moment accepts just 1 peptide as an entry and can be run through the terminal with the following code:

```
python get_peptide_score.py -i <peptide>
```

*Example:* python get_peptide_score.py -i FEKDAD

## Additional functions


The following scripts must be executed IN ORDER to re-generate all outputs:<br>
Source files are required and directory tree CANNOT be changed.

**extract_counts.py:** Generate an additional dataframe with peptide properties to be used in R scripts and the next function.
Packages concurrent.futures and biopython are required. Parallelization is necessary because the regular process was too slow.

**get_covariance_aa.py:** Generate the matrices with the theta score per position (5 matrices 4x16 and 10 matrices 16x16).
If executed, the matrices used in the previous version will be changed.

## Extra scripts

The following files are mostly scripts generated to support the analysis and some extra functions to make reporting easier:

**data_management.R:** Source file that summarizes the data and removes invalid amino acids and nucleic acids for data analysis.

**analyse_counts.R:** Script used to evaluate singleton impact, normalization strategies and initial amino acids covariance (later implemented in get_covariance_aa.py).

**contact_map_analysis.R:** Uses protein/DNA contact maps data to define contact regions. Covariance between amino acids and DNA was calculated here but as it did not show good results, implementation was not transfered to python. Additionally, terminal notation for PYMOL is available to facilitate structure representation.

**logos/automate_logo.py:** Implemented a web scraper to automate the logos generation using [enoLOGOS](http://www.benoslab.pitt.edu/cgi-bin/enologos/enologos.cgi). Packages selenium and webdriver_manager are required. The execution is too slow at the moment because of the response times of the website.

**logos/visualize_logos_chunk.R, visualize_logos.Rmd, visualize_logos_allChunks.Rmd and visualize_logos_allChunks2.Rmd:** Nested source files that generate an HTML output of the logos for comparison. The nesting is required because R markdown chunks do not support loops.

## License

MIT License: http://adampritchard.mit-license.org/ or see the [LICENSE](LICENSE) file.
