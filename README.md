# refract-chile
Repository to deposit relevant code generated during my trip to Santiago through the REFRACT program between Jan/April 2023.

The following files are available:

**get_peptide_score.py:** Main function to calculate the score for the provided peptide with 6 standard amino acids.
In position 0 only C, F, W and Y are accepted and they are not accepted in the other 5 positions.
At the moment accepts just 1 peptide as an entry and can be run through the terminal with the following code:

```
python get_peptide_score.py -i <peptide>
```

*Example:* python get_peptide_score.py -i FEKDAD

The following scripts must be executed IN ORDER to re-generate all outputs:
Source files are required and directory tree CANNOT be changed.

**extract_counts.py:** Generate an additional dataframe with peptide properties to be used in R scripts and the next function.
Packages concurrent.futures and biopython are required.

**get_covariance_aa.py:** Generate the matrices with the theta score per position (5 matrices 4x16 and 10 matrices 16x16).
If executed, the matrices used in the previous version will be changed.
