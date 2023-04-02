# To develop a DNA identification service:

The following steps could be taken:

Acquire a sequence database: In order to compare and identify the closest sequence to the provided sequence, a sequence database needs to be acquired. The database should contain a diverse range of DNA sequences and be adapted from GEO, which is a public functional genomics data repository.

Preprocessing the database: The acquired sequence database needs to be preprocessed to remove any redundancy or irrelevant sequences. This can be achieved by filtering out sequences that have low quality or are incomplete.

Align the sequences: Once the database has been preprocessed, the next step is to align the sequences. This involves comparing the test sequence to the sequences in the database to identify the closest sequence(s). This can be done using different alignment algorithms such as BLAST, ClustalW, or MAFFT, MEGA11 software.

Calculate the differences: After aligning the sequences, the differences between the test sequence and the closest sequence(s) can be calculated. This can be done by computing the edit distance or Hamming distance, which measures the number of changes required to transform one sequence into another.

Output the results: The final step is to output the closest sequence and the differences between the test sequence and the closest sequence(s).
Stretch goals:

Probabilities across the database: To calculate the probability of finding a similar sequence to the test sequence in the database, a statistical approach can be used. This involves calculating the p-value, which is a measure of the probability of obtaining a result as extreme as or more extreme than the observed result, assuming that the null hypothesis is true. The p-value can be calculated using various statistical methods such as Fisher's exact test or chi-squared test.

Reconstructed phylogeny: To gain insights into the evolutionary relationships between the sequences in the database, a phylogenetic analysis can be performed. This involves reconstructing the evolutionary history of the sequences and creating a phylogenetic tree that illustrates the relationships between the sequences. This can be done using different methods such as maximum likelihood or Bayesian inference. The reconstructed phylogeny can provide valuable information about the evolutionary origins of the sequences and their functional and structural relationships.

​
Gene Expression Omnibus
GEO is a public functional genomics data repository supporting MIAME-compliant data submissions. Array- and sequence-based data are accepted.

The GEO database of Series GSE90441 has been used to obtain preprocessed fasta sequence file, dog_breeds.fa which is used to identify closest breed associated with query sequence.

# Alignment of sequences

class PairwiseAligner(_algorithms.PairwiseAligner)
 |  PairwiseAligner(scoring=None, **kwargs)
 |  
 |  Performs pairwise sequence alignment using dynamic programming.
 |  
 |  This provides functions to get global and local alignments between two
 |  sequences.  A global alignment finds the best concordance between all
 |  characters in two sequences.  A local alignment finds just the
 |  subsequences that align the best.
 |  
 |  To perform a pairwise sequence alignment, first create a PairwiseAligner
 |  object.  This object stores the match and mismatch scores, as well as the
 |  gap scores.  Typically, match scores are positive, while mismatch scores
 |  and gap scores are negative or zero.  By default, the match score is 1,
 |  and the mismatch and gap scores are zero.  Based on the values of the gap
 |  scores, a PairwiseAligner object automatically chooses the appropriate
 |  alignment algorithm (the Needleman-Wunsch, Smith-Waterman, Gotoh, or
 |  Waterman-Smith-Beyer global or local alignment algorithm).
 |  
 |  Calling the "score" method on the aligner with two sequences as arguments
 |  will calculate the alignment score between the two sequences.
 |  Calling the "align" method on the aligner with two sequences as arguments
 |  will return a generator yielding the alignments between the two
 |  sequences.
 
 Some examples:
 |  
 |  >>> from Bio import Align
 |  >>> aligner = Align.PairwiseAligner()
 |  >>> alignments = aligner.align("TACCG", "ACG")
 |  >>> for alignment in sorted(alignments):
 |  ...     print("Score = %.1f:" % alignment.score)
 |  ...     print(alignment)
 |  ...
 |  Score = 3.0:
 |  target            0 TACCG 5
 |                    0 -|-|| 5
 |  query             0 -A-CG 3
 |  <BLANKLINE>
 |  Score = 3.0:
 |  target            0 TACCG 5
 |                    0 -||-| 5
 |  query             0 -AC-G 3
 |  <BLANKLINE>
    
---------------------------------------------------------------------------------------------------------------------------------    
    
# Alternate Method: Using MEGA11 software
    
MEGA11: Molecular Evolutionary Genetics Analysis Version 11
    
:: DESCRIPTION

MEGA （Molecular Evolutionary Genetics Analysis）is an integrated tool for automatic and manual sequence alignment, inferring phylogenetic trees, mining web-based databases, estimating rates of molecular evolution, and testing evolutionary hypotheses.
    
An aligned output text file can be saved after the run 'output.txt'.
    
# Calculate the Difference
    
After aligning the sequences, the PairwiseAligner package can be used to calculate the differences between the test sequence and the closest sequence(s).
    
# Probabilities across the database:
    
# Statistical functions (scipy.stats)
This module contains a large number of probability distributions, summary and frequency statistics, correlation functions and statistical tests, masked statistics, kernel density estimation, quasi-Monte Carlo functionality, and more.

Statistics is a very large area, and there are topics that are out of scope for SciPy and are covered by other packages. Some of the most important ones are:

statsmodels: regression, linear models, time series analysis, extensions to topics also covered by scipy.stats.

scipy.stats.shapiro
scipy.stats.shapiro(x)[source]
Perform the Shapiro-Wilk test for normality.

The Shapiro-Wilk test tests the null hypothesis that the data was drawn from a normal distribution.

Parameters:
xarray_like
Array of sample data.

Returns:
statisticfloat
The test statistic.

p-valuefloat
The p-value for the hypothesis test.
    
# Matplotlib is a comprehensive library for creating static, animated, and interactive visualizations in Python.

Installation

Install using pip:

pip install matplotlib

Install using conda:

conda install -c conda-forge matplotlib
    
matplotlib.pyplot.figure
matplotlib.pyplot.figure(num=None, figsize=None, dpi=None, *, facecolor=None, edgecolor=None, frameon=True, FigureClass=<class 'matplotlib.figure.Figure'>, clear=False, **kwargs)
Create a new figure, or activate an existing figure.
    
    
# Reconstructed phylogeny:
--------------------------    
    
Bio.Phylo.TreeConstruction module
Classes and methods for tree construction.

classBio.Phylo.TreeConstruction.DistanceMatrix(names, matrix=None)
Bases: Bio.Phylo.TreeConstruction._Matrix

Distance matrix class that can be used for distance based tree algorithms.

All diagonal elements will be zero no matter what the users provide.

__init__(self, names, matrix=None)
Initialize the class.

__setitem__(self, item, value)
Set Matrix’s items to values.

format_phylip(self, handle)
Write data in Phylip format to a given file-like object or handle.

The output stream is the input distance matrix format used with Phylip programs (e.g. ‘neighbor’). See: http://evolution.genetics.washington.edu/phylip/doc/neighbor.html

Parameters
handlefile or file-like object
A writeable file handle or other object supporting the ‘write’ method, such as StringIO or sys.stdout. On Python 3, should be open in text mode.
    

# Bio.Align.Applications package

Module contents
Alignment command line tool wrappers.

classBio.Align.Applications.MuscleCommandline(cmd='muscle', **kwargs)
Bases: Bio.Application.AbstractCommandline

Command line wrapper for the multiple alignment program MUSCLE.

http://www.drive5.com/muscle/
    
    
The final output is:
    
# Input files: 'dog.breeds.fa', 'mystery.fa'
# Closest sequence : gb|AY656744.1|
  mystery sequence:  gb|KM061522.1|
    
# the difference: 24
    
# p-value and phylogenic tree