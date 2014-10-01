# bbcontacts
** Prediction of protein beta-beta contacts at the residue level using direct coupling patterns **

bbcontacts is a Python program predicting residue-level contacts between beta-strands by detecting patterns in matrices of predicted couplings.
bbcontacts can make use of a secondary structure assignment or a secondary structure prediction.


## Requirements

### Necessary software

To run bbcontacts, you will need

   * Python version 2.6 or 2.7,
   * the Python packages for scientific computing [NumPy](http://www.numpy.org/) and [SciPy](http://www.scipy.org/scipylib/index.html).

bbcontacts has been tested on Ubuntu 14.04 and Scientific Linux 6.5.

### Recommended software to create input for bbcontacts

The following software is not necessary to run bbcontacts itself, but it is useful to generate input files for bbcontacts

   * [HHblits](http://toolkit.genzentrum.lmu.de/hhblits/) [1] for generation of multiple sequence alignments (takes a single sequence as input)
   * [CCMpred](https://bitbucket.org/soedinglab/ccmpred/) [2] for prediction of direct couplings (takes a multiple sequence alignment as input)
   * [DSSP](http://swift.cmbi.ru.nl/gv/dssp/) [3] for secondary structure assignment (takes a PDB structure as input)
   * [PSIPRED](http://bioinf.cs.ucl.ac.uk/psipred/) [4] for secondary structure prediction (takes a multiple sequence alignment as input)


## Installation
bbcontacts is written in Python, therefore no installation is required. Simply clone the git repository and run bbcontacts.py.


## License
bbcontacts is released under the GNU Affero General Public License v3 or later. See LICENSE for more details.


## Usage

     Usage: bbcontacts.py [options] couplingmatrix diversityvalue outputprefix (-d DSSPsecstructfile | -p PSIPREDsecstructfile)

### Input and options

To run bbcontacts for a given protein, you will need

   * a matrix of predicted couplings,
   * the diversity value of the alignment used to compute the contact predictions, defined as sqrt(N)/L where N is the number of sequences in the alignment and L is the protein length,
   * a 3-state secondary structure assignment (-d option) or prediction (-p option).

More specifically, for bbcontacts to give the best possible result, you should ensure that

   * the matrix of predicted couplings has been generated with CCMpred [2] using the default options: 50 iterations, sequence reweighting factor of 0.8, pairwise regularization pre-factor of 0.2 and with average product correction (APC). Note that because when CCMpred performs the APC step, the minimum coupling value gets subtracted from all coupling values, you should either make sure to use a smoothing range when running bbcontacts (-s option with a non-0 value; default is -s 10) or you should generate the raw prediction matrix from CCMpred and recompute the APC-corrected couplings.
   * the secondary structure assignment/prediction is given in the required (FASTA-like) format, with only 3 states (H, E and C). A DSSP [3] assignment mapped to 3 states or a PSIPRED [4] prediction are recommended, since bbcontacts was trained from DSSP assignments and PSIPRED predictions.

Additionally, 

   * the -c option allows you to specify a configuration file (by default, the config file is assumed to be the bbcontacts.conf file in the working directory),
   * the -s option allows you to specify a smoothing range (by default, -s 10 is used since it was found to give the best performance; the use of a larger value for the smoothing range (e.g. 20 or 30) is discouraged since it degrades the results and increases the runtime; -s 0 turns off the smoothing),
   * the -l option allows you to turn off the prediction-shortening mode (this mode is on by default but will only apply when long predictions occur; in this case, it will shorten the predictions to a reasonable length which improves the overall performance, but it will also slow down the bbcontacts prediction process),
   * the -n option allows you to specify an identifier for the run, which will be used in the output files,
   * the -e option allows you to provide an evaluation file, in the format provided for either dataset BetaSheet916 [5] or dataset BetaSheet1452 [6]. In this case, the identifier provided with the -n option is required and will be searched for in the evaluation file.

Here are two example command lines to run bbcontacts for the provided example:
       
       ./bbcontacts.py example/1nz0D.mat 0.376 1nz0D -p example/1nz0D.psipred
       ./bbcontacts.py example/1nz0D.mat 0.376 1nz0D -d example/1nz0D.dssp -c bbcontacts.conf -s 10 -l -n 1nz0D

### Output

When you run bbcontacts, you have to specify the output prefix as the third positional argument. Several output files will be generated:

   * outputprefix.filteredoutput.txt (or, if a DSSP assignment is specified, outputprefix.DSSP.filteredoutput.txt): this file contains the beta-contact predictions filtered to retain only predictions satisfying the topological constraints of beta-contacts (see paper for details),
   * outputprefix.filteredoutput.eval.txt (or, if a DSSP assignment is specified, outputprefix.DSSP.filteredoutput.eval.txt): this file is only generated if you gave an evaluation file (-e option) and an identifier (-n option) that was found in this file. This file has the same content as the filtered output file, plus two additional columns: a T/F column telling whether the contact is a TP or FP according to the evaluation file and a column containing the total number of true residue-level contacts for this identifier.


## Citation
If you use bbcontacts, please cite:

J. Andreani and J. Soeding (2014) "bbcontacts: predicting beta-strand interactions from direct coupling patterns" (submitted)


## Contact

[Johannes Soeding](mailto:johannes.soeding@mpibpc.mpg.de) (corresponding author)

[Jessica Andreani](mailto:jessica.andreani@mpibpc.mpg.de)


## References

   [1] M. Remmert, A. Biegert, A. Hauser and J. Soeding (2012) HHblits: lightning-fast iterative protein sequence searching by HMM-HMM alignment. Nature Methods, 9, 173-175.

   [2] S. Seemayer, M. Gruber and J. Soeding (2014) CCMpred --- fast and precise prediction of protein residue-residue contacts from correlated mutations. Bioinformatics, doi: 10.1093/bioinformatics/btu500.

   [3] W. Kabsch and C. Sander (1983). Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. Biopolymers, 22, 2577-2637.

   [4] D. T. Jones (1999). Protein secondary structure prediction based on position-specific scoring matrices. J. Mol. Biol., 292, 195-202.

   [5] J. Cheng and P. Baldi (2005). Three-stage prediction of protein beta-sheets by neural networks, alignments and graph algorithms. Bioinformatics, 21 Suppl 1, 75-84. [Link to the BetaSheet916 dataset file](http://contact.ics.uci.edu/betadata/betasheet.dat) (last accessed 12 September 2014)

   [6] C. Savojardo, P. Fariselli, P. L. Martelli and R. Casadio (2013). BCov: a method for predicting beta-sheet topology using sparse inverse covariance estimation and integer programming. Bioinformatics, 29, 3151-3157. [Link to the BetaSheet1452 dataset file](http://biocomp.unibo.it/savojard/bcov/BetaSheet1452.dat) (last accessed 12 September 2014)

