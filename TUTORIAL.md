# Step-by-step tutorial to run bbcontacts from a single sequence input

This is a step-by-step tutorial to run bbcontacts starting from a single sequence input.

The tutorial covers the following steps:

   * STEP 0: Install required software
   * STEP 1: Build a multiple sequence alignment
   * STEP 2: Run a direct coupling prediction
   * STEP 3: Generate a secondary structure prediction
   * STEP 4: Run bbcontacts
   * Troubleshooting / FAQ

**At each step, the expected output is available (for comparison purposes) in the exampleresults/ directory of the git repository.**

## STEP 0: Install required softare

### STEP 0.1: Install necessary software to run bbcontacts

To run bbcontacts, you will need

   * Python version 2.6 or 2.7 or 3.4,
   * the Python packages for scientific computing [NumPy](http://www.numpy.org/) and [SciPy](http://www.scipy.org/scipylib/index.html).

bbcontacts has been tested on Ubuntu 14.04 and Scientific Linux 6.5.

On Debian-based systems (including Ubuntu), you can install Python and the scientific computing packages using

          sudo apt-get install python
          sudo apt-get install python-numpy python-scipy

On Fedora and related systems, you can install Python and the scientific computing packages using

          sudo yum install python
          sudo yum install numpy scipy

### STEP 0.2: Install the software used to create input for bbcontacts

The following software needs to be installed to generate input files for bbcontacts

   * [the HHsuite](ftp://toolkit.genzentrum.lmu.de/pub/HH-suite/hhsuite-latest.tar.gz) [1] for generation of multiple sequence alignments (takes a single sequence as input)
   * to generate multiple sequence alignments for CCMpred and bbcontacts, the [uniprot20 database for HHblits](ftp://toolkit.genzentrum.lmu.de/pub/HH-suite/databases/hhsuite_dbs/uniprot20_latest.tar.gz) is also required
   * [CCMpred](https://bitbucket.org/soedinglab/ccmpred/) [2] for prediction of direct couplings (takes a multiple sequence alignment as input)
   * [PSIPRED](http://bioinf.cs.ucl.ac.uk/psipred/) [3] for secondary structure prediction (takes a multiple sequence alignment as input)
   * (only if you want to run bbcontacts using DSSP annotation as an input, e.g. for benchmarking purposes) [DSSP](http://swift.cmbi.ru.nl/gv/dssp/) [4] for secondary structure assignment (takes a PDB structure as input)
   
### STEP 0.3: Clone the git repository for bbcontacts.

Simply clone the [git repository for bbcontacts](https://bitbucket.org/soedinglab/bbcontacts).


## STEP 1: Build a multiple sequence alignment

### STEP 1.1: HHblits

First, run HHblits against the uniprot20 database, avoiding any filtering in order to retrieve as many homologous sequences as
possible (adapting the path to the uniprot database as needed):

          hhblits -i example/1nz0D.fas -oa3m example/1nz0D.a3m -all -maxfilt 100000 -realign_max 100000 -B 100000 -Z 100000 -d HHblits_databases/uniprot20_2013_03/uniprot20_2013_03

### STEP 1.2: HHfilter

Then, perform a filtering step using HHfilter:

          hhfilter -i example/1nz0D.a3m -o example/1nz0D.filt.a3m -id 90 -neff 15 -qsc -30

### STEP 1.3: Reformat the output alignment

In order to run CCMpred, the output alignment must be reformatted to the "PSICOV" format used by CCMpred,
using first the reformat.pl script from the HHsuite (to get an alignment in fasta format)
and then the convert_alignment.py script provided in the scripts directory of the CCMpred git repository:

          reformat.pl example/1nz0D.filt.a3m example/1nz0D.filt.fas -r
          convert_alignment.py example/1nz0D.filt.fas fasta example/1nz0D.filt.psc


## STEP 2: Run a direct coupling prediction

Direct coupling predictions are obtained with CCMpred run with the default options:

          ccmpred example/1nz0D.filt.psc example/1nz0D.filt.mat

The exact numerical values obtained for the couplings might depend on your system characteristics,
in particular whether you have a GPU that CCMpred can use to compute direct coupling predictions.


## STEP 3: Generate a secondary structure prediction

The secondary structure predictions are obtained with PSIPRED [3], as implemented in the addss.pl script from the HHsuite:

          addss.pl example/1nz0D.filt.a3m  example/1nz0D.filt.withss.a3m -a3m
          head -2 example/1nz0D.filt.withss.a3m > example/1nz0D.filt.psipred


## STEP 4: Run bbcontacts

### Input and options

To run bbcontacts for a given protein, you will need

   * the matrix of predicted couplings: example/1nz0D.filt.mat
   * the diversity value of the alignment used to compute the contact predictions, defined as sqrt(N)/L where N is the number of sequences in the alignment and L is the protein length: in the case of 1nz0D, the protein chain has length L=111 and the alignment contains 1741 sequences so the diversity value is 0.376
   * a 3-state secondary structure assignment (-d option) or prediction (-p option): in this case the PSIPRED prediction example/1nz0D.filt.psipred

Here is an example command line to run bbcontacts for the provided example from the input files we just computed:
       
          bbcontacts example/1nz0D.filt.mat 0.376 1nz0D -p example/1nz0D.filt.psipred

You will obtain an output file 1nz0D.filteredoutput.txt containing the bbcontacts predictions.
To compare the output you obtained to the expected output, you can then run:

          diff 1nz0D.filteredoutput.txt exampleresults/1nz0D.filteredoutput.txt
          vimdiff 1nz0D.filteredoutput.txt exampleresults/1nz0D.filteredoutput.txt

The Viterbi score values might differ slightly due to the numerical differences in the CCMpred computed couplings.

As the example provided here (1nz0D) is a chain from the BetaSheet916 dataset [5], you can also download the [dataset annotation file](http://contact.ics.uci.edu/betadata/betasheet.dat) and run:
       
          bbcontacts example/1nz0D.filt.mat 0.376 1nz0D -p example/1nz0D.filt.psipred -n 1nz0D -e betasheet.dat

In addition to the 1nz0D.filteredoutput.txt output file, you will obtain an additional output file 1nz0D.filteredoutput.eval.txt that contains two extra columns: a T/F column telling whether the contact is a TP or FP according to the evaluation file and a column containing the total number of true residue-level contacts for this PDB chain (useful for recall calculation).

Note: if the output file already exists (e.g. if you run bbcontacts.py twice on the same input, without removing the output in between), then the new output will be appended at the end of the .filteredoutput.txt file.

Note: If no pattern is detected above the Viterbi score thresholds specified in the configuration file, the output file contains only one line where all output values are NA.


## Troubleshooting

### I get errors when trying to run hhblits, hhfilter or addss.pl.

The most likely cause for this is that you have not properly installed and set up the HHsuite.

   * Make sure you have set the environment variable HHLIB to $INSTALL_DIR/lib/hh
   * Make sure you have specified the paths to BLAST and PSIPRED (used by addss.pl) in $INSTALL_DIR/lib/hh/scripts/HHPaths.pm

Note: In order to run the present bbcontacts pipeline, you do not need the BLAST databases, since addss.pl uses only some tools from the BLAST suite. However, you need the ["legacy" version of BLAST](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/) (i.e. not BLAST+), which is also the preferred version for PSIPRED.

### My CCMpred prediction is different from the one given in exampleresults.

The exact numerical values obtained for the couplings might depend on your system characteristics, in particular whether you have a GPU that CCMpred can use to compute direct coupling predictions.
However, the values you obtain should be of the same order as those in exampleresults and the matrix should be of the same size.

### My PSIPRED prediction is different from the one given in exampleresults.

This is probably because you are using a different version of PSIPRED. The files provided in exampleresults have been obtained with PSIPRED v2.6 as the parameters in addss.pl were optimized for PSIPRED V2.

### When trying to run bbcontacts, I get the error "Problem with query: The size of the coupling matrix does not match the length of the secondary structure sequence!"

The most likely cause for this error is that you have run CCMpred without first converting the output alignment from the HHsuite (in a3m format) to the "PSICOV" format required by CCMpred. Please make sure you have executed the commands listed in STEP 1.3 to convert the alignment.

### My bbcontacts output file is different from the one given in exampleresults.

The bbcontacts results (in particular the Viterbi scores) depend on the exact values of the couplings obtained with CCMpred. Therefore, if your CCMpred prediction is a little different from the one in exampleresults, your bbcontacts output file should also be a little different from the one in exampleresults. Most likely, this will affect only the Viterbi score values and not the predicted pairs of residues in contact or the order of the predictions.


## Citation
If you use bbcontacts, please cite:

J. Andreani and J. Soeding (2015) "bbcontacts: prediction of beta-strand pairing from direct coupling patterns", Bioinformatics. [doi:10.1093/bioinformatics/btv041](http://dx.doi.org/10.1093/bioinformatics/btv041)


## Contact

[Jessica Andreani](mailto:jessica.andreani@mpibpc.mpg.de)

[Johannes Soeding](mailto:johannes.soeding@mpibpc.mpg.de)


## References

   [1] M. Remmert, A. Biegert, A. Hauser and J. Soeding (2012) HHblits: lightning-fast iterative protein sequence searching by HMM-HMM alignment. Nature Methods, 9, 173-175.

   [2] S. Seemayer, M. Gruber and J. Soeding (2014) CCMpred --- fast and precise prediction of protein residue-residue contacts from correlated mutations. Bioinformatics, doi: 10.1093/bioinformatics/btu500.

   [3] D. T. Jones (1999). Protein secondary structure prediction based on position-specific scoring matrices. J. Mol. Biol., 292, 195-202.

   [4] W. Kabsch and C. Sander (1983). Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. Biopolymers, 22, 2577-2637.

   [5] J. Cheng and P. Baldi (2005). Three-stage prediction of protein beta-sheets by neural networks, alignments and graph algorithms. Bioinformatics, 21 Suppl 1, 75-84. [Link to the BetaSheet916 dataset file](http://contact.ics.uci.edu/betadata/betasheet.dat) (version used in the paper accessed 12 September 2014)


