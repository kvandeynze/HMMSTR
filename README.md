# HMMSTR
## A modified profile-HMM for tandem repeat copy-number determination from long-reads
HMMSTR calls tandem repeat copy number from raw, long-read, sequencing reads and reports copy numbers in both a read and sample specific format. While designed to model Nanopore sequencing errors in repetitive regions, HMMSTR can be applied to PacBio data as well and has flexible arguments to allow for custom error rates.

HMMSTR is optimized for targeted sequencing experiments and can be run with a single or multiple target regions/sequences in a global-alignment and reference free format.

## Dependencies
* Python >= 3.8
* colorama
* numpy
* pandas
* pickleshare
* scikit-learn
* scipy
* seaborn
* importlib-resources
* mappy

## Installation
HMMSTR is available on Pypi and Conda**
```
pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple HMMSTR==0.1.3
```
(conda install here)
(git clone install here)

## Usage
```
usage: hmmstr [-h] [--background BACKGROUND] [--E_probs E_PROBS] [--A_probs A_PROBS] [--custom_RM CUSTOM_RM] [--hmm_pre HMM_PRE] [--output_hist] [--max_peaks MAX_PEAKS] [--cpus CPUS] [--flanking_size FLANKING_SIZE]
              [--mode MODE] [--cutoff CUTOFF] [--k K] [--w W] [--use_full_read] [--peakcalling_method PEAKCALLING_METHOD] [--bandwidth BANDWIDTH] [--kernel KERNEL] [--bootstrap] [--call_width CALL_WIDTH]
              [--resample_size RESAMPLE_SIZE] [--allele_specific_CIs] [--allele_specific_plots] [--discard_outliers] [--filter_quantile FILTER_QUANTILE] [--cluster_only] [--save_intermediates]
              {targets_tsv,coordinates} ... out inFile
```

HMMSTR has 2 input modes:
1. ```targets_tsv```: Directly uses an input tsv with the following columns:
    1. name: the names of all targets for a given run
    2. prefix: the sequence directly upstream of the target repeat (200bp recommended)
    3. repeat: repeat motif for given target (must be on the same strand as prefix and suffix sequences)
    4. suffix: the sequence directly downstream of the target repeat (200bp recommended)
2. ```coordinates```: Uses bedfile to create the targets tsv from the following columns:
    1. Chromosoms
    2. Start coordinate
    3. End coordinate
    4. Repeat motif (on same strand as reference genome)
    5. Target name (optional, if not given name will be assigned as chr:start-end)

```coordinates``` also requires the following additional positional arguments:
1. chrom_sizes: path to chromosome sizes file corresponding to the reference genome used
2. ref: path to the reference genome to get flanking sequences from
3. input_flank_length: Lenth of the prefix and suffix to get from the reference genome, must be longer than 30bp (Default) or the ```--flanking_size``` optional parameter, (optional, default: 200)

### Required Positional Arguments
|  Argument &nbsp; &nbsp; &nbsp; | Description |
|---|---|
|out| Path to output directory with prefix for all results in this run|
|inFile| Sequence file in fasta or fastq format, can be gzipped (must end with .gz)|

### Optional Arguments
| Argument &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;| Description |
|---|---|
|--cpus| Maximum number of CPUs to use during read processing step (default: half of available CPUs)|
|--use_full_read| If passed, HMMSTR will use the full read sequence to predict copy number instead of subsetting each read based on flanking sequence alignment. Optimal for runs where the repeat is close to the end or start of the reads consistently (ie when running on PCR products where primers are relatively close to the repeat of interest)|
| --save_intermediates | Flag designating to save intermediate files including model inputs, raw count files, and state sequence files. NOTE: raw count files are required to recall alleles without rerunning the counting algorithm, see ```--cluster_only```|
| --cluster_only | Only run peak calling step on existing raw repeat copy counts data ```'out''target_name'_counts.txt```. NOTE: Must use the same output and target names as the run that produced the counts files.|

### Model Size
|  Argument &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;  &nbsp; &nbsp; &nbsp;  &nbsp; &nbsp; &nbsp;| Description |
|---|---|
|--flanking_size| Integer designating the number of bases flanking the repeat to encode in the model. Must be shorter or equal in length to given prefix and suffix. Note: significant increases in flanking size will increase runtime but may increase accuracy. Longer flanking sequences are recommended for more repetitive flanking regions or regions with high similarity with respect to sequence directly flanking the repeat (default: 30, 50-100 recommended for highly repetitive regions)|


### Alignment Options
|  Argument &nbsp; &nbsp; &nbsp; | Description |
|---|---|
|--mode| Mode used by mappy. map-ont (Nanopore), pb (PacBio), or sr (short accurate reads, use for accurate short flanking sequence input) (default: map-ont)|
|--cutoff| MapQ cutoff for prefix and suffix (default: 30, range: 0-60)|

### Peak-calling Options
|  Argument &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;| Description |
|---|---|
|--max_peaks| Integer designating the maximum number of alleles to call for a given run (default: 2)|
|--peakcalling_method| Used to override the default peak calling pipeline. Options include: gmm, kde, kde_throw_outliers (default:auto, HMMSTR chooses the best method based on the distribution of copy number per read)|
|--discard_outliers| If passed, outliers will be discarded based on quantile. If ```--filter_quantile``` not set, reads exceding the top and bottom quantile (0.25) will be discarded and marked as outliers in outputs|
|--filter_quantile| Float designating quantile of count frequency to discard when filtering outliers (default: 0.25)

#### KDE Options
|  Argument &nbsp; &nbsp; &nbsp; | Description |
|---|---|
| --bandwidth | Bandwidth to use for KDE. It is recommended to use the default scott method, especially when there is no expectation for the distribution of repeat lengths.|
| --kernel | Kernel to use for the KDE. Default is gaussian, allows for other kernels if testing different distributions is desired.|


### Output Options
|  Argument &nbsp; &nbsp; &nbsp; | Description |
|---|---|
| --output_hist | output supporting reads histogram showing how many reads were assigned to each repeat copy number per target in a single run|
| --bootstrap | Boolean designating to output bootstraped confidence intervals for allele calls. By default, the samples are drawn from the full dataset regardless of allele.|

#### Bootstrapping Options
|  Argument &nbsp; &nbsp; &nbsp;  &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;| Description |
|---|---|
|--call_width| Decimal percentage designating confidence interval width to calculate in bootstrapping (default: 0.95)|
|--resample_size| Number of times to resample the repeat copy number distribution during bootstrapping (default:x)|
|--allele_specific_CIs| Output allele-specific bootstrapped confidence intervals. This process separates data by assigned alleles before sampling.|
|--allele_specific_plots| Output allele-specific histograms with model of best fit|

### Advanced Options
#### Custom Model Parameter Options
Optional tsv inputs to set custom model parameters.
|  Argument &nbsp; &nbsp; &nbsp; | Description |
|---|---|
|--background| TSV with custom background frequencies to encode in genome states (see "background_example.tsv")|
|--E_probs| TSV with custom emission probabilities to be encoded in match states. These should correspond to the expected mismatch rate (see "emission_example.tsv")|
|--A_probs| TSV with custom transition probibilities to be encoded in the model. Column names in "P_xy" format such that 'x' is the first state type and 'y' is the state type 'x' transitions to (see "transition_example.tsv")|
|--custom_RM| TSV with columns corresponding to a given postion in the repeat motif and rows corresponding to possible nucleotides (and deletion character ''). This is used to designate custom nucleotide occupancy per position in a given motif in case of known mosaicism (ie AAGGG vs AAAAG at the CANVAS locus). Note: this matrix will be applied to all models in a given run, it is advised you only use it in single target runs (see "RM_example.tsv")|

#### Advanced Alignment Options
Parameters to pass to Mappy during alignment step
|  Argument &nbsp; &nbsp; &nbsp; | Description |
|---|---|
|--k| Integer designating kmer parameter to be passed to mappy (see mappy documentation)|
|--w| Window parameter to be passed to mappy (see mappy documentation)|

## Example Use Cases
### Basic Use: Single Plasmid Target
Here, we run HMMSTR on a sequence file containing nanopore reads from a plasmid construct with variable copies of an AAAAG repeat motif. Since these are plasmid contructs, we wrote our input tsv file ```AAAAG_input.txt``` by setting the prefix column to the 200bp upstream sequnce directly flanking the AAAAG repeat from the known backbone sequence and set the suffix column with the downstream flanking sequence. For this example, we will use all default parameters with the exception of ```--output_hist``` and ```--max_peaks```.
```
hmmstr targets_tsv AAAAG_input.txt ./tutorial_1 AAAAG_11012021_3000_sample.fasta --max_peaks 3 --output_hist
```
Outputs:
1. ```tutorial_1_genotype_calls.tsv```: TSV containing final allele calls per target (see detailed output section for column descriptions)
2. ```tutorial_1read_assignments.tsv```: TSV containing read level copy number predictions and allele assignments
3. ```tutorial_1_AAAAG_final_out.tsv```: TSV containing additional read-level stats reported from viterbi algorithm, can be used to make custom plots if desired. This file is produced for each input target.
4. ```tutorial_1_AAAAG_context_labeled.txt```: Text file contianing repeat sequence and flanking context sequence colored by the optimal state path along with the read name and strand. This can be viewed on the command line. This is helpful when determining if the prefix/suffix you inputted are well fit to the repeat of interest and can help in debugging your inputs. This file is produced for each input target.
5. ```tutorial_1AAAAGpeaks.pdf```: (Optional) Supporting read histogram displayed with the model of best fit as a density plot -- GMM or KDE depending on the peak caller chosen.
6. ```tutorial_1AAAAGAIC_BIC.pdf```: (Optional) If GMM chosen, the AIC and BIC are plot and outputted here. These metrics are used to determine the most likely number of clusters.
7. ```tutorial_1AAAAG_supporting_reads_hist.pdf```: (Optional) Raw supporting read histogram, copy number by number of supporting reads.
Below is an example of the *context_labeled.txt files:
![context labeled example](images/AAAAG_example_context_labelled.jpg)
* Red rectangles represent deletions, green represents insertions, bases labeled as in the repeat sequence are white and the prefix and suffix are in grey
The following plots are produced by the given command:
* Supporting read histogram (7)
![AAAAG example supporting read histogram](images/tutorial_1AAAAG_supporting_reads_hist.jpg)
* Model of best fit -- GMM (5)
![AAAAG example model of best fit](images/tutorial_1AAAAGpeaks.jpg)
* AIC/BIC plot (6)
![AAAAG example AIC/BIC](images/tutorial_1AAAAGAIC_BIC.jpg)

If the same command is run with the KDE ```--peakcalling_method``` option, the model of best fit plot would be the following:
```
hmmstr targets_tsv AAAAG_input.txt ./tutorial_1 AAAAG_11012021_3000_sample.fasta --max_peaks 3 --output_hist --peakcalling_method kde
```
![KDE model of best fit](images/tutorial_1_kdeAAAAG_KDE.jpg)

### Including allele specific output plots and confidence intervals
HMMSTR also includes options to visualize per-read copy number prediction distributions in an allele-specific format. Below is how we would use HMMSTR to output these plots as well as allele-specific confidence intervals. Note: these confidence intervals are produced by bootstrapping the median of a given allele with 100 resamples.
```
hmmstr targets_tsv AAAAG_input.txt ./tutorial_1 AAAAG_11012021_3000_sample.fasta --output_hist --max_peaks 3 --bootstrap --resample_size 100 --allele_specific_CIs --allele_specific_plots
```
Allele 1           |  Allele 2           |  Allele 3
:-------------------------:|:-------------------------:|:-------------------------:
![](images/tutorial_1_allele_specificAAAAGallele_1.jpg)  |  ![](images/tutorial_1_allele_specificAAAAGallele_2.jpg)  |  ![](images/tutorial_1_allele_specificAAAAGallele_3.jpg)
(30.0, 30.0) | (58.0, 59.0) | (16.0, 16.0)
