# Predictability of antigen binding based on short motifs in the antibody CDRH3

This GitHub repository contains code, analysis specifications and documentation for the manuscript titled "Predictability of antigen binding based on short motifs in the antibody CDRH3".

It is often assumed that the antigen binding status of an immune receptor (antibody/TCR) may be determined by the _presence_ of a short motif in the CDR3. 
We directly test this assumption using a new method for the discovery of high-precision, highly generalizable motifs. Our method is verified using experimental as well as simulated data, and is integrated into the AIRR-ML platform [immuneML](https://github.com/uio-bmi/immuneML). 


## Data preparation

### Peparing experimental data

Two experimental datasets are used: a previously publised dataset ([Mason, 2021](https://www.nature.com/articles/s41551-021-00699-9)) and the yet-to-be published Mehta dataset (manuscript in preparation).
The folder [preprocessing](source/preprocessing) contains two python scripts for preprocessing the data (removing duplicates and removing overlap between the binder/non-binder classes) and splitting the data into 50% training, 25% validation and 25% test data.

### Simulating data with and without ground-truth motifs

To test the performance of our method, two simulated datasets are used.
One simulated dataset contains ground-truth motifs for predicting antigen binding, meaning that a sequence assigned to the 'binder' class if (and only if) at least one ground-truth motif is present ([simulated_data](data/simulated_data)). 
The second dataset contains no ground-truth motifs, meaning binding is random and cannot be learned ([simulated_no_motifs](data/simulated_no_motifs)). 
The code for simulating synthetic datasets can be found in the folder [simulation](source/simulation). 



## Main analysis with immuneML

Six immuneML runs were performed to generate all manuscript figures and data. 
The YAML analysis specifications are located in the folder [yaml_files](immuneml_yaml_files).
The original immuneML results folders (excluding data) can be found at the following doi: [10.5281/zenodo.10047590](https://doi.org/10.5281/zenodo.10047590).

To rerun any of the immuneML analyses, make sure to [install immuneML](https://docs.immuneml.uio.no/latest/installation/install_with_package_manager.html). 
When applying the CNN proposed by Mason et al., (KerasSequenceCNN in immuneML) please ensure that the dependencies `keras` and `tensorflow` are installed (e.g., version 2.11.0).
These dependencies are optional for immuneML and are not automatically installed when installing immuneML.

immuneML can be run as follows:

```bash
immune-ml yaml_file.yaml results_folder
```

For mode details, see the [immuneML documentation](https://docs.immuneml.uio.no/latest/). 


### Simulated data analysis

First, a set of motifs with high precision on the training set are learned, and the validation set is used to determine optimal recall thresholds ([determine_recall_simulated.yaml](immuneml_yaml_files/determine_recall_simulated.yaml)). 
In the dataset without ground-truth motifs, the optimal recall threshold cannot be determined. 

The learned recall thresholds are used as input parameters for the next immuneML run ([motif_performance_simulated.yaml](immuneml_yaml_files/motif_performance_simulated.yaml)). 
Here, a new set of motifs is learned with high precision on the training + validation set, and recall scores that exceed the learned thresholds. 
In this immuneML run, only the dataset with ground-truth implanted motifs is used, and several reports are run to gain additional insight into motif performance and properties, such as plotting how well the learned motifs overlap with the ground-truth motifs. 

The positional amino acid distribution of the simulated dataset with ground-truth motifs is furthermore plotted in a separate immuneML run ([aa_distributions.yaml](immuneml_yaml_files/aa_distributions.yaml)).

### Experimental data analysis

Similarly to the simulated data analysis, a separate immuneML run is first done to determine the optimal recall thresholds on the Mason dataset ([determine_recall_mason.yaml](immuneml_yaml_files/determine_recall_mason.yaml)).
These thresholds are used as input parameters in the following immuneML run ([ml_mason.yaml](immuneml_yaml_files/ml_mason.yaml)), where several motif-based classifiers, sequence similarity-based classifiers and the CNN by Mason et al. are trained to predict antigen binding on the Mason dataset. This immuneML run furthermore includes the performance of individual learned motifs on the Mason training + validation set. 

Lastly, a separate immuneML run was made to assess the performance of each trained classifier, as well as individual motifs, on the separate Mehta test set ([ml_application_mehta.yaml](immuneml_yaml_files/ml_application_mehta.yaml)). 
Note that these classifiers and motifs were trained exclusively on the Mason training+validation dataset, not the Mehta dataset. 

The positional amino acid distribution of both Mason and Mehta datasets are plotted in a separate immuneML run ([aa_distributions.yaml](immuneml_yaml_files/aa_distributions.yaml)).

## Collecting results

A Python script was created to parse the immuneML results folders and retrieve only the necessary results files ([get_manuscript_figures_to_keep.py](source/figures_and_tables/get_manuscript_figures_to_keep.py)), 
in order to keep only the figures and tables that were used in the manuscript. The results of this script are stored in the folder [sorted_results](sorted_results).
Some additional supplementary results are generated by other scripts in the folder [figures_and_tables](source/figures_and_tables). 



