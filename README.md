# HLAProfiler

## Synopsis

HLAProfiler uses the k-mer content of next generation sequencing reads to call HLA types in a sample. Based on the k-mer content each each read pair is assigned to an HLA gene and the aggregate k-mer profile for the gene is compared to reference k-mer profiles to determin the HLA type. Currently HLAProfiler only supports paired-end RNA-seq data.

## Code Example

Database creation:

```
perl ~/opt/scripts/HLAProfiler/bin/HLAProfiler.pl build -t path/to/transcript.fa.gz -g /path/to/transcript_annotations.gtf.gz -e /path/to/hla_exclue_regions.bed -r /path/to/IMGT_reference.fasta -cwd /path/to/hla_cwd_alleles.txt -o /path/to/output_directory/ -db database_name -kp /path/to/kraken/ -c 12 
```

HLA calling:

```
perl HLAProfiler.pl predict -fastq1 sequencedata_1.fastq.gz -fastq2 sequencedata_2.fastq.gz -database_name database -database_dir /path/to/HLAProfiler -reference /path/to/HLAProfiler/database/data/reference/hla.ref.merged.fa -threads 8 -output_dir /path/to/output/dir -allele_refinement all -kraken_path /path/to/kraken/ -if -l sample.HLAProfiler.log
```

## Motivation

Accurate HLA calling in RNA-seq data, extends the utility and increases the information of this data. HLAProfiler fills the gap of existing tools by improving accuracy for some HLA genes, identifying rare alleles, and correctly identifying novel alleles and alleles with incomplete reference data.

## Installation

For automatic installation, download and unzip the HLAProfiler GitHub release tarball or clone the git repository into the directory where you want the package to be installed ("/path/to/HLAProfiler")
```
cd /path/to/HLAProfiler
perl bin/install.pl –d /path/to/install_dir/ -m wget –b /bin/directory/in/PATH/
```
The installer will retrieve the dependencies and place them in /bin/directory/in/PATH

Database download:
The database is available as a part of the repository. Additionally two binary files are included with the latest GitHub release and need to be downloaded and moved into the database folder.

## Tests
```
cd /path/to/HLAProfiler/ test 
perl run_tests.pl –t all –kp /path/to/kraken
```
If using the automated installer kraken (/path/to/kraken) will be located in bin/kraken-version/ under the HLAProfiler directory. 

## Contact
Martin Buchkovich (martin.buchkovich@q2labsolutions.com)

## License
Use of this software is limited to non-commercial use under Q Square Solutions Expression Analysis user license. See the [LICENSE](/LICENSE.md) file for more details. 
