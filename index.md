## What is it?

HLAProfiler uses the *k*-mer content of next generation sequencing reads to call HLA types in a sample. Based on the *k*-mer content each each read pair is assigned to an HLA gene and the aggregate *k*-mer profile for the gene is compared to reference *k*-mer profiles to determin the HLA type. Currently HLAProfiler only supports paired-end RNA-seq data.

## How does it work?
HLAProfiler utilizes the *k*-mer content of RNA-sequencing reads in two different ways in order to assign HLA types to a sample. First, it utilizes a modifed version of the Kraken taxonomic classifier along with a custom HLA taxonomy to assign reads to a HLA gene. Reads which cannot be assigned exclusively to a single gene are excluded from analysis. Second, HLAProfiler counts the observed *k*-mers in the gene specific reads and compares the *k*-mer profile to a database of expected reference *k*-mer profiles. Using this profile comparison, HLAProfiler identifies the HLA allele pair for each gene with the most supporting evidence in the observed sequence reads.
by simulating paired-end reads from each HLA reference allele without errors, filtering these reads based on *k*-mer content, and tallying *k*-mer counts of the reads filtering to the expected genes. 

Using HLAProfiler requires two steps, 1) database creation and 2) HLA calling.

### Database Creation
For each allele, HLAProfiler simulates paired-end reads without errors, assigns each read to an HLA gene and calculates the *k*-mer profile using reads from the expected gene. As there are thoussand of alleles in the IPD/IMGT reference database, this step is very time intensive. 
For convenience, a ready-to-use, downloadable database is available from one of three sources:
1) [GitHub HLAProfiler database only release](https://github.com/ExpressionAnalysis/HLAProfiler/releases/tag/v1.0.0-db_only). To access the full database, download this GitHub release (including the binaries) and make sure the associated binary files (database.idx and database.kdb) are located in the database directory
2) [Q2 Lab Solutions Genomics Bioinformatics website](http://www.q2labsolutions.com/genomics-laboratories/bioinformatics) (Request data at bottom of the page)
3) Request to the authors by [email](mailto:martin.buchkovich@q2labsolutions.com) 

To create a new database, an HLA reference fasta, transcriptome-wide transcript fa and gtf, and an exclusion bed are required. The transcript files and exclusion bed are used to create the distractome, which helps control for homology between HLA genes and other transcripts. The exclusion bed denotes genomic regions to exclude from the distractome. Any reads assigned to the distractome will be excluded from analysis. While we recommend the IPD/IMGT HLA database and GENCODE as the source of these references, any files can be used as long as they adhere to the following naming and format conventions.

#### Required Naming Conventions

###### Reference HLA fasta:
The sequence identifier must follow the format" >ID<tab>AlleleName
i.e.
```
HLA00001	A*01:01:01:01
```

###### Transcript fasta:
Fields of the sequence identifier must be separated by '|' and one of the fields must match the gene_name tag of the GTF annotation field.
i.e.
```
>ENST00000473358.1|ENSG00000243485.3|OTTHUMG00000000959.2|OTTHUMT00000002840.1|RP11-34P13.3-001|RP11-34P13.3|712|lincRNA|
```

###### Transcript GTF:
Standard [GTF](http://www.ensembl.org/info/website/upload/gff.html) format. The annotation column must contain the gene_name tag.

###### Exclusion Bed:
Standard [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format

####  Example build commands
Once the proper reference files are downloaded or created the database can be built using the following command:
```
perl ~/opt/scripts/HLAProfiler/bin/HLAProfiler.pl build -t path/to/transcript.fa.gz -g /path/to/transcript_annotations.gtf.gz -e /path/to/hla_exclue_regions.bed -r /path/to/IMGT_reference.fasta -cwd /path/to/hla_cwd_alleles.txt -o /path/to/output_directory/ -db database_name -kp /path/to/kraken/ -c number_of_threads 
```
HLA alleles can be processed in parallel and given the thousands of alleles that need processing we recommend specifying multiple threads using the -c option.

For more details on usage run:
```
perl HLAProfiler.pl build -h
```

### HLA calling:
HLAProfiler can be used for HLA calling by running the following example command:
```
perl HLAProfiler.pl predict -fastq1 sequencedata_1.fastq.gz -fastq2 sequencedata_2.fastq.gz -database_name database -database_dir /path/to/HLAProfiler -reference /path/to/HLAProfiler/database/data/reference/hla.ref.merged.fa -threads 8 -output_dir /path/to/output/dir -allele_refinement all -kraken_path /path/to/kraken/ -if -l sample.HLAProfiler.log
```

For more details on usage run:
```
perl HLAProfiler.pl predict -h
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
The database is included in the initial GitHub release.

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
