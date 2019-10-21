# Tourette Syndrome Candidate Gene Analysis (TSCGA)
This repository holds the script that gathers information on mutations and genes from multiple database files and writes to an output file.  

## Types of input files     
### Mutation file  
This is a tab delimited file that would contain any of the following information:  
* Mutation position  
* Type of mutation: ins, del, snp, or mnp
* Gene that the mutation is found in  
* Segregation pattern of the mutation if running pedigree analysis  

Multiple mutation files can be run at once by specifying them on the command line explicitly or running a folder of mutation files with the `folder/path/*` notation on the command line.  
Mutation files are named by the family and the inheritance pattern if pedigree information is being run. Ex. Family1_Dominant.txt

### Database Files  
#### Reflist files  
This input if for running pedigree analysis.  
Reflist files are the family segregation information and is used to compare the family segregation with the mutation segregation.  
Multiple families can be run at once by reading in the reflist files in the same manner as the mutation files.  
The program reads the file names and matches the segregations off of that so the family names in the file name must match for this reason. 

#### Annotation file  
This input is in the form of a annotated vcf file. The mutation position from the mutation file is matched to the corresponding line in the annotation file and the information is taken from there.  
The program takes up to the follwoing information from the annotated file:  
* Sift score  
* Polyphen-2 HVAR and HDIV scores  
* ExAC allele frequency  
* gnomAD allele frequency  
* Clinvar information  
* Interpro_domain information

#### Gene tolerance score file  
This input contains the following information in a tab delimited file:
* Gene symbol to match with gene symbol in the mutation file  
* mis_z score    
* pLI score  

#### Disease association file and diease key file  
The disease association file input has all the diseases that have been linked to a gene based on publication text mining.
The format of the lines are: Gene|disease_1;disease_2;disease_3 etc.   
The file has 1 gene per line  

The disease key file is a text file that has one element per line. The elements are divided by the headings
* general:  
  * element_1
  * element_2
* specific:  
  * element_1
* exclude:  
  * element_1

The heading tells the program whether to look for diseases with that keyword case insensitive (general), case sensitive (specific), or to exclude diseases with that keyword (exclude).  
Any combination of these categories are allowed and if this file is not supplied at all then all dieases are printed to the output.

#### GTEx tpm file  
This file has tpm numbers for all genes in over 40 tissues, collected by the Genome Tissue Expression Consortium project.
The file is a .gct file and holds the median tpm for all genes in each of the tissues.  
The program collects tissue expression levels and outputs the following:  
* Top tissue expression for:
  * Brain tissue
  * Non-Brain tissue  
  * All tissues
* Median tissue expression for:
  * Brain tissue
  * Non-Brain tissue  
  * All tissues
  
This allows for a quick look at if the gene is expressed higher in the brain than in other tissues  

#### Other Gene information files  
These files are user generated to add additional evidence for the genes.  
These files are one gene name per line and the program reports whether the gene in the mutation file was found in this given file  
Ex. For the TS project there were genes found oin a GWAS study and we wanted to add that information into the mutation files. We put the GWAS found genes into a file and the program reported in the mutation file if that gene was found in this GWAS study.  

Any number of these gene files can be supplied for additional information of genes in the mutation file.  
The program reads in the file and a string given, the string is used as a column header to tell if the gene was found in the provided file.  


## Output files  
The output file is the mutation file with the database information added to the end of the row.  
The program adds column headers to the end of the header row and fills in the information for the mutations in those columns.  
The program also adds a tag to the end of the file name and before the .txt extension to let the user know what type of information was added to the file.  
