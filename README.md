# PLINK

PLINK is a useful tool to quickly process and manipulate GWAS data formats. For more information, visit a [developper's website](http://zzz.bwh.harvard.edu/plink/index.shtml)

## PLINK data format
PLINK has its own format for genotype encoding, which is widely adopted by many GWAS data sources. 

* PED/MAP (non-binary format): .ped is basically a space-delimited text file that writes genotypes for each person in one row. There are 6 columns that precedes genotypes, which are reserved for family ID/individual ID/paternal ID/maternal ID/sex/phenotype(usually indicating case/control status) in this order.  

.ped file has to accompany .map file which contains genomic locations (chromosome number/base-pair coordinates in the chromosome) for the SNPs contained in the genotype array.

* BED (binary format): .bed is a compressed version of .ped and is recommended because it saves space and time. When you are received a .ped file, it can be converted to .bed file using the --make-bed option in PLINK::

```
./plink --noweb --file toy --make-bed --out test_convert
```

Unlike .ped, which contains data pertaining family relations and phenotypes, .bed purely contains genotypes. Thus, when .ped is converted to .bed, the family data is separately written in a .fam file. The binary format accompanies a genomic map file .bim file which is equivalent to the .map file in non-binary file system.  

This shows what's in the .bim file from the Hapmap phase 2 dataset:

```
head -n 10 hapmap_JPT_CHB_r23a.bim
```

More information on file formats can be found [here](https://www.cog-genomics.org/plink2/formats#map) 

## Filtering
The following command will extract only chromosome 8 data from the hapmap data:

```
./plink --noweb --bfile hapmap_JPT_CHB_r23a --chr 8 --make-bed --out hapmap_JPTCHB_chr8 
```

Explanation on the command above:
* --bfile loads a binary .bed file as an input. A file name without extension should be provided.
* --chr 8 selects only chromosome 8
* --make-bed --out creates a new .bed file from the chromosome 8 data loaded from the original .bed file

Every command in PLINK will automatically generate a .log file with a file name as provided by --out option.

For more information on extracting specific chromosomes/SNP ids, see [here] (http://zzz.bwh.harvard.edu/plink/dataman.shtml)

## Summary statistics on SNPs, and quality control (QC) 
The first step in a typical GWAS is to remove errorneous genotypes or individuals that can lead to false association results. There are 3 widely used criteria:

* Remove genotypes with a frequency of a minor allele in the array, or minor allele frequency (MAF) less than 5% (of 1%) (also known as "rare alleles") 
* Remove individuals with genotype calling rate less than 95%
* Remove genotypes with calling rate less than 95% 
* Remove SNPs that deviates from Hardy-Weinberg Equilibrium (HWE). Under the assumption that the individuals in the dataset are unrelated, deviation from HWE can indicate genotype error (HWE is broken when there is strong evolutionary selection/inbreeding)

Consult [this reference](https://www.ncbi.nlm.nih.gov/pubmed/21085122) for 
more details.

* Allele frequency
```
./plink --noweb --bfile hapmap_JPTCHB_chr8 --freq --out hapmap_chr8_freq
```
The allele frequency statistics per SNP are now stored in hapmap_chr8_freq.frq.

* Genotype calling(missing) rate
```
./plink --noweb --bfile hapmap_JPTCHB_chr8 --missing --out hapmap_chr8_missing
```
Per-individual and per-SNP missing rates are stored in .imiss and .lmiss, respectively.

* HWE
```
./plink --noweb --bfile hapmap_JPTCHB_chr8 --hardy --out hapmap_chr8_hardy
```
HWE p-values for each SNP are stored in .hwe file (lower p-value indicates deviation from the equilibrium).

* Quality control of the genotype dataset
The command below creates a new genotype dataset from the original data by excluding the individuals/SNPs that fails the quality cutoffs as given in command options:

```
./plink --noweb --bfile hapmap_JPTCHB_chr8 --mind 0.1 --geno 0.05 --maf 0.01 --hwe 0.00001 --make-bed --out hapmap_JPTCHB_chr8_QC 
```

In this example, the maximum per-individual missing rate 10%, maximum per-SNP missing rate 5%, minimum MAF of 1%, and minimum HWE p-value of 1E-5 are applied as thresholds.

## Linkage disequilibrium (LD) and pruning
Linkage disequilibrium refers to non-random association between two loci. LD occurs because the SNPs close to each other tend to be inherited together as a unit (haplotype). LD has several implications to keep in mind in GWAS: For example, a SNP that was discovered from one GWAS (lead SNP) is not necessarily causal to the phenotype of interest: it may be that the SNP is simply a proxy SNP which is in LD with a SNP that truly causes the effect (causal SNP). Some studies remove highly correlated SNPs in order to reduce the number of association tests so that less stringent multiple testing correction can be applied. 

Using plink, the genotype data can be 'pruned' by removing SNPs that are in LD above certain level of correlation:

```
./plink --noweb --bfile hapmap_JPTCHB_chr8_QC --indep-pairwise 50 5 0.2 --out LDsnps
```
The command above only creates a list of the SNPs that are in linkage equilibrium (high statistical independence). In order to actually create a pruned dataset, use --extract option:

```
./plink --noweb --bfile hapmap_JPTCHB_chr8_QC --extract LDsnps.prune.in --make-bed --out hapmap_JPTCHB_chr8_QC_pruned 
```
# GenABEL

Genotype data in plink format can be imported to R for custom analyses. GenABEL is an R package that facilitates GWAS workflow by fast data import/export/manipulation and many useful functions.

GenABEL can be installed and loaded in a R session:

```
install.packages("GenABEL")
```
## Importing PLINK to R using GenABEL 

GenABEL cannot recognize a binary .bed format, so the data first has to be converted to .ped using plink first (for the sake of running speed, we will use the pruned chromosome 8 data from the previous example)

```
./plink --noweb --bfile hapmap_JPTCHB_chr8_QC_pruned --recode --tab --out hapmap_JPTCHB_chr8_QC_pruned
``` 
 
This will create hapmap_JPTCHB_chr8_QC_pruned.ped/map.

Now go back to the R session for an extra step of converting the data to .raw format, which is a final form of genotype data that can be recognized by GenABEL:

```
library(GenABEL)
convert.snp.ped('path_to_ped_file','path_to_map_file','hapmap_JPTCHB_chr8_QC_pruned.raw')
```
Now import the converted genotype file, along with the phenotype file (containing any information other than genotypes) into the R session by:

```
dt <- load.gwaa.data('path_to_raw_file','path_to_phenotype_file')
```

The phenotype file is included in the directory /ptdata/ under the github clone.

 
















 


