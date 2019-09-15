# Introduction {#intro}

## Types of microbiota data

Next-generation sequencing has greatly increased our ability to survey the gut microbiota over the past 10-15 years. The most widely used technologies to assess microbial community data are 16S rRNA amplicon sequencing and shotgun metagenomics sequencing.  Both have been globaly referred to as "metagenomics", but we will make sure to differentiate the technologies.

### 16S rRNA amplicon sequencing

Sequences a single gene, the 16S ribosomal RNA (rRNA) gene.  The 16S rRNA gene is ubiquitous in the bacterial kingdom, but is highly variable among different types of bacteria.  About 80% of bacterial RNA is made up of rRNA, so this technology can identify rare bacteria groups at a high sensitivity ([Li, 2015](https://www.annualreviews.org/doi/pdf/10.1146/annurev-statistics-010814-020351)).

The advantage of 16S rRNA sequencing is that it is fairly inexpensive, routinely done, and has good sensitivey at the community level. However, there are 2 major limitations associated with this technology: 

1. It is not sensitive at species level, and 

2. It does not provide any information regarding genetic potential or bacterial functions.

16S rRNA sequencing sometimes gets a bad wrap because of these limitations, but there good reasons to use this technology.  Again, its much less expensive than shotgun metagenomics and is still very sensitive at what it can read. It is a good option if you have many samples, a limited budget, and if you just want to identify whether there is an alteration in the composition of the microbiota.  

### Shotgun metagenomic sequencing

Sequences all microbial genomic DNA.  Provides species level compositional data and bacterial genetic data.  The compositional data is exactly similar to 16S rRNA data, but superior because it gives species level depth.  However, it is much more expensive on a per sample basis and requires significantly greater computational power to process.  
  
You can label chapter and section titles using `{#label}` after them, e.g., we can reference Chapter \@ref(intro). If you do not manually label them, there will be automatic labels anyway, e.g., Chapter \@ref(methods).

## Microbial composition data from sequencing technologies

The overall structure of 16S rRNA and shotgun metagenomic sequencing data are similar. Lets quickly describe how the data is aquired to understand why this is.  Note, this workshop is set up for statistical analysis, so we will keep this section at a very, very basic level.

For both technologies, bacterial DNA must be isolated from the sample (e.g., feces, intestinal content, skin, soil, water, dog hair, etc.).  Once the bacterial DNA is isolated, a subsample of the DNA is extracted for PCR, a subset of the resulting amplicon is pooled into a library, and a subset of the library is then sequenced ([Morton et al., 2019](https://www.nature.com/articles/s41467-019-10656-5)). Main difference between the two technologies is the libraries.  Again, 16S rRNA libraries are focused on the 16S rRNA gene, while shotgun metagenomics is focused on all bacterial genetic material.  

In both cases, you end up with raw sequences that have to be mapped to referenced phylogenetic tree.  Generally, the sequences are aligned to reference gene until a set of sequences fill a percentage of the gene at a certain threshold.  Once this threshold or probability is reached, then these clusters of sequences are assigned to the phylogeny.  This is now considered a single "read" in your data.  The accumulation of these read counts provide us the data that can then be statistically analyzed.

There are some considerations you should consider when working with next-generation microbial sequencing data.  

1. It is essential to reach a minimum sampling depth.

2. It is semi-quantitative data, i.e., counts of sequencing reads are not actually measuring a single bacteria.  Data is commonly referred to as "relative abundance."

3. It does not take into consideration the microbial load of each sample.  

## Statistical analysis of microbial composition data

Microbial composition data is multi-level and high-dimensional. What do I mean by these terms?  It is multi-level because it can be analyzed at various taxonomy levels, e.g., species or phylum level.  It is high-dimensional because you will likely have hundreds to tens of thousands of variables (i.e., bacterial taxa) to analyze.  Spreadsheet or GUI (graphical user interface) based statistical software may be cumbersome or difficult to work with for these reasons.  Therefore, it is highly recommended to use language based stastical software, like R, that can be tailored to meet the needs of these factors. 

R is a freely available statistical language and environment that is built around an extremely large and active user group.  In fact, this web-book was created using R.  While there are many, many benefits to using R, those without progamming experience may have a steep learning curve.  The next chapter will be dedicated to detailing the basic knowledge required to analyze microbial sequencing data in R.

