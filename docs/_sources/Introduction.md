# IronMatt (Matt Larson Foundation) Single Cell Sequencing Analysis for Undergraduates

## Performed in R using the Seurat package

## Understanding RNA 

### Structure

Ribonucleic Acid, RNA, is a single strand of nucleotides that is composed of four nitrogenous bases: adenine (A), guanine (G), cytosine (C), and uracil (U). The bases can interact with eachother through hydrogen bonding and enable RNA to form complex structures. 

The central dogma explains the flow of genetic information. It starts with deoxyribonucleic acid (DNA); DNA is then transcribed into a specific RNA. Once the specific RNA is transcribed, RNA can either be functional on its own or proceed to be translated into a protein. One common misconception is that all DNA encodes for specific proteins. While some RNA converted from DNA is used to encode a specific protein, not all RNAs are translated into proteins. A very small percentage of the human genome actually contains protein coding DNA while the rest of the genome is non-coding DNA (ncDNA). DNA and RNA’s structure affect the functional use of both nucleic acids. DNA’s double helix structure is essential to maintain stability to ensure the preservation of genetic code. RNA’s single stranded structure allows RNA to be more dynamic, including its ability to be folded into specific structures.  

### Functions

Once RNA is formed, one of its functions is to carry genetic information from DNA to ribosomes for protein synthesis. Throughout this process, there are many different types of RNA that work together to achieve this goal, primarily:  

- **Messenger RNA (mRNA)** carries genetic information from DNA to ribosomes where protein synthesis occurs 

- **Ribosomal RNA (rRNA)** helps to form the structure of ribosomes and catalyzes the formation of peptide bonds to connect the amino acids to create the proteins. 


RNA molecules that do not assist in protein synthesis can perform other functions such as gene regulation, catalyze reactions, or provide structural support. 

- **Long Non-Coding RNA (lncRNA)** help in gene expression/regulation (cell responses) 

 
RNA is an important area of research because it can be utilized to study gene expression, diagnose diseases and therapy responses based on RNA levels within a sample, or develop therapies to target RNA or other molecules within a cell.  

![image.png](https://assets.technologynetworks.com/production/dynamic/images/content/296719/what-are-the-key-differences-between-dna-and-rna-296719-1280x720.webp?cb=13322129)
Image Source Technology Netwrks


 

# RNA Sequencing

## RNA Sequencing  

Sequencing is the process of determining the order of nucleotides within a molecule. Although this may sound simple and straightforward, there is a multistep process that is fairly complex. Understanding the process of sequencing is crucial for accurately analyzing genetic information, identifying mutations, and advancing research in medicine and biotechnology.  

Definition from NHGRI 

“RNA-seq is a method for sequencing an entire set of RNA molecules. This method involves isolating RNA molecules in the tissue or cell sample, copying those RNAs into DNA, and then sequencing the resulting DNA molecules. RNA-seq helps researchers determine which segments of DNA have been transcribed into RNA and how much a gene is expressed.” 



# Single Cell RNA Sequencing Platforms

# 10X Genomics 

Single cell transcriptomics examines the gene expression of individual cells from the same sample. Each cell from the sample is isolated and combined with the necessary reagents to separate the RNA to convert it into cDNA within a droplet. 

![Description](https://cdn.10xgenomics.com/image/upload/f_auto,q_auto,w_900,dpr_2.0/v1709753672/blog/GEM-X%20Launch%20blog/vertical_flipped_GEM96_H1975_high_4_S0001.gif)
![image.png](https://cdn.10xgenomics.com/image/upload/v1709930681/blog/GEM-X%20Launch%20blog/Figure_1.png)
Image Source 10X Genomics

After various steps, including amplification (to make more copies for better resolution), tagging, and barcoding (a unique code that tells you which droplet it came from since you will eventually combine the data), the cDNA is sequenced in a standard sequencing platform such as an Illumina Sequencer and processed through computational pipelines. 

![image.png](https://cdn.10xgenomics.com/image/upload/v1710206114/blog/GEM-X%20Launch%20blog/Figure_5.png)
Image Source 10X Genomics

The resulting data can be used by researchers to visualize the genes expressed within various cells in a sample, and provides an in-depth exploration of cellular functions and their intricate interactions within complex tissues.

The power of single cell sequencing has transformed science and medicine. For example, a researcher using bulk RNA-seq to study immune responses in an inflamed tissue might miss a small subset of hyperactivated immune cells driving excessive inflammation. This could lead to an incomplete understanding of disease mechanisms and hinder the identification of precise therapeutic targets.



This tutorial was created by the WiseLab of Salve Regina Univeristy with the support of the Matt Larson Foundation
https://ironmatt.org/ and RI-INBRE 

