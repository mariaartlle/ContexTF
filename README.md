# ContexTF: functional annotation of prokaryotic regulatory elements

### Table of Content
- [Biological framework](#1-biological-framework)
- [Pipeline overview](#2-pipeline-overview)
- [Installation and execution](#3-installation-and-execution)
- [Tutorial](#4-tutorial)
- [Bibliography](#4-Bibliography)

## 1. Biological framework
ContexTF uses the manually curated transcription factor database from model organisms to perform a homology-based search for putative transcription factors in unannotated prokaryotic genomes. It combines the results of three different homology searches: sequence, structure and genomic context. It is widely assumed that proteins with high sequence and, especially, structural similarity have high probability of performing the  same function (Krissinel, 2007). For this reason, most bioinformatic tools take advantage of these characteristic to annotate homologous proteins throughout species. In prokaryotic genomes, where most TFs are found in context with the proteins that they regulate, the conservation of the genes near the regulatory protein could also be used to annotate proteins in prokaryotic organisms. ContexTF integrates a classical search of homologous candidates and provides a snapshot of the genomic contexts of the involved proteins to evaluate the conservation and similarity of the flanking genes of the candidates.

## 2. Pipeline overview
![ContexTF overview](/ContexTF_overview.svg "ContexTF overview").
ContexTF is a standalone python pipeline that consists of three executable modules (Fig. 1). The first module relies on HMMER for performing a PFAM domain search of the model TFs present in the database into an unannotated genome. In this step, the PFAM architecture of each reference TF is established, removing the overlapping domains with lowest E-value. Next, these domain architectures are searched in the provided organism predicted proteome, ensuring that the TF candidates found share the same domain composition.

The second module checks the structural similarity between the reference TF and its candidates, using TM-align. This module execution is optional, as the Uniprot identifiers need to be mapped from the protein identifiers in the provided genome to download its AlphaFold predicted structure. If indicated, the HMMER results are used to superimpose the reference protein with each one of its candidates and the TM-score (Zhang & Skolnick, 2004), normalised by the reference protein’s length, is assessed. The TM-score threshold can be adjusted by the user; the default is set to consider proteins structurally similar when the TM-score is higher than 0.5, as it generally indicates that both proteins share the same fold in SCOP/CATH (Andreeva et al., 2020). As the structures are predicted, there might be unordered regions that interfere with the protein superimpositions. To avoid the loss of potential candidates with predicted regions with low confidence (pLDDT < 70) (Mariani et al., 2013), unordered re-
gions that represent more than 10% of the protein are removed before superimposing.

The final module of ContexTF yields the genomic context comparisons between the reference TF and the candidates obtained by HMMER or TM-align. All the reference TFs have been mapped to a genome or assembly from which the target protein will be plotted, along with the designated number of upstream and downstream coding sequences. This allows the user to compare a snapshot of the regions and assess the similarity between the reference and candidates. To facilitate the visual comparison, the flanking proteins are coloured by protein families found by an all-against-all search by MMSeqs2 or PSIBLAST, depending on the user’s preferences. The final output of the pipeline is a collection of figures, one for each reference TF and its corresponding candidates.

## 3. Installation

ContexTF is a standalone pipeline developed in Python 3.10+. To install ContexTF in your local computer execute the following commands: 

```shell
git clone https://github.com/mariaartlle/ContexTF.git
cd ContexTF
pip3 install .
```

To ensure the portability of ContexTF, a [Conda](https://docs.conda.io/en/latest/) environment is provided with all the python packages needed to execute it and their correspondant versions. To install, execute:
```shell
conda env create -n ContexTF_env --file ContexTF_env.yml
conda activate ContexTF_env
```

ContexTF relies on the local installation of three bioinformatic softwares:

- [HMMER](http://hmmer.org/) (version 3.3.2)
- [TMAlign](https://zhanggroup.org/TM-align/) (version 20190822)
- [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) (version 13-45111+ds-2) or [MMseqs2](https://github.com/soedinglab/MMseqs2) (version 2.12.0+)

## 4. Tutorial

### Basic execution example
The only required arguments to execute the pipeline are a [genomic FASTA file](/example/input/GCA_000008345.1_ASM834v1_genomic.fna) and a [GFF3 annotation file](/example/input/GCA_000008345.1_ASM834v1_genomic.gff): 
```shell 
ContexTF -fna example/input/GCA_000008345.1_ASM834v1_genomic.fna -gff3 example/input/GCA_000008345.1_ASM834v1_genomic.gff
```
With this we will obtain a tabular file summary of the HMMER candidates found in the provided genome and a genomic context figure for each reference TF found.

### More advanced example
ContexTF has more arguments that allow the user to refine the output. In this example, we are using the same genomic files but this time, the structural homology module will be executed and only the candidates with a TM-score higher than 0.7 will be considered as candidates. Regarding the genomic context figures, the protein families will be found with MMseqs2 and for each TF 8 genes will be plotted in the 5' site, and 4 in the 3'.

```shell
ContexTF -fna example/input/GCA_000008345.1_ASM834v1_genomic.fna -gff3 example/input/GCA_000008345.1_ASM834v1_genomic.gff --tmalign --tmalign_min_threshold 0.7 --GC_method mmseqs --n_flanking5 8 --n_flanking3 4 

```

## Bibliography
Andreeva, A., Kulesha, E., Gough, J., & Murzin, A. G. (2020). The SCOP database in 2020: expanded classification of representative family and superfamily domains of known protein structures. Nucleic Acids Research, 48(D1), D376–D382. https://doi.org/10.1093/nar/gkz1064

Krissinel, E. (2007). On the relationship between sequence and structure similarities in proteomics. Bioinformatics, 23(6), 717–723. https://doi.org/10.1093/bioinformatics/btm006

Mariani, V., Biasini, M., Barbato, A., & Schwede, T. (2013). IDDT: A local superposition-free score for comparing protein structures and models using distance difference tests. Bioinformatics, 29(21), 2722–2728. https://doi.org/10.1093/bioinformatics/btt473

Zhang, Y., & Skolnick, J. (2005). TM-align: A protein structure alignment algorithm based on the TM-score. Nucleic Acids Research, 33(7), 2302–2309. https://doi.org/10.1093/nar/gki524

