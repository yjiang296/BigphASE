# BigphASE

Bi-parental graph ASE analyzer
<br>
This python pakage was devised to analyze hybrid RNA-seq using bi-parental graph strategy.
<br>
version: 0.1.2
<br>
author: Yi Jiang (<jiangy296@mail2.sysu.edu.cn>); Muran Yi; Zhule Liu

<h3>Enviroment required:</h3>
xxx

<h3>Installation:</h3>
`pip install BigphASE`

<h3>Usage:</h3>
BigphASEtools [-options] <parent_A genome sequence fa> <parent_B genome sequence fa> <parent_A genome gff3> <clean RNAseq R1 fastq> <clean RNAseq R2 fastq> <prefix>

<h3>Main steps:</h3>

(1) Graph builder
Graph genome was built using one parental genome as backbone, with variations between two parental genomes as alternate paths.
(2) RNAseq mapping
Clean RNA-seq data were mapped to graph index.
(3) ASE genes analysis
After mapping, reads were assigned to parental origin based on their genotypes at SNP sites, then maternal and paternal origin reads count of each gene were summed.

![pipeline](https://github.com/yjiang296/BigphASE/assets/115337217/e85cd6dc-2704-4e55-a3f2-aaab0b281fce)
