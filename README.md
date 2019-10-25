[![HitCount](http://hits.dwyl.io/adigenova/wengan.svg)](http://hits.dwyl.io/adigenova/wengan)

# Wengan
An accurate and ultra-fast genome assembler

# SYNOPSIS

    # Assembling Oxford nanopore and illumina reads with WenganM
     wengan.pl -x ontraw -a M -s lib1.fwd.fastq.gz,lib1.rev.fastq.gz -l ont.fastq.gz -p asm1 -t 20 -g 3000

    # Assembling PacBio reads and illumina reads with WenganA
     wengan.pl -x pacraw -a A -s lib1.fwd.fastq.gz,lib1.rev.fastq.gz -l pac.fastq.gz -p asm2 -t 20 -g 3000

    # Assembling ultra-long nanopore reads and BGI reads with WenganM
     wengan.pl -x ontlon -a M -s lib2.fwd.fastq.gz,lib2.rev.fastq.gz -l ont.fastq.gz -p asm3 -t 20 -g 3000

    # Non-hybrid assembly of PacBio Circular Consensus Sequence data with WenganM
     wengan.pl -x pacccs -a M -l ccs.fastq.gz -p asm4 -t 20 -g 3000

    # Assembling ultra-long nanopore reads and Illumina reads with WenganD (need a high memory machine 600Gb)
     wengan.pl -x ontlon -a D -s lib2.fwd.fastq.gz,lib2.rev.fastq.gz -l ont.fastq.gz -p asm5 -t 20 -g 3000

    # Assembling pacraw reads wiht pre-assembled short-read contigs from Minia3
     wengan.pl -x pacraw -a M -s lib1.fwd.fastq.gz,lib1.rev.fastq.gz -l pac.fastq.gz -p asm6 -t 20 -g 3000 -c contigs.minia.fa

    # Assembling pacraw reads wiht pre-assembled short-read contigs from Abyss
     wengan.pl -x pacraw -a A -s lib1.fwd.fastq.gz,lib1.rev.fastq.gz -l pac.fastq.gz -p asm7 -t 20 -g 3000 -c contigs.abyss.fa

    # Assembling pacraw reads wiht pre-assembled short-read contigs from DiscovarDenovo
     wengan.pl -x pacraw -a D -s lib1.fwd.fastq.gz,lib1.rev.fastq.gz -l pac.fastq.gz -p asm8 -t 20 -g 3000 -c contigs.disco.fa

# DESCRIPTION

**Wengan** is a new genome assembler that unlike most of the current long-reads assemblers avoid entirely the all-vs-all read comparison.
The key idea behind **Wengan** is that long-read alignments can be **inferred by building paths** on a sequence graph. To achieve this, **Wengan** build a new sequence graph called the Synthetic Scaffolding Graph. The SSG is build from a spectrum of synthetic mate-pair libraries extracted from raw long-reads. Then, longer alignments are build by peforming a transitive reduction of the edges.
Another distinct feature of **Wengan** is that perform **self-validation** by following the read information. **Wengan** identify miss-assemblies at differents steps of the assembly process. For more information about the algorithmic ideas behind **Wengan** please read the preprint available on bioRxiv.

# SHORT-READ ASSEMBLY

**Wengan** uses a de bruijn graph assembler to build the assembly backbone from short-read data.
Currently, **Wengan** can use **Minia3**, **Abyss2** or **DiscoVarDenovo**.  The recomended short-read coverage
is **50-60X** of 150bp x 2 or 250bp x 2 short reads.

## WenganM \[M\]

This **Wengan** mode use the **Minia3** short-read assembler, this is the fastest mode of **Wengan** and can assemble a complete human genome
in less than 210 CPU hours (~50Gb of RAM).

## WenganA \[A\]

This **Wengan** mode use the **Abyss2** short-read assembler, this is the lowest memory mode of **Wengan** and can assemble a complete human genome
in less than 40Gb of RAM (~900 CPU hours). This assembly mode takes  ~2 days when using 20 CPUs on a single machine.

## WenganD \[D\]

This **Wengan** mode use the **DiscovarDenovo** short-read assembler, this is the greedier memory mode of **Wengan** and for assembling a complete human genome need about 600Gb of RAM (~900 CPU hours).
This assembly mode takes ~2 days when using 20 CPUs on a single machine.

# LONGREADS PRESETS

The presets define several variables of the wengan pipeline execution and depends on the long-read technology used to sequence the genome.
The recommended long-read coverage is 30X.

## ontlon

preset for raw ultra-long-reads from Oxford Nanopore, typically with an  N50 > 50kb.

## ontraw

preset for raw long-reads Nanopore reads typically with an  N50 ~\[15kb-40kb\].

## pacraw

preset for raw long-reads from Pacific Bioscience (PacBio) typically with an  N50 ~\[8kb-60kb\].

## pacccs (experimental)

preset for Circular Consensus Sequences from Pacific Bioscience (PacBio) typically with an N50 ~\[15kb\]. This type of data is not fully supported yet.

# Wengan benchmark

| Genome | Long reads| Short reads|Wengan Mode| NG50 (Mb) | CPU (h) | RAM (Gb) | Fasta file|
|:------:|:------:|:------:|:------:|:------:|:------:|:------:|:------:|
|  |  | 2x150bp 50X (GIAB:[rs1][g1] , [rs2][g2])| WenganA      | 23.08   | 671     | 45  | [asm][NA12878.WenganA.ONT-ul-rel5.fa.gz]
| NA12878| ONT 35X ([rel5][rel5]) | 2x150bp 50X (GIAB:[rs1][g1] , [rs2][g2])| WenganM     | 16.67   | 185     | 53  | [asm][NA12878.WenganM.ONT-ul-rel5.fa.gz]
|  |  |  2x250bp 60X (ENA:[rs1][wdsna1] , [rs2][wdsna2])| WenganD      | 33.13   |  550    | 622  | [asm][NA12878.WenganD.ONT-ul-rel5.fa.gz]
| HG00073   | PAC 90X (ENA:[rl1][ena])| 2x250bp 63X (ENA:[rs1][wdhg1] , [rs2][wdhg2])| WenganD     | 29.2   | 800     | 644  | [asm][HG00733.WenganD.PAC-SequelI.fa.gz]
| NA24385   | ONT 60X (GIAB:[rl1][giab]) | 2x250bp 70X (GIAB:[rs1][g3])| WenganD     | 48.8   | 910     | 650  | [asm][NA24385.WenganD.ONT-ul-final.fa.gz]
| CHM13   | ONT 50X (T2T:[rel2][t2t])| 2x250bp 66X (ENA:[rs1][wdch1] , [rs2][wdch2])| WenganD     | 57.4   | 1027     | 647  | [asm][CHM13.WenganD.ONT-T2T-rel2.fa.gz]

[rel5]: https://github.com/nanopore-wgs-consortium/NA12878/blob/master/nanopore-human-genome/rel5.md
[t2t]: https://github.com/nanopore-wgs-consortium/CHM13
[ena]: https://www.ebi.ac.uk/ena/data/view/SRX4480530
[giab]: https://cutt.ly/BekQW1l

[wdsna1]: https://www.ebi.ac.uk/ena/data/view/SRR891258
[wdsna2]: https://www.ebi.ac.uk/ena/data/view/SRR891259
[wdch1]: https://www.ebi.ac.uk/ena/data/view/SRR3189742
[wdch2]: https://www.ebi.ac.uk/ena/data/view/SRR3189741
[wdhg1]: https://www.ebi.ac.uk/ena/data/view/SRR5534476
[wdhg2]: https://www.ebi.ac.uk/ena/data/view/SRR5534475

[g1]: https://cutt.ly/iekQmOn
[g2]: https://cutt.ly/EekQWrz
[g3]: https://cutt.ly/lekQWmZ

[NA12878.WenganA.ONT-ul-rel5.fa.gz]: https://zenodo.org/record/12598666/files/NA12878.WenganA.ONT-ul-rel5.fa.gz?download=1
[NA12878.WenganM.ONT-ul-rel5.fa.gz]: https://zenodo.org/record/12598666/files/NA12878.WenganM.ONT-ul-rel5.fa.gz?download=1
[NA12878.WenganD.ONT-ul-rel5.fa.gz]: https://zenodo.org/record/12598666/files/NA12878.WenganD.ONT-ul-rel5.fa.gz?download=1
[CHM13.WenganD.ONT-T2T-rel2.fa.gz]: https://zenodo.org/record/12598666/files/CHM13.WenganD.ONT-T2T-rel2.fa.gz?download=1
[NA24385.WenganD.ONT-ul-final.fa.gz]: https://zenodo.org/record/12598666/files/NA24385.WenganD.ONT-ul-final.fa.gz?download=1
[HG00733.WenganD.PAC-SequelI.fa.gz]: https://zenodo.org/record/12598666/files/HG00733.WenganD.PAC-SequelI.fa.gz?download=1

The assemblies generated using Wengan can be downloaded from [Zenodo](https://zenodo.org/record/12598666).
All the assemblies were ran as described in the wengan preprint. NG50 was caculated using a genome size of 3.14Gb.

# Wengan components
+ A bruijn graph assembler ([Minia](https://github.com/GATB/minia), [Abyss](https://github.com/bcgsc/abyss) or [DiscovarDenovo](https://software.broadinstitute.org/software/discovar/blog/))
+ [FastMIN-SG](https://github.com/adigenova/fastmin-sg)
+ [IntervalMiss](https://github.com/adigenova/intervalmiss)
+ [Liger](https://github.com/adigenova/liger)

<img src="./wengan-diagram.svg">

# About the name
**Wengan** is a [Mapudungun](https://en.wikipedia.org/wiki/Mapuche_language) word. The Mapudungun is the language of the [**Mapuche**](https://en.wikipedia.org/wiki/Mapuche) people, the largest indigenous inhabitants of south-central Chile. **Wengan** means "***Making the path***".

