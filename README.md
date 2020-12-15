[![HitCount](http://hits.dwyl.io/adigenova/wengan.svg)](http://hits.dwyl.io/adigenova/wengan)


# Wengan
An accurate and ultra-fast genome assembler
### Version: 0.2 (18/05/2020)

Table of Contents
=================

   * [SYNOPSIS](#synopsis)
   * [Description](#description)
   * [Short-read assembly](#short-read-assembly)
      * [WenganM [M]](#wenganm-m)
      * [WenganA [A]](#wengana-a)
      * [WenganD [D]](#wengand-d)
   * [Long-read presets](#long-read-presets)
      * [ontlon](#ontlon)
      * [ontraw](#ontraw)
      * [pacraw](#pacraw)
      * [pacccs (experimental)](#pacccs-experimental)
   * [Wengan demo](#wengan-demo)
   * [Wengan benchmark](#wengan-benchmark)
   * [Wengan components](#wengan-components)
   * [Getting the latest source code](#getting-the-latest-source-code)
      * [Instructions](#instructions)
         * [Containers](#containers)
         * [Building Wengan from source](#building-wengan-from-source)
            * [Requirements](#requirements)
   * [Limitations](#limitations)
   * [About the name](#about-the-name)
   * [Citation](#citation)

# SYNOPSIS

    # Assembling Oxford Nanopore and Illumina reads with WenganM
     wengan.pl -x ontraw -a M -s lib1.fwd.fastq.gz,lib1.rev.fastq.gz -l ont.fastq.gz -p asm1 -t 20 -g 3000

    # Assembling PacBio reads and Illumina reads with WenganA
     wengan.pl -x pacraw -a A -s lib1.fwd.fastq.gz,lib1.rev.fastq.gz -l pac.fastq.gz -p asm2 -t 20 -g 3000

    # Assembling ultra-long Nanopore reads and BGI reads with WenganM
     wengan.pl -x ontlon -a M -s lib2.fwd.fastq.gz,lib2.rev.fastq.gz -l ont.fastq.gz -p asm3 -t 20 -g 3000

    # Hybrid long-read only assembly of PacBio Circular Consensus Sequence and Nanopore data with WenganM
     wengan.pl -x ccsont -a M -l ont.fastq.gz -b ccs.fastq.gz -p asm4 -t 20 -g 3000
 
    # Assembling ultra-long Nanopore reads and Illumina reads with WenganD (need a high memory machine 600GB)
     wengan.pl -x ontlon -a D -s lib2.fwd.fastq.gz,lib2.rev.fastq.gz -l ont.fastq.gz -p asm5 -t 20 -g 3000

    # Assembling pacraw reads with pre-assembled short-read contigs from Minia3
     wengan.pl -x pacraw -a M -s lib1.fwd.fastq.gz,lib1.rev.fastq.gz -l pac.fastq.gz -p asm6 -t 20 -g 3000 -c contigs.minia.fa

    # Assembling pacraw reads with pre-assembled short-read contigs from Abyss
     wengan.pl -x pacraw -a A -s lib1.fwd.fastq.gz,lib1.rev.fastq.gz -l pac.fastq.gz -p asm7 -t 20 -g 3000 -c contigs.abyss.fa

    # Assembling pacraw reads with pre-assembled short-read contigs from DiscovarDenovo
     wengan.pl -x pacraw -a D -s lib1.fwd.fastq.gz,lib1.rev.fastq.gz -l pac.fastq.gz -p asm8 -t 20 -g 3000 -c contigs.disco.fa

# Description

**Wengan** is a new genome assembler that, unlike most of the current long-reads assemblers, avoids entirely the all-vs-all read comparison.
The key idea behind **Wengan** is that long-read alignments can be **inferred by building paths** on a sequence graph. To achieve this, **Wengan** builds a new sequence graph called the Synthetic Scaffolding Graph (SSG). The SSG is built from a spectrum of synthetic mate-pair libraries extracted from raw long-reads. Longer alignments are then built by performing a transitive reduction of the edges.
Another distinct feature of **Wengan** is that it performs **self-validation** by following the read information. **Wengan** identifies miss-assemblies at different steps of the assembly process. For more information about the algorithmic ideas behind **Wengan**, please read the preprint available in bioRxiv.


# Short-read assembly

**Wengan** uses a de Bruijn graph assembler to build the assembly backbone from short-read data.
Currently, **Wengan** can use **Minia3**, **Abyss2** or **DiscoVarDenovo**.  The recommended short-read coverage
is **50-60X** of 2 x 150bp or 2 x 250bp reads.

## WenganM \[M\]

This **Wengan** mode uses the **Minia3** short-read assembler. This is the fastest mode of **Wengan** and can assemble a complete human genome in less than 210 CPU hours (~50GB of RAM).

## WenganA \[A\]

This **Wengan** mode uses the **Abyss2** short-read assembler. This is the lowest memory mode of **Wengan** and can assemble a complete human genome with less than 40GB of RAM (~900 CPU hours). This assembly mode takes ~2 days when using 20 CPUs on a single machine.

## WenganD \[D\]

This **Wengan** mode uses the **DiscovarDenovo** short-read assembler. This is the greedier memory mode of **Wengan** and for assembling a complete human genome needs about 600GB of RAM (~900 CPU hours). This assembly mode takes ~2 days when using 20 CPUs on a single machine.

# Long-read presets

The presets define several variables of the Wengan pipeline execution and depend on the long-read technology used to sequence the genome.
The recommended long-read coverage is 30X.

## ontlon

preset for raw ultra-long-reads from Oxford Nanopore, typically with an N50 > 50kb.

## ontraw

preset for raw Nanopore reads typically with an N50 ~\[15kb-40kb\].

## pacraw

preset for raw long-reads from Pacific Bioscience (PacBio) typically with an N50 ~\[8kb-60kb\].

## pacccs (experimental)

preset for Circular Consensus Sequences from Pacific Bioscience (PacBio) typically with an N50 ~\[15kb\]. This type of data is not fully supported yet.

# Wengan demo
The repository [wengan_demo](https://github.com/adigenova/wengan_demo) contains a small dataset and instructions to test [Wengan v0.2](https://github.com/adigenova/wengan/releases/tag/v0.2).

```
#fetch the demo dataset
git clone https://github.com/adigenova/wengan_demo.git
```

# Wengan benchmark

| Genome | Long reads| Short reads|Wengan Mode| NG50 (Mb) | CPU (h) | RAM (GB) | Fasta file|
|:------:|:------:|:------:|:------:|:------:|:------:|:------:|:------:|
|  |  | 2x150bp 50X (GIAB:[rs1][g1] , [rs2][g2])| WenganA      | 25.99   | 725     | 45  | [asm][NA12878.WenganA.ONT-ul-rel5.fa.gz]
| NA12878| ONT 35X ([rel5][rel5]) | 2x150bp 50X (GIAB:[rs1][g1] , [rs2][g2])| WenganM     | 17.23   | 203     | 53  | [asm][NA12878.WenganM.ONT-ul-rel5.fa.gz]
|  |  |  2x250bp 60X (ENA:[rs1][wdsna1] , [rs2][wdsna2])| WenganD      | 35.31   |  589    | 622  | [asm][NA12878.WenganD.ONT-ul-rel5.fa.gz]
| HG00073   | PAC 90X (ENA:[rl1][ena])| 2x250bp 63X (ENA:[rs1][wdhg1] , [rs2][wdhg2])| WenganD     | 32.35   | 936   | 644  | [asm][HG00733.WenganD.PAC-SequelI.fa.gz]
| NA24385   | ONT 60X (GIAB:[rl1][giab]) | 2x250bp 70X (GIAB:[rs1][g3])| WenganD     | 50.59  | 963     | 651  | [asm][NA24385.WenganD.ONT-ul-final.fa.gz]
| CHM13   | ONT 50X (T2T:[rel3][t2t])| 2x250bp 66X (ENA:[rs1][wdch1] , [rs2][wdch2])| WenganD     | 69.72   | 1198     | 646  | [asm][CHM13.WenganD.ILL_UL_R3.fa.gz]

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

[NA12878.WenganA.ONT-ul-rel5.fa.gz]: https://zenodo.org/record/3779515/files/NA12878.WenganA.ONT-ul-rel5.fa.gz?download=1
[NA12878.WenganM.ONT-ul-rel5.fa.gz]: https://zenodo.org/record/3779515/files/NA12878.WenganM.ONT-ul-rel5.fa.gz?download=1
[NA12878.WenganD.ONT-ul-rel5.fa.gz]: https://zenodo.org/record/3779515/files/NA12878.WenganD.ONT-ul-rel5.fa.gz?download=1
[CHM13.WenganD.ILL_UL_R3.fa.gz]: https://zenodo.org/record/3779515/files/CHM13.WenganD.ILL_UL_R3.fa.gz?download=1
[NA24385.WenganD.ONT-ul-final.fa.gz]: https://zenodo.org/record/3779515/files/NA24385.WenganD.ONT-ul-final.fa.gz?download=1
[HG00733.WenganD.PAC-SequelI.fa.gz]: https://zenodo.org/record/3779515/files/HG00733.WenganD.PAC-SequelI.fa.gz?download=1

The assemblies generated using Wengan (v0.2) can be downloaded from [Zenodo](https://zenodo.org/record/3779515).
All the assemblies were ran as described in the Wengan manuscript. NG50 was computed using a genome size of 3.08Gb.

# Wengan components
+ A de Bruijn graph assembler ([Minia](https://github.com/GATB/minia), [Abyss](https://github.com/bcgsc/abyss) or [DiscovarDenovo](https://software.broadinstitute.org/software/discovar/blog/))
+ [FastMIN-SG](https://github.com/adigenova/fastmin-sg)
+ [IntervalMiss](https://github.com/adigenova/intervalmiss)
+ [Liger](https://github.com/adigenova/liger)

<img src="./wengan-diagram.svg">

# Getting the latest source code

## Instructions
It is recommended to use/download the latest binary release (Linux) from :
https://github.com/adigenova/wengan/releases

### Containers
To facilitate the execution of Wengan, we provide docker/singularity containers.
Wengan images are hosted on [Dockerhub](https://hub.docker.com/repository/docker/adigenova/wengan) and can be downloaded with the command:

```
docker pull adigenova/wengan:v0.2
```

Alternatively, using singularity:

```
export TMPDIR=/tmp
singularity pull docker://adigenova/wengan:v0.2
```

#### Run WenganM using singularity
```
#using singularity
CONTAINER=/path_to_container/wengan_v0.2.sif

#location of wengan in the container
WENGAN=/wengan/wengan-v0.2-bin-Linux/wengan.pl

#run WenganM with singularity exec
singularity exec $CONTAINER perl ${WENGAN} \
 -x pacraw \
 -a M \
 -s short.R1.fastq.gz,short.R2.fastq.gz \
 -l pacbio.clr.fastq.gz \
 -p asm_wengan -t 20 -g 3000
```


### Building Wengan from source
To compile Wengan run the following command:

```bash
#fetch Wengan and its components
git clone --recursive https://github.com/adigenova/wengan.git wengan
```

There are specific instructions for each Wengan component. 
After compilation you have to copy the binaries to wengan-dir/bin. 

#### Requirements
c++ compiler; compilation was tested with gcc version GCC/7.3.0-2.30 (Linux) and clang-1000.11.45.5 (Mac OSX).
cmake 3.2+.

#### Specific component source code versions used to build Wengan v0.2

1. abyss commit [d4b4b5d](https://github.com/bcgsc/abyss/tree/d4b4b5d3091d90a4967180d987bd7168dbf04585)
2. discovarexp-51885 commit  [f827bab](https://github.com/adigenova/discovarexp-51885/tree/f827bab9bd0e328fee3dd57b7fefebfeebd92be4)
3. minia commit [017d23e](https://github.com/GATB/minia/tree/017d23e60d56db183c499bb2241345e95514ebbe)
4. fastmin-sg commit [861b061](https://github.com/adigenova/fastmin-sg/tree/861b061ecf53cb1aff1c485bac9aca17c399bb7f)
5. intervalmiss commit [11be8b42](https://github.com/adigenova/intervalmiss/tree/11be8b428ef31ae7ad0b94d20c216f4a3d608685)
6. liger commit [63a044b0](https://github.com/adigenova/liger/tree/63a044b04c910d41f28c6196d7d027d9cee3353a)
7. seqtk commit [2efd0c8](https://github.com/adigenova/seqtk/tree/2efd0c85767b2e8ae2366d7ea7edb8041adb0eb1)


# Limitations

    1.- Genomes larger than 4Gb are not supported yet.
    
# About the name
**Wengan** is a [Mapudungun](https://en.wikipedia.org/wiki/Mapuche_language) word. Mapudungun is the language of the [**Mapuche**](https://en.wikipedia.org/wiki/Mapuche) people, the largest indigenous inhabitants of south-central Chile. **Wengan** means "***Making the path***".

# Citation
Di Genova, A., Buena-Atienza, E., Ossowski, S. and Sagot,M-F. **Efficient hybrid de novo assembly of human genomes with WENGAN**. Nature Biotechnology (2020), [link](https://doi.org/10.1038/s41587-020-00747-w)
