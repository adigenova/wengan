# Wengan
An accurate and ultra-fast genome assembler

# SYNOPSIS

    # Assembling Oxford nanopore and illumina reads wiht WenganM
     wengan.pl -x ontraw -a M -s lib1.fwd.fastq.gz,lib1.rev.fastq.gz -l ont.fastq.gz -p asm1 -t 20 -g 3000

    # Assembling PacBio reads and illumina reads with WenganA
     wengan.pl -x pacraw -a A -s lib1.fwd.fastq.gz,lib1.rev.fastq.gz -l pac.fastq.gz -p asm2 -t 20 -g 3000

    # Assembling ultra-long nanopore reads and BGI reads with WenganM
     wengan.pl -x ontlon -a M -s lib2.fwd.fastq.gz,lib2.rev.fastq.gz -l ont.fastq.gz -p asm3 -t 20 -g 3000

    # Non-hybrid assembly of PacBio Circular Consensus Sequence data with WenganM
     wengan.pl -x pacccs -a M -l ccs.fastq.gz -p asm4 -t 20 -g 3000

    # Assembling ultra-long nanopore reads and Illumina reads with WenganD (requires a high memory machine 600Gb)
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

preset for raw ultra-long-reads from Oxford Nanopore, tipically having an  N50 > 50kb.

## ontraw

preset for raw long-reads Nanopore reads tipically having an  N50 ~\[15kb-40kb\].

## pacraw

preset for raw long-reads from Pacific Bioscience (PacBio) tipically having an  N50 ~\[8kb-60kb\].

## pacccs

preset for Circular Consensus Sequences from Pacific Bioscience (PacBio) tipically having an  N50 ~\[15kb\].

# WENGAN ADVANCED OPTIONS

The following options allows to override the presets of **Wengan** components.
Don't change this variables if you are not sure.

## FastMin-SG

An alignment-free algorithm for ultrafast scaffolding graph construction from short or long reads.

## IntervalMiss

IntervalMiss detect miss-assembled contigs and correct them when necessary.

## Liger

Liger use the Synthetic Scaffoding Graph to compute overlap among long reads,
order and orient short contigs, validate scaffols sequences, fill the gaps and
polishing of the assembly.

# Wengan components
+ A bruijn graph assembler ([Minia](https://github.com/GATB/minia), [Abyss](https://github.com/bcgsc/abyss) or [DiscovarDenovo](https://software.broadinstitute.org/software/discovar/blog/))
+ [FastMIN-SG](https://github.com/adigenova/fastmin-sg)
+ [IntervalMiss](https://github.com/adigenova/intervalmiss)
+ [Liger](https://github.com/adigenova/liger)

<img src="./wengan-diagram.svg">

# About the name
**Wengan** is a [Mapudungun](https://en.wikipedia.org/wiki/Mapuche_language) word. The Mapudungun is the language of the [**Mapuche**](https://en.wikipedia.org/wiki/Mapuche) people, the largest indigenous inhabitants of south-central Chile. **Wengan** means "***Making the path***".

