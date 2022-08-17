# crRNA library (2021)
Repository containing scripts and the pipeline followed to select a collection of *C. elegans* CRISPR RNAs (crRNAs) targets with multiple purposes.

## Basic gist
The main use of this repository is to produce collections of *C. elegans* crRNAs for a given use, *e.g.* N or C-terminus tagging, considering their possible off-target sequences, cutting efficiency, location within a gene, and possible cut between protein domains (as predicted by [Alphafold Protein Structure Database](https://alphafold.ebi.ac.uk/)) as illustrated in the following diagram and section.

![Screenshot](img/crRNA_library-Main_gist.png)

This repository contains only the scripts for the selection while its usage is described elsewhere, for example, in the [crRNA library track selection app](https://github.com/AmhedVargas/CRISPR_library_track).

## Considerations
### crRNA selection
All possible *C. elegans* crRNA target sequences were scrutinized, filtered, and assigned to at least one of the following categories:

**1. ATG_500:** crRNAs targeting the 500 bp upstream region of a genebody with homology arms next to the expected cut site.
    
**2. ATG_250:** crRNAs targeting the 250 bp upstream region of a genebody with homology arms around the expected cut site.
    
**3. ATG:** crRNAs targeting the ATG of each isoform in around a 10bp window with homology arms starting after the start codon (in frame).

**4. CDS:** crRNAs targeting each possible CDS of a given gene. For CDS predicted to encode for a stable protein domain with high confidence, i.e. Alphafold's pLDTT median score above 90, we prioritize the selection of crRNA with the lowest pLDTT score in its cutting site. With this criteria we hoped to avoid disruption of sensitive regions for protein activity, however, note that studies in *drosophila* have revealed that efficient internal protein fluorescent tagging can occur irrespective of the insertion site (*e.g.* [Nagarkar-Jaiswal *et al.* 2015](https://elifesciences.org/articles/5338)).

**4a. CDS guidescan:** crRNAs obtained from the program guidescan. We prioritized these ones as scores for cutting efficiency and on-site activity were available. Note that we obtained these scores by running a [docker container](https://hub.docker.com/layers/xerez/guidescan/latest/images/sha256-ea5c5ed0b873205243babb26a49f85f14f2c05fd992e66f6ff13722842df9ef7) with its software pre-installed and querying all possible sites across the *C. elegans* genome.
    
**4b. CDS hamming:** For genes with no unique crRNA site (*i.e.* 0 off-targets) such as recent gene duplications, we selected crRNAs with low off-targets as calculated by the hamming distance (*i.e.* a metric that calculate similarity of sequences by determining how many base substitutions are required to go from one to the another) to other sites to the genome. 
    
**5. Stop:** crRNAs targeting the end of a genebody with homology arms adjacent to the beginning of the stop codon.
    
**6. lincRNA:** crRNAs targeting long no coding RNAs.
    
**7. miRNAs:** crRNAs targeting microRNA sites.

Please note that for genes with multiple start or stop sites, we also tried to select corresponding crRNAs (*i.e.* ATG_500, ATG_250, ATG, and Stop).

Then for each crRNA in the selection, adjacent sequences were added and treated as homology arms for CRISPR experiments. In these arms, restriction sites as re-appearing crRNA sites were removed by altering their sequences (with synonym mutations if located in coding regions).

### Ordering schema within plates
The crRNA library was designed to be composed of **four** 384 plates each of them with 120 oligos. Read Al-johani *et al.* for further information. The final arrangement of our plates had the following considerations:
    
* Last well of each plate (384) contains sequences used as co-CRISPR.
    
* On plate 4, well 379 to 382 contains a set of selected crRNAs on genes that are used as controls.
    
* Well 383 is a duplication of co-CRIPSR oligos, and so well 383 and 384 are the same for plate 4.

## Library generation
To ease and automate the crRNA selection process, we figured out that a good strategy was to produce genomic tracks from where we could select crRNAs and optimize further our selections. Therefore we focused in the creation of:

**1. Genebodies and CDS tracks:** For each *C. elegans* protein coding gene annotated in the WS283 version of [WormBase](https://wormbase.org/) we extracted every annotated CDS, collapsed the repeated CDS, and created a genomic track with them. Similarly, we annotated the location of their ATGs and stop codons. 

**2. An Alphafold confidence track:** Predicted confidence scores were downloaded for all available *C. elegans* proteins from the [Alphafold Protein Structure Database](https://alphafold.ebi.ac.uk/). Corresponding WormBase transcripts with the same coding sequence length were matched and location extracted from WormBase annotations. For each codon a corresponding confidence score and genomic location was assigned and written into a single bedgraph per transcript. A final bedgraph with all *C. elegans* annotated transcripts was produced and used for further analysis such as, defining if CDS was "structured" (confidence score above 90) and what score had the cutting site of crRNAs.

**3. crRNAs tracks:** crRNA genomic tracks were obtained from two different sources. The first set comes directly from [guidescan](https://guidescan.com/) calculations for ce11 and obtained trhough running a batch query within its [docker distribution](https://hub.docker.com/layers/xerez/guidescan/latest/images/sha256-ea5c5ed0b873205243babb26a49f85f14f2c05fd992e66f6ff13722842df9ef7). The second set was consisted of every possible crRNA with NGG pam site within *C. elegans* ce11/WS235 genome with additional info of their uniqueness within the genome (calculated by a [hamming distance](https://github.com/AmhedVargas/CelegansHammingAlignments)).

**Scripts to produce each track can be found in the `tracks` folder of this repository.**

Once we had our different genomic tracks, we proceeded to intersect them to obtain selection metrics for each crRNA. Selection of crRNAs to target as many *C. elegans* genes as possible followed these guidelines:

* Instead of working out genebodies, we aimed to select a single crRNA per unique CDS of a given gene.

* For each crRNA with a cutting site within coding regions, we assigned them an alphafold confidence score. If a CDS had a median confidence score above 90, we picked the crRNA with the lowest confidence score. Otherwise, the crRNA with the highest guidescan on-site efficiency and lowest off-site targeting was selected.

* If a gene had zero crRNAs targets from the guidescan set, we proceeded to use crRNAs with hamming distances instead prioritizing those with lower off-target effects as calculated via hamming distance.

* Also, we selected the closest crRNA to the start and stop codon of each gene isoform but only if crRNAs with a cutting site closer to 10 bp to either the start or stop codon existed. Similarly, we took the crRNA with the closest cutting site to 250 and 500 base pairs to the start of a gene (which we called ATG_250 and ATG_500).

* Additionally, we selected crRNAs that targeted the genomic locations of long non-coding RNAs and microRNAs.

For every crRNA, we obtained adjacent sequences to their cutting site which could be used as homology arms for CRISPR experiments. In order to avoid re-cutting of the crRNA and facilitate our cloning experiments with restriction enzymes (RE), we modified the sequence of the homology arms by adding three synonymous mutations if the crRNA or RE target sequence was in a coding region and one at the PAM or at the RE recognition site if not. 

**The above steps were performed using the code located in the `processing` folder of this repository.**

Finally, we proceeded to assembly the crRNA selection with their homology arms into a scaffold, and fill 384 plates with these scaffold using the crRNA target genomic location as ordering schema. From the resulting library, a bed file that was used to develop a [shiny app](https://wormbuilder.dev/crRNALib/). 

**The code for the library assembly can be found at the `assembly` folder of this repository.**

## Troubleshoot

Please feel free to [e-mail me](mailto:amhed.velazquez@kaust.edu.sa) for any question, doubt or error in the code.

## Citation
For now there is no available publication
