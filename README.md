# crRNA library (2021)
Repository containing scripts and the pipeline followed to select a collection of *C. elegans* CRISPR RNAs (crRNAs) targets with multiple purposes.

## Basic gist
The main use of this repository is to produce collections of *C. elegans* crRNAs for a given use, *e.g.* N or C-terminus tagging, considering their possible off-target sequences, cutting efficiency, location within a gene, and possible cut between protein domains (as predicted by [Alphafold Protein Structure Database](https://alphafold.ebi.ac.uk/)) as ilustraded in the following diagram.

![Screenshot](img/crRNA_library-Main_gist.png)

This repository contains only the scripts for the selection while its usage is described elsewhere, for example, in the [crRNA library track selection app](https://github.com/AmhedVargas/CRISPR_library_track).

## Considerations
### crRNA selection
1. We got different categories for the selection, the main being:
    
    -ATG_500: crRNAs targeting the 500bp upstream region of a genebody with homology arms next to the expected cut site
    
    -ATG_250: crRNAs targeting the 250bp upstream region of a genebody with homology arms around the expected cut site
    
    -ATG: crRNAs targeting the ATG of each isoform in around a 10bp window with homology arms starting after the ATG (in frame)
    
    -CDS guidescan: primary source of crRNAs with info and cutting scores defined by the program guidescan. We obtained these scores by running a container with its software pre-installed and querying all posible sites accross the *C. elegans* genome.
    
    -CDS hamming: secondary source of cRNAs where the scoring system was the least number of offtargets. This tracks was mostly used in case of gene duplications.
    **Note as well that for both CDSs crRNA selection in "Structured CDS" (CDS with median alfaphold score values over 90) the resulting cRNA selected was the one with the lower alphafold value predicted at cutting site.
    Also, note that we targetted CDSs instead of isoforms. That made easier the computation needed to select the crRNAs, but made difficult to track repeated crRNA targets (which ultimately were collapsed if the type or their homology arms were the same)**
    
    -Stop: crRNAs targetting the end of a genebody with homology arms splitting before Stop codon
    
    -lincRNA: crRNAs targeting long no coding RNAs
    
    -miRNAs: crRNAs targetting microRNA sites
    
2. crRNAs can belong to different categories as they were somewhat collapsed by position. 

3. ATGs, STOPs, and ATG_250s and ATG_500s crRNA were selected per isoform.

### Synthesis order
1. The first iteration did not filled 4*384 plates, so we arranged the last plate in other ways to complete the order:
    
    -Last well of each plate (384) contains sequences used as co-CRISPR
    
    -On plate 4, well 379 to 382 contains a set of selected crRNAs on genes that are used as controls.
    
    -Well 383 is a duplication of co-CRIPSR oligos, and so well 383 and 384 are the same for plate 4.



