# crRNA library (2021)
Repository containing scripts and the pipeline followed to select a collection of *C. elegans* CRISPR RNAs (crRNAs) targets with multiple purposes.

## Basic gist
The main use of this repository is to produce collections of *C. elegans* crRNAs for a given use, *e.g.* N or C-terminus tagging, considering their possible off-target sequences, cutting efficiency, location within a gene, and possible cut between protein domains (as predicted by [Alphafold Protein Structure Database](https://alphafold.ebi.ac.uk/)) as illustrated in the following diagram and section.

![Screenshot](img/crRNA_library-Main_gist.png)

This repository contains only the scripts for the selection while its usage is described elsewhere, for example, in the [crRNA library track selection app](https://github.com/AmhedVargas/CRISPR_library_track).

## Considerations
### crRNA selection
All possible *C. elegans* crRNA target sequences were scrutinized, filtered and assigned to at least one of the following categories:

**1. ATG_500:** crRNAs targeting the 500bp upstream region of a genebody with homology arms next to the expected cut site
    
**2. ATG_250:** crRNAs targeting the 250bp upstream region of a genebody with homology arms around the expected cut site
    
**3. ATG:** crRNAs targeting the ATG of each isoform in around a 10bp window with homology arms starting after the ATG (in frame)
    
**4a. CDS guidescan:** primary source of crRNAs with info and cutting scores defined by the program guidescan. We obtained these scores by running a container with its software pre-installed and querying all posible sites accross the *C. elegans* genome.
    
**4b. CDS hamming:** secondary source of cRNAs where the scoring system was the least number of offtargets. This tracks was mostly used in case of gene duplications.
    **Note as well that for both CDSs crRNA selection in "Structured CDS" (CDS with median alphafold score values over 90) the resulting cRNA selected was the one with the lower alphafold value predicted at cutting site.
    Also, note that we targetted CDSs instead of isoforms. That made easier the computation needed to select the crRNAs, but made difficult to track repeated crRNA targets (which ultimately were collapsed if the type or their homology arms were the same)**
    
**5. Stop:** crRNAs targetting the end of a genebody with homology arms splitting before Stop codon
    
**6. lincRNA:** crRNAs targeting long no coding RNAs
    
**7. miRNAs:** crRNAs targetting microRNA sites

Please note that for genes with multiple start or stop sites, we also tried to select corresponding crRNAs (*i.e.* ATG_500, ATG_250, ATG, and Stop).

Then for each crRNA in the selection, adjacent sequences were added and treated as homology arms for CRISPR experiments. In these arms, restriction sites as re-appearing crRNA sites were removed by altering their sequences (with synonym mutations if located in coding regions).

### Oligo ordering
1. The first iteration did not filled 4*384 plates, so we arranged the last plate in other ways to complete the order:
    
    -Last well of each plate (384) contains sequences used as co-CRISPR
    
    -On plate 4, well 379 to 382 contains a set of selected crRNAs on genes that are used as controls.
    
    -Well 383 is a duplication of co-CRIPSR oligos, and so well 383 and 384 are the same for plate 4.



