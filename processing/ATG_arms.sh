awk -F"\t" '{OFS="\t"; strand= $6; split($4,coso,"WBGene"); split(coso[2],data,";"); if(strand=="+"){pos=($3-3)}else{pos=($2+3)}; print $1,(pos-45),pos,$4,$5,$6,"WBGene"data[1]}' Collapsed_crRNAs.ATG.250.bed | sort -k1,1 -k2,2n > Collapsed_crRNAs.leftArm.ATG.250.bed

awk -F"\t" '{OFS="\t"; strand= $6; split($4,coso,"WBGene"); split(coso[2],data,";"); if(strand=="+"){pos=($3-3)}else{pos=($2+3)}; print $1,(pos-45),pos,$4,$5,$6,"WBGene"data[1]}' Collapsed_crRNAs.ATG.500.bed | sort -k1,1 -k2,2n > Collapsed_crRNAs.leftArm.ATG.500.bed

awk -F"\t" '{OFS="\t"; strand= $6; split($4,coso,"WBGene"); split(coso[2],data,";"); if(strand=="+"){pos=($3-3)}else{pos=($2+3)}; print $1,(pos),(pos+45),$4,$5,$6,"WBGene"data[1]}' Collapsed_crRNAs.ATG.500.bed | sort -k1,1 -k2,2n > Collapsed_crRNAs.rightArm.ATG.500.bed

awk -F"\t" '{OFS="\t"; strand= $6; split($4,coso,"WBGene"); split(coso[2],data,";"); if(strand=="+"){pos=($3-3)}else{pos=($2+3)}; print $1,(pos),(pos+45),$4,$5,$6,"WBGene"data[1]}' Collapsed_crRNAs.ATG.250.bed | sort -k1,1 -k2,2n > Collapsed_crRNAs.rightArm.ATG.250.bed

awk -F"\t" '{if(NF==2){hash[$1]=$2}else{OFS="\t"; print $1,$2,$3,$4,$5,hash[$7]}}' Genes_strands Collapsed_crRNAs.rightArm.ATG.500.bed > tmp

mv tmp Collapsed_crRNAs.rightArm.ATG.500.bed

awk -F"\t" '{if(NF==2){hash[$1]=$2}else{OFS="\t"; print $1,$2,$3,$4,$5,hash[$7]}}' Genes_strands Collapsed_crRNAs.rightArm.ATG.250.bed > tmp

mv tmp Collapsed_crRNAs.rightArm.ATG.250.bed

awk -F"\t" '{if(NF==2){hash[$1]=$2}else{OFS="\t"; print $1,$2,$3,$4,$5,hash[$7]}}' Genes_strands Collapsed_crRNAs.leftArm.ATG.250.bed > tmp

mv tmp Collapsed_crRNAs.leftArm.ATG.250.bed

awk -F"\t" '{if(NF==2){hash[$1]=$2}else{OFS="\t"; print $1,$2,$3,$4,$5,hash[$7]}}' Genes_strands Collapsed_crRNAs.leftArm.ATG.500.bed > tmp

mv tmp Collapsed_crRNAs.leftArm.ATG.500.bed


paste Collapsed_crRNAs.leftArm.ATG.500.bed Collapsed_crRNAs.rightArm.ATG.500.bed | awk -F"\t" '{key1=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6; strand=$6; key2=$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; if(strand=="+"){if($2 < $8){print key1 > "tmpleft"; print key2 > "tmpright";}else{ print key2 > "tmpleft"; print key1 > "tmpright";}}else{ if($2 > $8){print key1 > "tmpleft"; print key2 > " tmpright";}else{ print key2 > "tmpleft"; print key1 > "tmpright";}}}'

mv tmpleft ATG.500.leftArm.bed

mv tmpright ATG.500.rightArm.bed

paste Collapsed_crRNAs.leftArm.ATG.250.bed Collapsed_crRNAs.rightArm.ATG.250.bed | awk -F"\t" '{key1=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6; strand=$6; key2=$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; if(strand=="+"){if($2 < $8){print key1 > "tmpleft"; print key2 > "tmpright";}else{ print key2 > "tmpleft"; print key1 > "tmpright";}}else{ if($2 > $8){print key1 > "tmpleft"; print key2 > " tmpright";}else{ print key2 > "tmpleft"; print key1 > "tmpright";}}}'

mv tmpright ATG.250.rightArm.bed

mv tmpleft ATG.250.leftArm.bed

awk -F"\t" '{OFS="\t"; if($2 < 0 ){$2=0; $3=45}; print $1,$2,$3,$4,$5,$6}' ATG.500.rightArm.bed > tmp
mv tmp ATG.500.rightArm.bed

awk -F"\t" '{OFS="\t"; if($2 < 0 ){$2=0; $3=45}; print $1,$2,$3,$4,$5,$6}' ATG.250.rightArm.bed > tmp
mv tmp ATG.250.rightArm.bed


bedtools getfasta -fi c_elegans.PRJNA13758.WS282.genomic.fa -bed ATG.500.rightArm.bed -s -tab > ATG.500.rightArm.fasta
bedtools getfasta -fi c_elegans.PRJNA13758.WS282.genomic.fa -bed ATG.500.leftArm.bed -s -tab > ATG.500.leftArm.fasta
bedtools getfasta -fi c_elegans.PRJNA13758.WS282.genomic.fa -bed ATG.250.rightArm.bed -s -tab > ATG.250.rightArm.fasta
bedtools getfasta -fi c_elegans.PRJNA13758.WS282.genomic.fa -bed ATG.250.leftArm.bed -s -tab > ATG.250.leftArm.fasta

paste Collapsed_crRNAs.ATG.250.bed ATG.250.leftArm.fasta ATG.250.rightArm.fasta | awk -F"\t" '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$8,$10}' > Collapsed_crRNAs.ATG.250.wArms.bed
paste Collapsed_crRNAs.ATG.500.bed ATG.500.leftArm.fasta ATG.500.rightArm.fasta | awk -F"\t" '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$8,$10}' > Collapsed_crRNAs.ATG.500.wArms.bed




 bedtools intersect -a Cel_CDs_medianAF.bed -b guideScan_allinfo.removed.simplify.cutting.AF.bed -wo | awk -F"\t" '{OFS="\t"; if(NR==1){key=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6; split("",hash,""); }; newkey=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6; if(newkey != key){split(key,info,"\t"); if(info[5]>85){score=0; lowscore=75; pick=""; for(cr in hash){split(cr,crinfo,";"); split(crinfo[5],pp,"="); split(crinfo[length(crinfo)],afscore,"="); if(afscore[2] < lowscore){lowscore=afscore[2]; pick=cr};}; if(pick ==""){score=0; pick="";  for(cr in hash){split(cr,crinfo,";"); split(crinfo[5],pp,"="); if(pp[2] > score){pick=cr; score = pp[2]};}}; res[key"\t"pick]++; used[pick]++; key=newkey; split("",hash,"");hash[$10]=$11}else{score=0; pick="";  for(cr in hash){split(cr,crinfo,";"); split(crinfo[5],pp,"="); if(pp[2] > score){pick=cr; score = pp[2]};}; res[key"\t"pick]++; used[pick]++; hash[$10]=$11}; key=newkey; split("",hash,""); hash[$10]=$11}else{hash[$10]=$11} } END{for(key in res){print key}}' | awk -F"\t" '{OFS="\t"; split($7,info,";"); split(info[1],start,"="); split(info[2],end,"="); split(info[3],strand,"="); print $1,start[2],end[2],$7";"$4,".",strand[2]}' | sort -k1,1 -k2,2n > First_crRNA_selection.tmp.bed






paste Collapsed_crRNAs.ATG.250.bed ATG.250.leftArm.fasta ATG.250.rightArm.fasta | awk -F"\t" '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$8,$10}' > Collapsed_crRNAs.ATG.250.wArms.bed



