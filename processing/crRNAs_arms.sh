for file in `ls Collapsed_crRNAs*`; do echo ${file%.bed}; bedtools intersect -a Genebodies_p11.bed -b ${file} -wo > ${file%.bed}.Intersections; perl Process_arms.pl Genebodies_data.hom.2.tsv ${file%.bed}.Intersections 4 > ${file%.bed}.arms.tsv; done

awk -F"\t" '{split($4,info,";"); split($10,data,";"); if(info[2]==data[(length(data)-1)]){print $0}}' ../Collapsed_crRNAs.ATG.Intersections > Collapsed_crRNAs.ATG.Intersections

perl ../Process_arms.pl ../Genebodies_data.hom.2.tsv Collapsed_crRNAs.ATG.Intersections 4 > Collapsed_crRNAs.ATG.arms.tsv

 paste Collapsed_crRNAs.ATG.Intersections Collapsed_crRNAs.ATG.arms.tsv | awk -F"\t" '{OFS="\t"; split($10,info,";"); crRNA=info[4]; if($14==info[4]){print $7,$8,$9,$10,$11,$12,$15,$16}}' - > Collapsed_crRNAs.ATG.wArms.bed


 awk -F"\t" '{split($4,info,";"); split($10,data,";"); if(info[2]==data[(length(data)-1)]){print $0}}' ../Collapsed_crRNAs.TAA.Intersections > Collapsed_crRNAs.TAA.Intersections

 perl ../Process_arms.pl ../Genebodies_data.hom.2.tsv Collapsed_crRNAs.TAA.Intersections 4 > Collapsed_crRNAs.TAA.arms.tsv

 paste Collapsed_crRNAs.TAA.Intersections Collapsed_crRNAs.TAA.arms.tsv | awk -F"\t" '{OFS="\t"; split($10,info,";"); crRNA=info[4]; if($14==info[4]){print $7,$8,$9,$10,$11,$12,$15,$16}}' - > Collapsed_crRNAs.TAA.wArms.bed

 awk -F"\t" '{split($4,info,";"); split($10,data,";"); if(info[2]==data[(length(data)-2)]){print $0}}' ../Collapsed_crRNAs.CDS.2.Intersections > Collapsed_crRNAs.CDS.2.Intersections

 perl ../Process_arms.pl ../Genebodies_data.hom.2.tsv Collapsed_crRNAs.CDS.2.Intersections 4 > Collapsed_crRNAs.CDS.2.arms.tsv

 paste Collapsed_crRNAs.CDS.2.Intersections Collapsed_crRNAs.CDS.2.arms.tsv | awk -F"\t" '{OFS="\t"; split($10,info,";"); crRNA=info[4]; if($14==info[4]){print $7,$8,$9,$10,$11,$12,$15,$16}}' - > Collapsed_crRNAs.CDS.2.wArms.bed

 awk -F"\t" '{split($4,info,";"); split($10,data,";"); if(info[2]==data[(length(data)-2)]){print $0}}' ../Collapsed_crRNAs.Hamming.CDSs.2.Intersections > Collapsed_crRNAs.Hamming.CDSs.2.Intersections

 perl ../Process_arms.pl ../Genebodies_data.hom.2.tsv Collapsed_crRNAs.Hamming.CDSs.2.Intersections 4 > Collapsed_crRNAs.Hamming.CDSs.2.arms.tsv










bedtools intersect -a Cel_CDs_medianAF.bed -b guideScan_allinfo.removed.simplify.cutting.AF.bed -wo | awk -F"\t" '{
    OFS="\t"; 
    if(NR==1){
        key=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6; split("",hash,""); 
    }; 
    newkey=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6; 
    if(newkey != key){
        split(key,info,"\t"); 
        if(info[5]>85){
            score=0; lowscore=90; pick=""; 
            for(cr in hash){
                split(cr,crinfo,";"); split(crinfo[5],pp,"="); split(crinfo[length(crinfo)],afscore,"="); 
                if(afscore[2] < lowscore){
                    lowscore=afscore[2]; pick=cr;
                };
            }; 
            if(pick ==""){
                score=0; pick="";  
                for(cr in hash){
                    split(cr,crinfo,";"); split(crinfo[5],pp,"="); 
                    if(pp[2] > score){
                        pick=cr; score = pp[2]
                    };
                }
            }; 
            res[key"\t"pick]++; used[pick]++; key=newkey; split("",hash,"");hash[$10]=$11
        }else{
            score=0; pick="";  for(cr in hash){
                split(cr,crinfo,";"); split(crinfo[5],pp,"=");
                if(pp[2] > score){
                    pick=cr; score = pp[2]
                };
            }; 
            res[key"\t"pick]++; used[pick]++; hash[$10]=$11
        }; 
        key=newkey; split("",hash,""); hash[$10]=$11
        }else{
            hash[$10]=$11
            }
        } 
        END{for(key in res){print key}}' | awk -F"\t" '{OFS="\t"; split($7,info,";"); split(info[1],start,"="); split(info[2],end,"="); split(info[3],strand,"="); print $1,start[2],end[2],$7";"$4,".",strand[2]}' | sort -k1,1 -k2,2n > First_crRNA_selection.90.bed


bedtools intersect -a Cel_CDs_medianAF.bed -b guideScan_allinfo.removed.simplify.cutting.AF.bed -wo | awk -F"\t" '{OFS="\t"; key=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6; struct=$5; crRNA=$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12;  split($10,info,";"); split(info[5],data,"="); cute=data[2]; split(info[length(info)],data,"="); afscore=data[2]; if(res[key]==""){res[key]=crRNA; eff[key]=cute; scor[key]=afscore}else{ if(struct > 85){ if(scor[key]> afscore){res[key]=crRNA; eff[key]=cute; scor[key]=afscore}}else{if(eff[key]<cute){res[key]=crRNA; eff[key]=cute; scor[key]=afscore}} }} END{for(key in res){print key"\t"res[key]}}' | sort -k1,1 -k2,2n > cRNAselection


###Now intersect for non coding RNAs
bedtools intersect -a miRNA.bed -b guideScan_allinfo.removed.simplify.cutting.AF.bed -wo | awk -F"\t" '{OFS="\t"; key=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6; struct=$5; crRNA=$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12;  split($10,info,";"); split(info[5],data,"="); cute=data[2]; split(info[length(info)],data,"="); afscore=data[2]; if(res[key]==""){res[key]=crRNA; eff[key]=cute; scor[key]=afscore}else{ if(struct > 85){ if(scor[key]> afscore){res[key]=crRNA; eff[key]=cute; scor[key]=afscore}}else{if(eff[key]<cute){res[key]=crRNA; eff[key]=cute; scor[key]=afscore}} }} END{for(key in res){print key"\t"res[key]}}' | sort -k1,1 -k2,2n > miRNA.crRNA.selection.tsv


bedtools intersect -a lincRNA.bed -b guideScan_allinfo.removed.simplify.cutting.AF.bed -wo | awk -F"\t" '{OFS="\t"; key=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6; struct=$5; crRNA=$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12;  split($10,info,";"); split(info[5],data,"="); cute=data[2]; split(info[length(info)],data,"="); afscore=data[2]; if(res[key]==""){res[key]=crRNA; eff[key]=cute; scor[key]=afscore}else{ if(struct > 85){ if(scor[key]> afscore){res[key]=crRNA; eff[key]=cute; scor[key]=afscore}}else{if(eff[key]<cute){res[key]=crRNA; eff[key]=cute; scor[key]=afscore}} }} END{for(key in res){print key"\t"res[key]}}' | sort -k1,1 -k2,2n > lincRNA.crRNA.selection.tsv


bedtools intersect -a ncRNA.bed -b guideScan_allinfo.removed.simplify.cutting.AF.bed -wo | awk -F"\t" '{OFS="\t"; key=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6; struct=$5; crRNA=$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12;  split($10,info,";"); split(info[5],data,"="); cute=data[2]; split(info[length(info)],data,"="); afscore=data[2]; if(res[key]==""){res[key]=crRNA; eff[key]=cute; scor[key]=afscore}else{ if(struct > 85){ if(scor[key]> afscore){res[key]=crRNA; eff[key]=cute; scor[key]=afscore}}else{if(eff[key]<cute){res[key]=crRNA; eff[key]=cute; scor[key]=afscore}} }} END{for(key in res){print key"\t"res[key]}}' | sort -k1,1 -k2,2n > ncRNA.crRNA.selection.tsv

###Now get bed files iwht spacers and proper strand to make arms 45
for file in `ls *.selection.tsv`; do echo ${file}; 
awk -F"\t" '{OFS="\t"; split($10,info,";"); split(info[1],start,"="); split(info[2],end,"="); split(info[3],strand,"="); print $1,start[2],end[2],$10";"$4";Type='${file%.crRNA.selection.tsv}'",$11,strand[2]}' ${file}  | sort -k1,1 -k2,2n > ${file%.tsv}.bed;
awk -F"\t" '{OFS="\t"; strand=$6; pos=$8; if(strand=="+"){print $1,(pos-45),pos,$10";"$4";Type='${file%.crRNA.selection.tsv}'",$5,strand > "'${file%.tsv}.left.arm.bed'"; print $1,(pos),(pos+45),$10";"$4";Type='${file%.crRNA.selection.tsv}'",$5,strand > "'${file%.tsv}.right.arm.bed'";}else{print $1,(pos-45),pos,$10";"$4";Type='${file%.crRNA.selection.tsv}'",$5,strand > "'${file%.tsv}.right.arm.bed'"; print $1,(pos),(pos+45),$10";"$4";Type='${file%.crRNA.selection.tsv}'",$5,strand > "'${file%.tsv}.left.arm.bed'";}}' ${file}
bedtools getfasta -fi c_elegans.PRJNA13758.WS282.genomic.fa -bed ${file%.tsv}.left.arm.bed -s -tab -name > left
bedtools getfasta -fi c_elegans.PRJNA13758.WS282.genomic.fa -bed ${file%.tsv}.right.arm.bed -s | grep -v ">" > right
paste left right | perl -pe 's/\([+-]\)//g' | awk -F"\t" '{if(NF==3){hash[$1]=$2"\t"$3}else{ OFS="\t"; print $0,hash[$4]}}' - ${file%.tsv}.bed > ${file%.tsv}.wArms.bed;
rm left right
done



bedtools intersect -a Genebodies_p11.bed -b First_crRNA_selection.90.bed -wao | awk -F"\t" '{split($4,info,";"); split($10,data,";"); if(info[2]==data[(length(data)-1)]){print $0}}' 


bedtools intersect -a Genebodies_p11.bed -b First_crRNA_selection.90.bed -wo | awk -F"\t" '{split($4,info,";"); split($10,data,";"); if(info[2]==data[(length(data)-1)]){print $0}}' > Collapsed_crRNAs.CDS.Intersections
perl Process_arms4Terminal.2.pl Genebodies_data.hom.2.tsv Collapsed_crRNAs.CDS.Intersections 4 CDS > Collapsed_crRNAs.CDS.arms.tsv
awk -F"\t" '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10";Type=CDSguiscan",$11,$12,$13}' Collapsed_crRNAs.CDS.Intersections > tmp
mv tmp Collapsed_crRNAs.CDS.Intersections
paste Collapsed_crRNAs.CDS.Intersections Collapsed_crRNAs.CDS.arms.tsv | awk -F"\t" '{OFS="\t"; split($10,info,";"); crRNA=info[4]; if($14==info[4]){print $7,$8,$9,$10,$11,$12,$15,$16}}' - > Collapsed_crRNAs.CDS.wArms.bed



