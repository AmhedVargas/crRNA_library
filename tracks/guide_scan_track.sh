sudo docker run --rm -it -v $(pwd):/data xerez/guidescan guidescan_guidequery -b /data/cas9_ce11_all_guides.bam -c chrX:1-1000000 --target within -o /data/res
sudo docker run --rm -it -v $(pwd):/data xerez/guidescan guidescan_guidequery -b /data/cas9_ce11_all_guides.bam -c chrX:1000000-2000000 --target within -o /data/res
sudo docker run --rm -it -v $(pwd):/data xerez/guidescan guidescan_guidequery -b /data/cas9_ce11_all_guides.bam -c chrX:2000000-3000000 --target within -o /data/res
sudo docker run --rm -it -v $(pwd):/data xerez/guidescan guidescan_guidequery -b /data/cas9_ce11_all_guides.bam -c chrX:3000000-4000000 --target within -o /data/res
sudo docker run --rm -it -v $(pwd):/data xerez/guidescan guidescan_guidequery -b /data/cas9_ce11_all_guides.bam -c chrX:4000000-5000000 --target within -o /data/res
sudo docker run --rm -it -v $(pwd):/data xerez/guidescan guidescan_guidequery -b /data/cas9_ce11_all_guides.bam -c chrX:5000000-6000000 --target within -o /data/res
sudo docker run --rm -it -v $(pwd):/data xerez/guidescan guidescan_guidequery -b /data/cas9_ce11_all_guides.bam -c chrX:6000000-7000000 --target within -o /data/res
sudo docker run --rm -it -v $(pwd):/data xerez/guidescan guidescan_guidequery -b /data/cas9_ce11_all_guides.bam -c chrX:7000000-8000000 --target within -o /data/res
sudo docker run --rm -it -v $(pwd):/data xerez/guidescan guidescan_guidequery -b /data/cas9_ce11_all_guides.bam -c chrX:8000000-9000000 --target within -o /data/res
sudo docker run --rm -it -v $(pwd):/data xerez/guidescan guidescan_guidequery -b /data/cas9_ce11_all_guides.bam -c chrX:9000000-10000000 --target within -o /data/res
sudo docker run --rm -it -v $(pwd):/data xerez/guidescan guidescan_guidequery -b /data/cas9_ce11_all_guides.bam -c chrX:10000000-11000000 --target within -o /data/res

sudo docker run --rm -it -v $(pwd):/data xerez/guidescan guidescan_guidequery -b /data/cas9_ce11_all_guides.bam -c chrX:11000000-12000000 --target within -o /data/res
sudo docker run --rm -it -v $(pwd):/data xerez/guidescan guidescan_guidequery -b /data/cas9_ce11_all_guides.bam -c chrX:12000000-13000000 --target within -o /data/res
sudo docker run --rm -it -v $(pwd):/data xerez/guidescan guidescan_guidequery -b /data/cas9_ce11_all_guides.bam -c chrX:13000000-14000000 --target within -o /data/res
sudo docker run --rm -it -v $(pwd):/data xerez/guidescan guidescan_guidequery -b /data/cas9_ce11_all_guides.bam -c chrX:14000000-15000000 --target within -o /data/res
sudo docker run --rm -it -v $(pwd):/data xerez/guidescan guidescan_guidequery -b /data/cas9_ce11_all_guides.bam -c chrX:15000000-16000000 --target within -o /data/res
sudo docker run --rm -it -v $(pwd):/data xerez/guidescan guidescan_guidequery -b /data/cas9_ce11_all_guides.bam -c chrX:16000000-17000000 --target within -o /data/res
sudo docker run --rm -it -v $(pwd):/data xerez/guidescan guidescan_guidequery -b /data/cas9_ce11_all_guides.bam -c chrX:17000000-17718942 --target within -o /data/res




awk -F"\t" '{if(NR==1){print $0}else{OFS="\t"; if($6=="+"){print $1,$2,($3-3),$4,$5,$6}else{print $1,($2+3),($3),$4,$5,$6}}}' cas9_ce11_all_guides.sfulike.bed | sort -k1,1 -k2,2n > tmp.bed

awk -F"\t" '{if(NR==1){print $0}else{OFS="\t"; if($6=="+"){ print $1,($3-3),($3-3),$4,$5,$6,$7,$8}else{print $1,($2+3),($2+3),$4,$5,$6,$7,$8}}}' NGG_guides_WS250.bed | sort -k1,1 -k2,2n > test.bed

## add trtacks on
echo track gffTags=\"on\" > on.txt

grep -v "track" cas9_ce11_all_guides.sfulike.bed | awk -F"\t" '{{OFS="\t"; if($6=="+"){print $1,$2,($3-3),$4,$5,$6}else{print $1,($2+3),($3),$4,$5,$6}}}' - | sort -k1,1 -k2,2n | cat on.txt - > cas9_ce11_all_guides.sfulike.fixed.bed

##Cutting site
awk -F"\t" '{if(NR==1){print $0}else{OFS="\t"; if($6=="+"){ print $1,($3-3),($3-3),$4,$5,$6,$7,$8}else{print $1,($2+3),($2+3),$4,$5,$6,$7,$8}}}' NGG_guides_WS250.bed | sort -k1,1 -k2,2n > test.bed

grep -v "track" cas9_ce11_all_guides. | awk -F"\t" '{OFS="\t"; if($6=="+"){ print $1,($3-3),($3-3),"Start="$2";End="$3";Strand="$6";"$4,$5,$6}else{print $1,($2+3),($2+3),"Start="$2";End="$3";Strand="$6";"$4,$5,$6}}' | sort -k1,1 -k2,2n | cat on.txt - > NGG_guides_WS250.cutting.bed

grep -v "track" cas9_ce11_all_guides.sfulike.fixed.bed | awk -F"\t" '{OFS="\t"; if($6=="+"){ print $1,($3-3),($3-3),"Start="$2";End="$3";Strand="$6";"$4,$5,$6}else{print $1,($2+3),($2+3),"Start="$2";End="$3";Strand="$6";"$4,$5,$6}}' | sort -k1,1 -k2,2n | cat on.txt - > cas9_ce11_all_guides.cutting.all.bed





