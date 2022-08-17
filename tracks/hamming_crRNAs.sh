##Directory in IBEX
#/ibex/scratch/velazqam/sgRNAs/parts\

cd /ibex/scratch/velazqam/sgRNAs/parts

## Right now there is files until 18


cat 00/Histogram.txt 01/Histogram.txt 02/Histogram.txt 03/Histogram.txt 04/Histogram.txt 05/Histogram.txt 06/Histogram.txt 07/Histogram.txt 08/Histogram.txt 09/Histogram.txt 10/Histogram.txt 11/Histogram.txt 12/Histogram.txt 13/Histogram.txt 14/Histogram.txt 15/Histogram.txt 16/Histogram.txt 17/Histogram.txt 18/Histogram.txt | awk -F" " 'BEGIN{flag=0} {if(flag==1){if(length($1) > length(old)){print $1""old" "$2" "$3" "$4" "$5" "$6" "$7; flag=0}else{print $1""old" "data; flag=0}}else{if(length($1)<20){flag=1; old=$1; data=($2" "$3" "$4" "$5" "$6" "$7);}else{print $0}}}' > Seq00_18.txt

###Directory at workstation
cd /home/velazqam/Documents/Projects/Ibex/sgRNAs/res

cp velazqam@dm:/home/velazqam/Desktop/Files/sgRNA_hamming/Seq00_18.txt .

head -n 1275069 ../20_mer/NGG_Cel_WS282.counts |  paste - Seq00_18.txt > Seq00_18.tsv

##Replace unique count with hash done by k-mers in perl
##redo
awk '{print $1}' Seq00_18.txt > Seq00_18.seqs

awk -F"\t" '{split($1,infoa,":"); split(infoa[1],chr,">"); split(infoa[2],infob,"("); split(infob[1],pos,"-"); split(infob[2],infoc,")"); print chr[2]"\t"pos[1]"\t"pos[2]"\t"$2"\t.\t"infoc[1]}' Seq00_18.tsv| head

## COmbine All
awk -F"\t" '{split($1,infoa,":"); split(infoa[1],chr,">"); split(infoa[2],infob,"("); split(infob[1],pos,"-"); split(infob[2],infoc,")"); split($4,mat," "); print chr[2]"\t"pos[1]"\t"pos[2]"\t"$2";"$3";"mat[3]";"mat[4]";"mat[5]";"mat[6]";"mat[7]"\t.\t"infoc[1]}' Seq00_18.tsv > Seq00_18-All.bed

awk -F"\t" '{split($4,info,";"); if(info[2]==1){print $0}}' Seq00_18-All.bed > Seq00_18-0MM.bed
awk -F"\t" '{split($4,info,";"); if(info[2]==1){if(info[3]==0){print $0}}}' Seq00_18-All.bed > Seq00_18-1MM.bed
awk -F"\t" '{split($4,info,";"); if(info[2]==1){if(info[3]==0){if(info[4]==0){print $0}}}}' Seq00_18-All.bed > Seq00_18-2MM.bed
awk -F"\t" '{split($4,info,";"); if(info[2]==1){if(info[3]==0){if(info[4]==0){if(info[5]==0){print $0}}}}}' Seq00_18-All.bed > Seq00_18-3MM.bed
awk -F"\t" '{split($4,info,";"); if(info[2]==1){if(info[3]==0){if(info[4]==0){if(info[5]==0){if(info[6]==0){print $0}}}}}}' Seq00_18-All.bed > Seq00_18-4MM.bed
awk -F"\t" '{split($4,info,";"); if(info[2]==1){if(info[3]==0){if(info[4]==0){if(info[5]==0){if(info[6]==0){if(info[7]==0){print $0}}}}}}}' Seq00_18-All.bed > Seq00_18-5MM.bed


##In IBEX
#Combine results of histograms
for num in {00..99}; do echo $num; cat ${num}/Histogram.txt >> ../All_Histogram.txt; done
#Parse them to combine splitted lines
cat All_Histogram.txt | awk -F" " 'BEGIN{flag=0} {if(flag==1){if(length($1) > length(old)){print $1""old" "$2" "$3" "$4" "$5" "$6" "$7; flag=0}else{print $1""old" "data; flag=0}}else{if(length($1)<20){flag=1; old=$1; data=($2" "$3" "$4" "$5" "$6" "$7);}else{print $0}}}' > Combined_histograms.tsv



awk -F"\t" '{split($4,info,";"); if(info[2]==1){print $0}}' Hamming_sgRNAs.bed > Hamming_sgRNAs-0MM.bed
awk -F"\t" '{split($4,info,";"); if(info[2]==1){if(info[3]==0){print $0}}}' Hamming_sgRNAs.bed > Hamming_sgRNAs-1MM.bed
awk -F"\t" '{split($4,info,";"); if(info[2]==1){if(info[3]==0){if(info[4]==0){print $0}}}}' Hamming_sgRNAs.bed > Hamming_sgRNAs-2MM.bed
awk -F"\t" '{split($4,info,";"); if(info[2]==1){if(info[3]==0){if(info[4]==0){if(info[5]==0){print $0}}}}}' Hamming_sgRNAs.bed > Hamming_sgRNAs-3MM.bed
awk -F"\t" '{split($4,info,";"); if(info[2]==1){if(info[3]==0){if(info[4]==0){if(info[5]==0){if(info[6]==0){print $0}}}}}}' Hamming_sgRNAs.bed > Hamming_sgRNAs-4MM.bed
awk -F"\t" '{split($4,info,";"); if(info[2]==1){if(info[3]==0){if(info[4]==0){if(info[5]==0){if(info[6]==0){if(info[7]==0){print $0}}}}}}}' Hamming_sgRNAs.bed > Hamming_sgRNAs-5MM.bed



