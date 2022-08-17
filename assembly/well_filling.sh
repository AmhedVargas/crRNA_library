##First attempt
##Working directory
cd ~/Documents/Projects/New_libraries_Oct_2021/CRISPR/2_12_21/warms

cat *.bed | sort -k1,1 -k2,2n > NewLib

perl Process_arms_in_beds.2.pl NewLib pattern > NewLib.modified.arms

awk -F"\t" '{split($4,info,";"); key=$9""info[4]""$10; split(info[length(info)],data,"="); typ=data[2]; if(hash[key]==""){hash[key]=$0; tipos[key]=""}else{tipos[key]=tipos[key]","typ}} END{for (key in hash){OFS="\t"; split(hash[key],coso,"\t"); print coso[1],coso[2],coso[3],coso[4]""tipos[key],coso[5],coso[6],coso[7],coso[8],coso[9],coso[10]}}' NewLib.modified.arms | sort -k1,1 -k2,2n > NewLib.modified.arms.remdup

awk -F"\t" '{OFS="\t"; split($4,info,";"); print $1,$2,$3,$4,$5,$6,$9,$10,"ATTC"substr($10,1,44)"TTCAGActcttgacactggtggccttgcggttagcgcgttcgattagcgtgccggcacccactggttcaattccCTG"substr(info[4],2,19)"GTTtgagagctagacggtttTTT"substr($9,2,44)"GGCT"}' NewLib.modified.arms.remdup | grep -v "TCTAGA" | grep -v "GAGACC" | grep -v "GGTCTC" | awk -F"\t" '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8}' | sort -k1,1 -k2,2n > crRNAlib.v2.rmdup.re.tsv

cd ../assembly/

cp ../warms/crRNAlib.v2.rmdup.re.tsv 

awk -F"\t" '{OFS="\t"; split($4,info,";"); print "GGTCTCtATTC"substr($8,1,44)"TTCTAGActcttgacactggtggccttgcggttagcgcgttcgattagcgtgccggcacccactggttcaattccctG"substr(info[4],2,19)"gtttgagagctagacggtttttt"substr($7,2,44)"GGCTtGAGACC"}' crRNAlib.v2.rmdup.re.tsv > Scafold_to_iterate.txt

R

lib=readLines("Scafold_to_iterate.txt")

fwd=read.table("Ready_primersF.tsv",sep="\t",header=F)

rev=read.table("Ready_primersR.tsv",sep="\t",header=F)

fwd[,2]=as.character(fwd[,2])

fwd[,1]=as.character(fwd[,1])

rev[,1]=as.character(rev[,1])

testis=paste(as.character(c(rep(fwd[1,2],10),rep(fwd[2,2],10),rep(fwd[3,2],10),rep(fwd[4,2],10),rep(fwd[5,2],10),rep(fwd[6,2],10),rep(fwd[7,2],10),rep(fwd[8,2],10),rep(fwd[9,2],10),rep(fwd[10,2],10),rep(fwd[11,2],10),rep(fwd[12,2],10))),paste(lib,rev[,2],sep=" "), sep=" ")

writeLines(testis,"first_attempt.txt")

testis=paste(as.character(c(rep(fwd[1,1],10),rep(fwd[2,1],10),rep(fwd[3,1],10),rep(fwd[4,1],10),rep(fwd[5,1],10),rep(fwd[6,1],10),rep(fwd[7,1],10),rep(fwd[8,1],10),rep(fwd[9,1],10),rep(fwd[10,1],10),rep(fwd[11,1],10),rep(fwd[12,1],10))),paste(lib,rev[,1],sep=" "), sep=" ")

writeLines(testis,"second_attempt.txt")

q()

paste first_attempt.txt second_attempt.txt | awk '{num=((NR-1) % 120)+1; print num"_"$4"_"$6"\tCTCACCGCTCTTGTAGCATG"$1""$2""$3"CCAGGAAGAGATTGCCGGTC"}' > Assembly.tsv

paste crRNAlib.v2.rmdup.re.tsv Assembly.tsv > crRNAlib.v3.tsv

awk -F"\t" 'BEGIN{well =1; plate=1} { print $0"\t"well"\t"plate; if((NR % 120)==0){print $0"\t"well"\t"plate; well++}; if((well % 384)==0){well =1; plate++}}' crRNAlib.v3.tsv | awk -F"\t" '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$12,$11,$9,$10}' > crRNAlib.v3.ordered.tsv

cat crRNAlib.v3.ordered.tsv | awk -F"\t" '{OFS="\t"; split($4,info,";"); type = info[length(info)]; rgb="128,128,128"; if(type =="Type=ATG"){rgb="255,51,51"}; if(type =="Type=ATG_250"){rgb="255,0,0"}; if(type =="Type=ATG_500"){rgb="204,0,0"}; if(type =="Type=CDSguiscan"){rgb="255,128,0"}; if(type =="Type=CDShamm"){rgb="0,128,255"}; print $1,$2,$3,$4,$5,$6,$2,$3,rgb}' > crRNAlib.v3.bed

awk -F"\t" '{split($4,info,";"); print info[length(info)]}' crRNAlib.v3.tsv | sort | uniq -c

###Again and final after removal of ncRNAs
cat *.bed | sort -k1,1 -k2,2n > NewLib
perl Process_arms_in_beds.2.pl NewLib pattern > NewLib.modified.arms
awk -F"\t" '{split($4,info,";"); key=$9""info[4]""$10; split(info[length(info)],data,"="); typ=data[2]; if(hash[key]==""){hash[key]=$0; tipos[key]=""}else{tipos[key]=tipos[key]","typ}} END{for (key in hash){OFS="\t"; split(hash[key],coso,"\t"); print coso[1],coso[2],coso[3],coso[4]""tipos[key],coso[5],coso[6],coso[7],coso[8],coso[9],coso[10]}}' NewLib.modified.arms | sort -k1,1 -k2,2n > NewLib.modified.arms.remdup
awk -F"\t" '{OFS="\t"; split($4,info,";"); print $1,$2,$3,$4,$5,$6,$9,$10,"ATTC"substr($10,1,44)"TTCAGActcttgacactggtggccttgcggttagcgcgttcgattagcgtgccggcacccactggttcaattccCTG"substr(info[4],2,19)"GTTtgagagctagacggtttTTT"substr($9,2,44)"GGCT"}' NewLib.modified.arms.remdup | grep -v "TCTAGA" | grep -v "GAGACC" | grep -v "GGTCTC" | awk -F"\t" '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8}' | sort -k1,1 -k2,2n > crRNAlib.v2.rmdup.re.tsv
cd ../assembly/
cp ../warms/crRNAlib.v2.rmdup.re.tsv .
awk -F"\t" '{OFS="\t"; split($4,info,";"); print "GGTCTCtATTC"substr($8,1,44)"TTCTAGActcttgacactggtggccttgcggttagcgcgttcgattagcgtgccggcacccactggttcaattccctG"substr(info[4],2,19)"gtttgagagctagacggtttttt"substr($7,2,44)"GGCTtGAGACC"}' crRNAlib.v2.rmdup.re.tsv > Scafold_to_iterate.txt

R
lib=readLines("Scafold_to_iterate.txt")
fwd=read.table("Ready_primersF.tsv",sep="\t",header=F)
rev=read.table("Ready_primersR.tsv",sep="\t",header=F)
fwd[,2]=as.character(fwd[,2])
fwd[,1]=as.character(fwd[,1])
rev[,1]=as.character(rev[,1])
testis=paste(as.character(c(rep(fwd[1,2],10),rep(fwd[2,2],10),rep(fwd[3,2],10),rep(fwd[4,2],10),rep(fwd[5,2],10),rep(fwd[6,2],10),rep(fwd[7,2],10),rep(fwd[8,2],10),rep(fwd[9,2],10),rep(fwd[10,2],10),rep(fwd[11,2],10),rep(fwd[12,2],10))),paste(lib,rev[,2],sep=" "), sep=" ")
writeLines(testis,"first_attempt.txt")
testis=paste(as.character(c(rep(fwd[1,1],10),rep(fwd[2,1],10),rep(fwd[3,1],10),rep(fwd[4,1],10),rep(fwd[5,1],10),rep(fwd[6,1],10),rep(fwd[7,1],10),rep(fwd[8,1],10),rep(fwd[9,1],10),rep(fwd[10,1],10),rep(fwd[11,1],10),rep(fwd[12,1],10))),paste(lib,rev[,1],sep=" "), sep=" ")
writeLines(testis,"second_attempt.txt")
q("no")

paste first_attempt.txt second_attempt.txt | awk '{num=((NR-1) % 120)+1; print num"_"$4"_"$6"\tCTCACCGCTCTTGTAGCATG"$1""$2""$3"CCAGGAAGAGATTGCCGGTC"}' > Assembly.tsv
paste crRNAlib.v2.rmdup.re.tsv Assembly.tsv > crRNAlib.v3.tsv
awk -F"\t" 'BEGIN{well =1; plate=1} { print $0"\t"well"\t"plate; if((NR % 120)==0){print $0"\t"well"\t"plate; well++}; if((well % 384)==0){well =1; plate++}}' crRNAlib.v3.tsv | awk -F"\t" '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$12,$11,$9,$10}' > crRNAlib.v3.ordered.tsv
cat crRNAlib.v3.ordered.tsv | awk -F"\t" '{OFS="\t"; split($4,info,";"); type = info[length(info)]; rgb="128,128,128"; if(type =="Type=ATG"){rgb="255,51,51"}; if(type =="Type=ATG_250"){rgb="255,0,0"}; if(type =="Type=ATG_500"){rgb="204,0,0"}; if(type =="Type=CDSguiscan"){rgb="255,128,0"}; if(type =="Type=CDShamm"){rgb="0,128,255"}; print $1,$2,$3,$4,$5,$6,$2,$3,rgb}' > crRNAlib.v3.bed



cat *.bed | sort -k1,1 -k2,2n > NewLib
perl Process_arms_in_beds.2.pl NewLib pattern > NewLib.modified.arms
awk -F"\t" '{split($4,info,";"); key=$9""info[4]""$10; split(info[length(info)],data,"="); typ=data[2]; if(hash[key]==""){hash[key]=$0; tipos[key]=""}else{tipos[key]=tipos[key]","typ}} END{for (key in hash){OFS="\t"; split(hash[key],coso,"\t"); print coso[1],coso[2],coso[3],coso[4]""tipos[key],coso[5],coso[6],coso[7],coso[8],coso[9],coso[10]}}' NewLib.modified.arms | sort -k1,1 -k2,2n > NewLib.modified.arms.remdup
awk -F"\t" '{OFS="\t"; split($4,info,";"); print $1,$2,$3,$4,$5,$6,$9,$10,"ATTC"substr($10,1,44)"TTCAGActcttgacactggtggccttgcggttagcgcgttcgattagcgtgccggcacccactggttcaattccCTG"substr(info[4],2,19)"GTTtgagagctagacggtttTTT"substr($9,2,44)"GGCT"}' NewLib.modified.arms.remdup | grep -v "TCTAGA" | grep -v "GAGACC" | grep -v "GGTCTC" | awk -F"\t" '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8}' | sort -k1,1 -k2,2n > crRNAlib.v2.rmdup.re.tsv
cd ../assembly/
cp ../warms/crRNAlib.v2.rmdup.re.tsv .
awk -F"\t" '{OFS="\t"; split($4,info,";"); print "GGTCTCtATTC"substr($8,1,44)"TTCTAGActcttgacactggtggccttgcggttagcgcgttcgattagcgtgccggcacccactggttcaattccctG"substr(info[4],2,19)"gtttgagagctagacggtttttt"substr($7,2,44)"GGCTtGAGACC"}' crRNAlib.v2.rmdup.re.tsv > Scafold_to_iterate.txt




R
lib=readLines("Scafold_to_iterate.txt")
fwd=read.table("Ready_primersF.tsv",sep="\t",header=F)
rev=read.table("Ready_primersR.tsv",sep="\t",header=F)
fwd[,2]=as.character(fwd[,2])
fwd[,1]=as.character(fwd[,1])
rev[,1]=as.character(rev[,1])

testis=paste(as.character(c(rep(fwd[1,2],10),rep(fwd[2,2],10),rep(fwd[3,2],10),rep(fwd[4,2],10),rep(fwd[5,2],10),rep(fwd[6,2],10),rep(fwd[7,2],10),rep(fwd[8,2],10),rep(fwd[9,2],10),rep(fwd[10,2],10),rep(fwd[11,2],10),rep(fwd[12,2],10))),paste(lib,rev[,2],sep=" "), sep=" ")
writeLines(testis,"first_attempt.txt")

testis=paste(as.character(c(rep(fwd[1,1],10),rep(fwd[2,1],10),rep(fwd[3,1],10),rep(fwd[4,1],10),rep(fwd[5,1],10),rep(fwd[6,1],10),rep(fwd[7,1],10),rep(fwd[8,1],10),rep(fwd[9,1],10),rep(fwd[10,1],10),rep(fwd[11,1],10),rep(fwd[12,1],10))),paste(lib,rev[,1],sep=" "), sep=" ")
writeLines(testis,"second_attempt.txt")
q()

paste first_attempt.txt second_attempt.txt | awk '{num=((NR-1) % 120)+1; print num"_"$4"_"$6"\tCTCACCGCTCTTGTAGCATG"$1""$2""$3"CCAGGAAGAGATTGCCGGTC"}' > Assembly.tsv

paste crRNAlib.v2.rmdup.re.tsv Assembly.tsv > crRNAlib.v3.tsv

awk -F"\t" 'BEGIN{well =1; plate=1} { print $0"\t"well"\t"plate; if((NR % 120)==0){print $0"\t"well"\t"plate; well++}; if((well % 384)==0){well =1; plate++}}' crRNAlib.v3.tsv | awk -F"\t" '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$12,$11,$9,$10}' > crRNAlib.v3.ordered.tsv

cat crRNAlib.v3.ordered.tsv | awk -F"\t" '{OFS="\t"; split($4,info,";"); type = info[length(info)]; rgb="128,128,128"; if(type =="Type=ATG"){rgb="255,51,51"}; if(type =="Type=ATG_250"){rgb="255,0,0"}; if(type =="Type=ATG_500"){rgb="204,0,0"}; if(type =="Type=CDSguiscan"){rgb="255,128,0"}; if(type =="Type=CDShamm"){rgb="0,128,255"}; print $1,$2,$3,$4,$5,$6,$2,$3,rgb}' > crRNAlib.v3.bed

wc -l crRNAlib.v3.bed

awk -F"\t" '{split($4,info,";"); print info[length(info)]}' crRNAlib.v3.tsv | sort | uniq -c




