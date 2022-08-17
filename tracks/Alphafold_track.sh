###With regards to Alphafold conversion
##Working directory
cd /home/velazqam/Documents/Projects/New_libraries_Oct_2021/Alphafold
##Make directory to store pDDTT errors
mkdir pDDTT
##Obtain pDDTT values for all AA
for file in `ls ../UP000001940_6239_CAEEL/*.cif.gz`; do zcat ../UP000001940_6239_CAEEL/${file} | perl -pe 's/\s+/\t/g' | perl -pe 's/ATOM/\nATOM/g' | awk -F"\t" '{print $9";"$25";"$6";"$15}' | uniq | tail -n +2 | perl -pe 's/\;/\t/g' > ${file%-model_v1.cif.gz}.tsv; done
##Move the errors two folders above
mv ../UP000001940_6239_CAEEL/*.tsv .
##Generate a list of them
ls *.tsv > files.txt
##Move to directory from where analisis will be performed
cd /home/velazqam/Documents/Projects/New_libraries_Oct_2021/CRISPR/AA_coordinates/
##Copy pDDTT errors
cp -r /home/velazqam/Documents/Projects/New_libraries_Oct_2021/Alphafold/pDDTT/ .
##Basic analisis with R
#Move to directory
cd pDDTT/
#Create a list of files to be analized
#ls *.tsv > files.txt
##Now R
R
library(IRanges)
files=readLines("files.txt")
vals=c()
for(file in files){tt=read.table(file,sep="\t",header=F); vals=append(vals,median(tt[,4]))}
pdf("pDDTT_distribution.pdf")
boxplot(vals,main="Median pDDTT value across C. elegans proteins", yaxt="n", xlab="pDDTT scores")
axis(2,at=(c(2:10)*10),labels=T)
dev.off()

pass=c()
for(file in files){
    tt=read.table(file,sep="\t",header=F); 
    tresh=median(tt[,4])
    if(tresh > 70){
        pass=append(pass,file)
        nfile=sub(".tsv","_loc.csv",file)
        up=70
        down=70
        st=min(which(tt[,4]>up))
        ed=max(which(tt[,4]>up))
        ##Define ranges
        rg=c()
        rg=which((tt[,4]< down) & (tt[,1]>=st & tt[,1]<=ed))
        if(length(rg)>0){
            a=IRanges(rg)
            write.table(x=(as.data.frame(a)[,1:2]),nfile,sep=",",col.names=F,row.names=F,quote=F)
        }
    }
    }
q()
##NOw let's do conversions
cd ..
#go to loc
mkdir loc
cd loc
#Copy files
cp ../pDDTT/*.csv .
##How many proteins
ls *.csv | wc -l
#12177
#copy proper file
cp ../Analisis/Cel_GeneswPROT.bed .
##Make exons
awk -F"\t" '{OFS="\t"; split($7,starts,";"); split($8,ends,";"); strand=$6; for(i=1; i<= length(starts); i++){if(strand=="+"){num=i}else{num=(length(starts)-i+1)}; print $1,starts[i],ends[i],$4";Exon;"num";",$5,strand}}' Cel_GeneswPROT.bed > Cel_GeneswPROT_exons.bed
##Now introns
awk -F"\t" '{OFS="\t"; split($7,starts,";"); split($8,ends,";"); strand=$6; for(i=1; i<= (length(starts)-1); i++){if(strand=="+"){num=i}else{num=(length(starts)-i)}; print $1,ends[i],starts[i+1],$4";Intron;"num";",$5,strand}}' Cel_GeneswPROT.bed  > Cel_GeneswPROT_introns.bed
##MAke genome
zcat ../Analisis/c_elegans.PRJNA13758.WS282.genomic.fa.gz > c_elegans.genome.fa
#make full data
cat Cel_GeneswPROT_exons.bed Cel_GeneswPROT_introns.bed | sort -k1,1 -k2,2n > full_data.bed
##NOw combinethem
bedtools getfasta -fi c_elegans.genome.fa -bed full_data.bed -name -tab -s | awk -F"\t" '{split($1,info,";"); print info[(length(info)-1)]"\t"info[(length(info)-2)]"\t"$0}' | sort -k1,1n -k2,2 -k3,3 | awk -F"\t" '{split($3,info,";"); key=info[1]";"info[2]";"info[3]; seq=$4; if($2=="Intron"){seq=tolower(seq)}; hash[key]=hash[key]""seq;} END{for(key in hash){print key"\t"hash[key]}}'  > Transcripts_seqs.tsv
###Then add names and coordinates
awk -F"\t" '{if(NF==2){seq[$1]=$2}else{split($4,info,";"); key=info[1]";"info[2]";"info[3]; print $0"\t"seq[key]}}' Transcripts_seqs.tsv Cel_GeneswPROT.bed > FullDB.tsv
##THen make DB
mkdir DB
cd DB
##Make files similar to "AF-X5LX93-F1_loc.csv" where the uniprot name is inbetween
awk -F"\t" '{split($4,info,";"); print $0 > "AF-"info[4]"-F1_seq.txt"}' ../FullDB.tsv
###Proof of concept
awk -F"\t" 'BEGIN{count=0}{if(NF==1){count = (count + 1); split($0,coors,","); start[count]=coors[1]; end[count]=coors[2]}else{if(count > 0){strand=$6; split($9,seq,""); for(i=1; i<= count; i++){aastart=(start[i]*3-3); aaend=(end[i]*3); ggstart=0; ggend=0; pos=0; for(j=1; j<=length(seq); j++){if(seq[j] ~ /[A-Z]/){pos=pos+1}; if(pos==aastart){ggstart=j}; if(pos==aaend){ggend=j} }; if(strand =="+"){print $1"\t"($2+ggstart)"\t"($2+ggend)"\t"$4"\t"$5"\t"strand}else{print $1"\t"($3-ggend)"\t"($3-ggstart)"\t"$4"\t"$5"\t"strand}} }}}' ../loc/AF-A0A060Q5Z9-F1_loc.csv  AF-A0A060Q5Z9-F1_seq.txt
##csr-1
awk -F"\t" 'BEGIN{count=0}{if(NF==1){count = (count + 1); split($0,coors,","); start[count]=coors[1]; end[count]=coors[2]}else{if(count > 0){strand=$6; split($9,seq,""); for(i=1; i<= count; i++){aastart=(start[i]*3-3); aaend=(end[i]*3); ggstart=0; ggend=0; pos=0; for(j=1; j<=length(seq); j++){if(seq[j] ~ /[A-Z]/){pos=pos+1}; if(pos==aastart){ggstart=j}; if(pos==aaend){ggend=j} }; if(strand =="+"){print $1"\t"($2+ggstart)"\t"($2+ggend)"\t"$4"\t"$5"\t"strand}else{print $1"\t"($3-ggend)"\t"($3-ggstart)"\t"$4"\t"$5"\t"strand}} }}}' ../loc/AF-H2KZD5-F1_loc.csv  AF-H2KZD5-F1_seq.txt
##Now big for loop
for file in `ls *.txt`; do awk -F"\t" 'BEGIN{count=0}{if(NF==1){count = (count + 1); split($0,coors,","); start[count]=coors[1]; end[count]=coors[2]}else{if(count > 0){strand=$6; split($9,seq,""); for(i=1; i<= count; i++){aastart=(start[i]*3-3); aaend=(end[i]*3); ggstart=0; ggend=0; pos=0; for(j=1; j<=length(seq); j++){if(seq[j] ~ /[A-Z]/){pos=pos+1}; if(pos==aastart){ggstart=j}; if(pos==aaend){ggend=j} }; if(strand =="+"){print $1"\t"($2+ggstart)"\t"($2+ggend)"\t"$4"\t"$5"\t"strand}else{print $1"\t"($3-ggend)"\t"($3-ggstart)"\t"$4"\t"$5"\t"strand}} }}}' ../loc/${file%_seq.txt}_loc.csv ${file} > ${file%_seq.txt}_regions.bed; done
##Now move to another directory
cd ..
mkdir final_intersection
cd final_intersection
cat ../DB/*.bed | sort -k1,1 -k2,2n > AA_linkers_alphafold.bed
###Didnt work
##redo by renaming files
cd ..
mkdir REDO
cd REDO
mkdir locs
mkdir trans
cd locs
cd ..
echo ""> coms.sh; for file in `ls ../pDDTT/*.csv `; do base=(${file#../pDDTT/}); num=(`wc -l ${file%_loc.csv}.tsv`); echo "cp "${file}" "${base%-F1_loc.csv}"si"${num}"ze-F1_loc.csv" >> coms.sh;  done
sh coms.sh
mv *.csv locs
cd trans
awk -F"\t" '{split($4,info,";"); count=0; split($9,seq,""); for(i=1; i<= length(seq); i++){if(seq[i] ~ /[A-Z]/){count++}}; nono=((count/3)-1); print $0 > "AF-"info[4]"si"nono"ze-F1_seq.txt"}' ../../FullDB.tsv
##Now big for loop
for file in `ls *.txt`; do awk -F"\t" 'BEGIN{count=0}{if(NF==1){count = (count + 1); split($0,coors,","); start[count]=coors[1]; end[count]=coors[2]}else{if(count > 0){strand=$6; split($9,seq,""); for(i=1; i<= count; i++){aastart=(start[i]*3-3); aaend=(end[i]*3); ggstart=0; ggend=0; pos=0; for(j=1; j<=length(seq); j++){if(seq[j] ~ /[A-Z]/){pos=pos+1}; if(pos==aastart){ggstart=j}; if(pos==aaend){ggend=j} }; if(strand =="+"){print $1"\t"($2+ggstart)"\t"($2+ggend)"\t"$4"\t"$5"\t"strand}else{print $1"\t"($3-ggend)"\t"($3-ggstart)"\t"$4"\t"$5"\t"strand}} }}}' ../locs/${file%_seq.txt}_loc.csv ${file} > ${file%_seq.txt}_regions.bed; done
cd ..
mkdir final_intersection
cd final_intersection
cat ../trans/*.bed | sort -k1,1 -k2,2n > AA_linkers_alphafold.bed

bedtools subtract -a AA_linkers_alphafold.bed -b ../../Cel_GeneswPROT_introns.bed -s | head



echo ""> coms.sh; for file in `ls ../../pDDTT/*.tsv `; do base=(${file#../../pDDTT/}); num=(`wc -l ${file}`); echo "cp "${file}" "${base%-F1.tsv}"si"${num}"ze-F1.tsv" >> coms.sh;  done


ls *.tsv > files.txt
##Now R
R
library(IRanges)
files=readLines("files.txt")

pass=c()
for(file in files){
    tt=read.table(file,sep="\t",header=F); 
    tresh=median(tt[,4])
    if(tresh > 70){
        pass=append(pass,file)
        nfile=sub(".tsv","_75_loc.csv",file)
        up=75
        down=75
        st=min(which(tt[,4]>up))
        ed=max(which(tt[,4]>up))
        ##Define ranges
        rg=c()
        rg=which((tt[,4]< down) & (tt[,1]>=st & tt[,1]<=ed))
        if(length(rg)>0){
            a=IRanges(rg)
            write.table(x=(as.data.frame(a)[,1:2]),nfile,sep=",",col.names=F,row.names=F,quote=F)
        }
    }
    }

pass=c()
for(file in files){
    tt=read.table(file,sep="\t",header=F); 
    tresh=median(tt[,4])
    if(tresh > 70){
        pass=append(pass,file)
        nfile=sub(".tsv","_80_loc.csv",file)
        up=80
        down=80
        st=min(which(tt[,4]>up))
        ed=max(which(tt[,4]>up))
        ##Define ranges
        rg=c()
        rg=which((tt[,4]< down) & (tt[,1]>=st & tt[,1]<=ed))
        if(length(rg)>0){
            a=IRanges(rg)
            write.table(x=(as.data.frame(a)[,1:2]),nfile,sep=",",col.names=F,row.names=F,quote=F)
        }
    }
    }

pass=c()
for(file in files){
    tt=read.table(file,sep="\t",header=F); 
    tresh=median(tt[,4])
    if(tresh > 70){
        pass=append(pass,file)
        nfile=sub(".tsv","_85_loc.csv",file)
        up=85
        down=85
        st=min(which(tt[,4]>up))
        ed=max(which(tt[,4]>up))
        ##Define ranges
        rg=c()
        rg=which((tt[,4]< down) & (tt[,1]>=st & tt[,1]<=ed))
        if(length(rg)>0){
            a=IRanges(rg)
            write.table(x=(as.data.frame(a)[,1:2]),nfile,sep=",",col.names=F,row.names=F,quote=F)
        }
    }
    }

pass=c()
for(file in files){
    tt=read.table(file,sep="\t",header=F); 
    tresh=median(tt[,4])
    if(tresh > 70){
        pass=append(pass,file)
        nfile=sub(".tsv","_90_loc.csv",file)
        up=90
        down=90
        st=min(which(tt[,4]>up))
        ed=max(which(tt[,4]>up))
        ##Define ranges
        rg=c()
        rg=which((tt[,4]< down) & (tt[,1]>=st & tt[,1]<=ed))
        if(length(rg)>0){
            a=IRanges(rg)
            write.table(x=(as.data.frame(a)[,1:2]),nfile,sep=",",col.names=F,row.names=F,quote=F)
        }
    }
    }

q()


for file in `ls *.txt`; do awk -F"\t" 'BEGIN{count=0}{if(NF==1){count = (count + 1); split($0,coors,","); start[count]=coors[1]; end[count]=coors[2]}else{if(count > 0){strand=$6; split($9,seq,""); for(i=1; i<= count; i++){aastart=(start[i]*3-3); aaend=(end[i]*3); ggstart=0; ggend=0; pos=0; for(j=1; j<=length(seq); j++){if(seq[j] ~ /[A-Z]/){pos=pos+1}; if(pos==aastart){ggstart=j}; if(pos==aaend){ggend=j} }; if(strand =="+"){print $1"\t"($2+ggstart)"\t"($2+ggend)"\t"$4"\t"$5"\t"strand}else{print $1"\t"($3-ggend)"\t"($3-ggstart)"\t"$4"\t"$5"\t"strand}} }}}' ../pDDTT/${file%_seq.txt}_75_loc.csv ${file} > ../75/${file%_seq.txt}_regions.bed; done

for file in `ls *.txt`; do awk -F"\t" 'BEGIN{count=0}{if(NF==1){count = (count + 1); split($0,coors,","); start[count]=coors[1]; end[count]=coors[2]}else{if(count > 0){strand=$6; split($9,seq,""); for(i=1; i<= count; i++){aastart=(start[i]*3-3); aaend=(end[i]*3); ggstart=0; ggend=0; pos=0; for(j=1; j<=length(seq); j++){if(seq[j] ~ /[A-Z]/){pos=pos+1}; if(pos==aastart){ggstart=j}; if(pos==aaend){ggend=j} }; if(strand =="+"){print $1"\t"($2+ggstart)"\t"($2+ggend)"\t"$4"\t"$5"\t"strand}else{print $1"\t"($3-ggend)"\t"($3-ggstart)"\t"$4"\t"$5"\t"strand}} }}}' ../pDDTT/${file%_seq.txt}_80_loc.csv ${file} > ../80/${file%_seq.txt}_regions.bed; done

for file in `ls *.txt`; do awk -F"\t" 'BEGIN{count=0}{if(NF==1){count = (count + 1); split($0,coors,","); start[count]=coors[1]; end[count]=coors[2]}else{if(count > 0){strand=$6; split($9,seq,""); for(i=1; i<= count; i++){aastart=(start[i]*3-3); aaend=(end[i]*3); ggstart=0; ggend=0; pos=0; for(j=1; j<=length(seq); j++){if(seq[j] ~ /[A-Z]/){pos=pos+1}; if(pos==aastart){ggstart=j}; if(pos==aaend){ggend=j} }; if(strand =="+"){print $1"\t"($2+ggstart)"\t"($2+ggend)"\t"$4"\t"$5"\t"strand}else{print $1"\t"($3-ggend)"\t"($3-ggstart)"\t"$4"\t"$5"\t"strand}} }}}' ../pDDTT/${file%_seq.txt}_85_loc.csv ${file} > ../85/${file%_seq.txt}_regions.bed; done

for file in `ls *.txt`; do awk -F"\t" 'BEGIN{count=0}{if(NF==1){count = (count + 1); split($0,coors,","); start[count]=coors[1]; end[count]=coors[2]}else{if(count > 0){strand=$6; split($9,seq,""); for(i=1; i<= count; i++){aastart=(start[i]*3-3); aaend=(end[i]*3); ggstart=0; ggend=0; pos=0; for(j=1; j<=length(seq); j++){if(seq[j] ~ /[A-Z]/){pos=pos+1}; if(pos==aastart){ggstart=j}; if(pos==aaend){ggend=j} }; if(strand =="+"){print $1"\t"($2+ggstart)"\t"($2+ggend)"\t"$4"\t"$5"\t"strand}else{print $1"\t"($3-ggend)"\t"($3-ggstart)"\t"$4"\t"$5"\t"strand}} }}}' ../pDDTT/${file%_seq.txt}_90_loc.csv ${file} > ../90/${file%_seq.txt}_regions.bed; done

##It worked despite my wrong codings
##Now do finalize by taking out intronics sequences

cd /home/velazqam/Documents/Projects/New_libraries_Oct_2021/CRISPR/AA_coordinates/Dif_pDDTT

cat 75/*.bed | sort -k1,1 -k2,2n > Alphfold_linkers.75.bed
cat 80/*.bed | sort -k1,1 -k2,2n > Alphfold_linkers.80.bed
cat 85/*.bed | sort -k1,1 -k2,2n > Alphfold_linkers.85.bed
cat 90/*.bed | sort -k1,1 -k2,2n > Alphfold_linkers.90.bed

##Substract 
bedtools subtract -a Alphfold_linkers.75.bed -b ../Cel_GeneswPROT_introns.bed -s > Alphafold_linkers.75.CDS.bed
bedtools subtract -a Alphfold_linkers.80.bed -b ../Cel_GeneswPROT_introns.bed -s > Alphafold_linkers.80.CDS.bed
bedtools subtract -a Alphfold_linkers.85.bed -b ../Cel_GeneswPROT_introns.bed -s > Alphafold_linkers.85.CDS.bed
bedtools subtract -a Alphfold_linkers.90.bed -b ../Cel_GeneswPROT_introns.bed -s > Alphafold_linkers.90.CDS.bed

