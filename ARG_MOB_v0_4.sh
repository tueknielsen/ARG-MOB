#ARG-MOB version 0.4
#Tue Kjærgaard Nielsen
#January 7th 2021

#Disclaimer:
#This bash script performs identification of ARGs and MGEs in RefSeq complete bacterial genomes and prepared for further analyses in R.
#It is provided as documentation but is not written with the intend of enabling other researchers to run it as a script or software. \
As such, it is not cleaned up or streamlined for wide distribution and implementation.
#The script is not meant to be run as a program. There are a few steps where it is important to wait for multiple parallel opearations to finish before the workflow is continued. \
It is therefore advised to run each step independently and monitor when it has finished, before moving on. 

#Many tasks are parallelized via functions that are set to run 100 tasks in parallel. This may not be optimal for all computers/servers. \
Analyses were performed on our server "thoth" with the following specifications:
#HPE ProLiant DL580 Gen10 Server equipped with four Intel Xeon Gold 6254 processors and 3 TB RAM running Ubuntu 18.04 LTS

#During development, the pipeline was named 'ProxARG' and this name may still occur in the script.

#Important note: some ARGs and IS elements are encoded by two or more ORFs. These will be overrepresented in final results, but should not overall affect results much.

#Non-exhaustive list of applied software and their versions
#DIAMOND 				v0.9.25.126
#ncbi-genome-download 	0.2.11
#FASTX toolkit 			0.0.14
#Prodigal 				V2.6.3
#Samtools 				1.10
#seqtk 					1.3-r106
#USEARCH 				v11.0.667_i86linux64
#Integron_finder 		1.5.1
#barrnap				0.9



#Databases required					path on our thoth server
#RefSeq complete genomes			/database/ncbi/ncbi_complete_genomes_121219/refseq_complete.fasta
#CARD								/home/tkn/work/ProxARG19/nucleotide_fasta_protein_homolog_model
#CARD ARO index						/DATA_1/tkn/ProxARG19/CARD_db/aro_index.tsv
#RGI 								/DATA_1/tkn/ProxARG19/localDB/proteindb.fsa
#NCBI acc to taxid					/DATA_1/tkn/ProxARG19/nucl_gb.accession2taxid.subset
#NCBI fullnamelineage				/DATA_1/tkn/ProxARG19/fullnamelineage_edit.dmp
##ISfinder db as included in prokka /opt/prokka/db/kingdom/Bacteria/IS

#Define the number of threads used for parallel processing
numthreads=100

#Paths to programs - not all are listed here
diamond=/home/tkn/Programs/diamond
usearch=usearch_64

### Database notes
#For faster grep'ing in the nucl_gb.accession2taxid file, subset to only accession numbers starting with NZ_ and NC_ (complete genomes)
#rm -f nucl_gb.accession2taxid.subset
#grep '^NC_' nucl_gb.accession2taxid > nucl_gb.accession2taxid.subset
#grep '^NZ_' nucl_gb.accession2taxid >> nucl_gb.accession2taxid.subset
#A few entries are not correctly formatted. Correct these in the subset
#grep 'CP035045' nucl_gb.accession2taxid | sed 's/CP/NZ_CP/g' >> nucl_gb.accession2taxid.subset
#grep 'CP035049' nucl_gb.accession2taxid | sed 's/CP/NZ_CP/g' >> nucl_gb.accession2taxid.subset
#grep 'CP035051' nucl_gb.accession2taxid | sed 's/CP/NZ_CP/g' >> nucl_gb.accession2taxid.subset
#grep 'CP035047' nucl_gb.accession2taxid | sed 's/CP/NZ_CP/g' >> nucl_gb.accession2taxid.subset
#grep 'CP035052' nucl_gb.accession2taxid | sed 's/CP/NZ_CP/g' >> nucl_gb.accession2taxid.subset


#Superphyla are confusing the taxonomy order, since they push the phylum one level down for bacteria with superphyla assigned.
#Remove supergroups from fullnamelineage.dmp
#sed -e 's/Terrabacteria group; //' -e 's/PVC group; //' \
-e 's/FCB group; Bacteroidetes\/Chlorobi group; //' \
-e 's/FCB group; //' -e 's/Cyanobacteria\/Melainabacteria group; //' fullnamelineage.dmp > fullnamelineage_edit.dmp

#All complete refseq bacterial genomes were downloaded from NCBI with
#ncbi-genome-download -l complete -p 20 -F fasta bacteria
#Path: /database/ncbi/ncbi_completef_genomes_121219/refseq/bacteria/
#Database: /database/ncbi/ncbi_complete_genomes_121219/refseq_complete
#CARD database: /home/tkn/work/ProxARG19/nucleotide_fasta_protein_homolog_model


#Download all gff annoations for extracting 16S control genes
#ncbi-genome-download -l complete -p 20 -F gff bacteria -o complete_gff -v
#ncbi-genome-download -l complete -p 20 -F faa bacteria -o complete_faa -v




###Start of analysis

#Blast all genomes against CARD. Proteins from RefSeq complete genomes were previously predicted with prodigal to ensure that gene calling was performed uniformly across all genomes.
#Prodigal loop to speed up gene calling on all complete genomes
#fasta_formatter -i refseq_complete.fasta -t -o refseq_complete.tab
#split -l 400 refseq_complete.tab
#for i in x??
#do
#awk -F'\t' '{print ">"$1"\n"$2}' $i > $i.fa
#prodigal -a $i.faa -d $i.fna -p meta -i $i.fa &
#done
#cat *faa > all_refseq_complete_prodigal.faa
#Format the faa fasta headers, so they include the gene positions. Also add "pr" before each protein number
#fasta_formatter -i all_refseq_complete_prodigal.faa -w0 | sed -e 's/_/_pr/2' -e 's/ \# /_/' -e 's/ \# /_/' > all_refseq_complete_prodigal_edited.faa
#sed 's/ # .*//' /DATA_1/tkn/ProxARG19/all_refseq_complete_prodigal_edited.faa > /DATA_1/tkn/ProxARG19/all_refseq_complete_prodigal_edited2.faa

$diamond blastp --query /DATA_1/tkn/temp_prodigal/all_refseq_complete_prodigal_single_meta.faa --threads $numthreads --db /DATA_1/tkn/ProxARG19/CARD_db/protein_fasta_protein_homolog_model.dmnd \
--max-target-seqs 1 --evalue 10E-10 --query-cover 80 --subject-cover 80 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp \
--out CARD_refseq.diamond --header

#wc -l CARD_refseq.diamond
#1341463 CARD_refseq.diamond 
# This is also counting the header lines that happened because of the '--header' flag


#Wait for diamond blastp to finish

#From CARD/RGI there are curated bitscore cutoffs for each database protein (RGI cutoffs). Extract these and filter hits.
grep -v "^#" CARD_refseq.diamond > temp.dmnd
rm -f x??
split --number=l/$numthreads -d temp.dmnd
rgi_cutoff_fun () {
    curr_sub=$1
    number=$2
    acc_list=$(grep -v "#" $curr_sub | awk '{print $2}' | sed -e 's/gb|//' -e 's/|.*//')
    rm -f rgi_cutoffs_$number
    for acc in $acc_list
    do
        model=$(grep "$acc" /DATA_1/tkn/ProxARG19/CARD_db/aro_index.tsv | awk '{print $4"_"$3}')
        #echo $model
        bitscore_cutoff=$(grep "$model" /DATA_1/tkn/ProxARG19/localDB/proteindb.fsa | awk '{print $7}')
        #echo $bitscore_cutoff
        echo $bitscore_cutoff >> rgi_cutoffs_$number
        done
}

rm -f pids
for a in x??
do
    number=$(echo $a | sed 's/x//')
    rgi_cutoff_fun $a $number &
    echo $! >> pids
done

#Wait for all processes to finish
all_pids=$(cat pids)
for pid in $all_pids; do
echo $pid
wait $pid
done

#Add rgi cutoffs back to diamond hits
for b in x??
do
    number=$(echo $b | sed 's/x//')
    paste $b rgi_cutoffs_$number > temp.$number
done

cat temp.[0-9][0-9] > diamond_hits_rgi_cutoff
rm temp.* x?? rgi_cutoffs_??
#wc -l diamond_hits_rgi_cutoff
#1341460 diamond_hits_rgi_cutoff
# A difference of three because of the headers!

#If the bitscore is higher than the curated RGI-defined cutoff OR if %ID is higher than 80%, print the diamond hits. Add the ratio of bitscore/cutoff
awk '{if ($12 >= $14 || $3 >= 80) print $0"\t"$12/$14}' diamond_hits_rgi_cutoff  > diamond_hits_rgi_cutoff.pass
awk '{if ($12 < $14 && $3 < 80) print $0"\t"$12/$14}' diamond_hits_rgi_cutoff  > diamond_hits_rgi_cutoff.fail

#176888 diamond_hits_rgi_cutoff.pass
#1164572 diamond_hits_rgi_cutoff.fail

cp diamond_hits_rgi_cutoff.pass CARD_refseq_filt.diamond


#Investigate the taxonomy of RGI score cutoff and %ID > 80% cutoff
awk '{if ($12 < $14 && $3 >= 80) print $0"\t"$12/$14}' diamond_hits_rgi_cutoff  > diamond_hits_rgi_fail_ID_pass
accs=$(awk -F'_' '{print $1"_"$2}' diamond_hits_rgi_fail_ID_pass)

tax_function () {
    acc_list=$1
    number=$2
    accessions=$(cat $acc_list)
    rm -f tax.$number
    for acc in $accessions
    do
        taxid=$(grep -m 1 -w "$acc" /DATA_1/tkn/ProxARG19/nucl_gb.accession2taxid.subset | awk '{print $3}' )
        taxonomy=$(grep -m 1 -w "^$taxid" /DATA_1/tkn/ProxARG19/fullnamelineage_edit.dmp | sed -e 's/ /_/g' | awk -F'|' '{print $3}' | awk -F';_' '{print $1";_"$2";_"$3";_"$4";_"$5";_"$6";_"$7}')
        echo $taxonomy | sed -e 's/ /\t/g' -e 's/;_/\t/g' >> tax.$number
    done
}

#Process diamond hits that fail RGI cutoff but pass min. 80% ID (and 80% query coverage)
mkdir tempdir
cd tempdir
awk -F'_' '{print $1"_"$2}' ../diamond_hits_rgi_fail_ID_pass > all_accessions
split --number=l/$numthreads -d all_accessions

#Start the function in a loop
rm -f pids
for i in x??
do
    num=$(echo $i | sed 's/x//')
    echo $i $num
    tax_function $i $num &
    echo $! >> pids
done

#Wait for all processes to finish
all_pids=$(cat pids)
for pid in $all_pids; do
    echo $pid
    wait $pid
done

cd ..
cat tempdir/tax.?? > rgi_fail_ID_pass.tax
rm -r tempdir
#wc -l rgi_fail_ID_pass.tax
#61620 rgi_fail_ID_pass.tax

##Process diamond hits that pass RGI cutoff
awk '{if ($12 >= $14) print $0"\t"$12/$14}' diamond_hits_rgi_cutoff  > diamond_hits_rgi_pass
#115268 diamond_hits_rgi_pass

mkdir tempdir
cd tempdir
awk -F'_' '{print $1"_"$2}' ../diamond_hits_rgi_pass > all_accessions
split --number=l/$numthreads -d all_accessions

rm -f pids
for i in x??
do
    num=$(echo $i | sed 's/x//')
    tax_function $i $num &
    echo $! >> pids
done

#Wait for all processes to finish
all_pids=$(cat pids)
for pid in $all_pids; do
    echo $pid
    wait $pid
done

cd ..
cat tempdir/tax.?? > rgi_pass.tax
rm -r tempdir
#115268 rgi_pass.tax

paste diamond_hits_rgi_pass rgi_pass.tax > diamond_hits_rgi_pass.Rdat
paste diamond_hits_rgi_fail_ID_pass rgi_fail_ID_pass.tax > diamond_hits_rgi_fail_ID_pass.Rdat

#Hit regions might extend beyond the end, but only until the end is extracted. The end coordinate will still be wrong (larger than the total length). 
#If start coordinate is negative, replace it with 1
halfdist=12170

awk '{print $1}' CARD_refseq_filt.diamond | sed 's/>//' | awk -F'_' -v var="$halfdist"  '{print $1"_"$2" "$3-var" "$4+var}' | awk '{if ($2 < 0) print $1":1-"$3; else print $1":"$2"-"$3}' | sed '/\#/d' > ARG_coordinates.txt
awk '{print $2}' CARD_refseq_filt.diamond | awk -F'\|' '{print $3}' | sed -e 's/ARO://' -e '/^$/d' > AROS
awk '{print $1}' CARD_refseq_filt.diamond | awk -F'_' '{print $5}' > prID
wc -l ARG_coordinates.txt AROS prID
#176888 ARG_coordinates.txt
#176888 AROS
#176888 prID

xargs samtools faidx /database/ncbi/ncbi_complete_genomes_121219/refseq_complete.fasta < ARG_coordinates.txt > temp0
grep -c '>' temp0
#176888

#Fasta to tab to append the AROs annotations
fasta_formatter -i temp0 -t -o temp
paste temp prID AROS | awk '{print ">"$1"_"$3"|"$4"\n"$2}' > ARG_regions.fasta
grep -c '>' ARG_regions.fasta
#176888

rm temp* prID

####

#Extract loci based on ARG hits and extract subset of prodigal-predicted amino acid sequences for ARG loci passing passing filters.
grep '>' ARG_regions.fasta | sed -e 's/>//' -e 's/:/_/' -e 's/-/_/' -e 's/|/_/' > temp.coord
rm -f c???
split --number=l/$numthreads -a 3 -d temp.coord c


round5k() {
    echo $(( ((${1%.*}+5000)/5000)*5000 ))
}

extract_fun () {
	file_in=$1
    num=$2
    coord_list=$(cat $file_in)
    for i in $coord_list
    do
        acc=$(echo $i | awk -F'_' '{print $1"_"$2}')
        start_coord=$(echo $i | awk -F'_' '{print $3}')
        end_coord=$(echo $i | awk -F'_' '{print $4}')
        prefix=$(echo $i | awk -F'_' '{print $1}')
        if [ $prefix == "NC" ]
        then
            grep "$acc"  /DATA_1/tkn/temp_prodigal/NC.faa.tab | sed 's/ /_/g' | \
            awk -F'_' -v start=$start_coord -v end=$end_coord -v locus=$i '{if ($3 >= start && $4 <= end) print locus":"$0}' >> temp.$num
        else
            prefix2=$(echo $i | awk -F'_' '{print $2}' | sed 's/[0-9]*\.[0-9]*//')
        if [ $prefix2 == "CP" ]
        then
            number=$(echo $i | awk -F'_' '{print $2}' | sed -e 's/CP//' -e 's/\.[0-9]*//' -e 's/^0*//')
            rnum=$(round5k $number)
            grep "$acc"  /DATA_1/tkn/temp_prodigal/NZ_CP.$rnum.tab | sed 's/ /_/g' | \
            awk -F'_' -v start=$start_coord -v end=$end_coord -v locus=$i '{if ($3 >= start && $4 <= end) print locus":"$0}' >> temp.$num
        else
            grep "$acc"  /DATA_1/tkn/temp_prodigal/NZ_$prefix2.faa.tab | sed 's/ /_/g' | \
            awk -F'_' -v start=$start_coord -v end=$end_coord -v locus=$i '{if ($3 >= start && $4 <= end) print locus":"$0}' >> temp.$num
        fi
        fi
    done
}

#Extract amino acid subsets from ARG loci. These are used for prediction of IS elements and integrons
rm -f temp.??? pids
for c in c???
do
    num=$(echo $c | sed 's/c//')
    extract_fun $c $num &
    echo $! >> pids
done

#Wait for all processes to finish
all_pids=$(cat pids)
for pid in $all_pids; do
    echo $pid
    wait $pid
done

cat temp.??? | sort -u | awk '{print ">"$1"\n"$2}' > ARG_loci.faa
grep -c '>' ARG_loci.faa
#4233790 fasta headers

#sed 's/_#_.*//' ARG_loci.faa > ARG_loci_edit.faa
rm temp.??? c???

#Get index of ARG proteins from faa file. ARGs passing filters are stored in CARD_refseq_filt.diamond
#The temp.coord and ARG_list is in the same order. They can be pasted together to get list of Locus:ARG combinations
awk '{print $1,$2}' CARD_refseq_filt.diamond | awk -F'|' '{print $1,$3}' | awk '{print $1,$3}' > ARG_list
paste temp.coord ARG_list | sed 's/\t/:/' > locus_arg_list



#Get chr/plasmid info
split --number=l/$numthreads -d locus_arg_list l

replicon_fun () {
    inlist=$1
    num=$2
    acc_list=$(cat $inlist | sed 's/\.[0-9]*.*//')
    rm -f templist$num
    for acc in $acc_list
    do
        cacc=$(grep "$acc" /database/ncbi/ncbi_complete_genomes_121219/refseq_complete_chr_headers.txt)
        if [ ! -z "$cacc" ]
        then
            replicon=Chr
            echo $acc $replicon >> templist.$num
        else
            replicon=Plasmid
            echo $acc $replicon >> templist.$num
        fi
    done
}
rm -f templist.?? pids
for l in l??
do
    number=$(echo $l | sed 's/l//')
    replicon_fun $l $number &
    echo $! >> pids
done


#Wait for all processes to finish
all_pids=$(cat pids)
for pid in $all_pids; do
    echo $pid
    wait $pid
done



cat templist.* > replicon_acc
paste locus_arg_list replicon_acc > locus_arg_list_repl
#176888 locus_arg_list_repl

rm -f temp templist.?? l??



awk '{print $1}' locus_arg_list > temp
seqtk subseq ARG_loci.faa temp | fasta_formatter -t | awk '{print ">"$1"\n"$2}' > loci_ARGs.faa
#loci_ARGs.faa should be same length as locus_arg_list (176888) 



#### Find the distance between ARGs and IS elements if both were distributed randomly ####
#Important note: This part of the analyses is not included in the current version of the manuscript. 
grep '>' ARG_regions.fasta | sed 's/:.*//' | sort -u > ARG_accs
#13298

#Do diamond blastp for all proteins against ISfinder db
#For IS prediction with PROKKA, cutoffs are qcov 90 and E-value 1E-30. Use the same cutoffs
$diamond blastp --db /opt/prokka/db/kingdom/Bacteria/IS --query /DATA_1/tkn/temp_prodigal/all_refseq_complete_prodigal_single_meta.faa \
--outfmt 6 --threads 70 --evalue 10E-30 --query-cover 90 --max-target-seqs 1 > prodigal_all_ISfinder.blastp 
#843626 prodigal_all_ISfinder.blastp
#790628

grep '>' /database/ncbi/ncbi_complete_genomes_121219/refseq_complete_plasmid_headers.txt | awk '{print $1}' | sed 's/>//' > all_refseq_plasmid_acc
grep '>' /database/ncbi/ncbi_complete_genomes_121219/refseq_complete_chr_headers.txt | awk '{print $1}' | sed 's/>//' > all_refseq_chr_acc
#14280 all_refseq_plasmid_acc
#16785 all_refseq_chr_acc

rm -f simulated_dist.txt
acc_fun () {
    file_in=$1
    acc_list=$(cat $file_in | sed 's/>//')
    rm -f $file_in.temp
    for i in $acc_list
    do
        #echo $i
        acc_p=$(grep "$i" all_refseq_plasmid_acc)
        if [ -z "$acc_p" ]
        then
            #echo "$i is a chromosome"
            acc=$(grep "$i" all_refseq_chr_acc)
            seqlength=$(grep -w "$acc" /database/ncbi/ncbi_complete_genomes_121219/refseq_complete_seqlengths | awk '{print $2}')
            acc2=$(echo $acc | sed 's/\.[0-9]*//')
            #num_is=$(grep "$i" protein_IS_all.blastp | wc -l)
            #num_arg=$(grep "$i" protein_ARO_all_filt.blastp | wc -l)
            num_arg=$(grep "$i" diamond_hits_rgi_cutoff.pass | wc -l)
            num_is=$(grep "$i" prodigal_all_ISfinder.blastp | wc -l)
            #echo $seqlength $num_is
            if [ $num_is = 0 ]
            then
                #echo $i $acc $seqlength $num_is $num_arg $seqlength "Chr" >> simulated_dist.txt
                echo $i $acc $seqlength $num_is $num_arg $seqlength "Chr" >> $file_in.temp
            else
                mean_dist=$(echo $seqlength $num_is | awk '{print $1/(2*($2+1))}')
                #echo $i $acc $seqlength $num_is $num_arg $mean_dist "Chr" >> simulated_dist.txt
                echo $i $acc $seqlength $num_is $num_arg $mean_dist "Chr" >> $file_in.temp
            fi
        else
            #echo "$i is a plasmid"
            acc=$acc_p
            seqlength=$(grep -w "$acc" /database/ncbi/ncbi_complete_genomes_121219/refseq_complete_seqlengths | awk '{print $2}')
            #num_is=$(grep "$i" protein_IS_all.blastp | wc -l)
            #num_arg=$(grep "$i" protein_ARO_all_filt.blastp | wc -l)
            num_arg=$(grep "$i" diamond_hits_rgi_cutoff.pass | wc -l)
            num_is=$(grep "$i" prodigal_all_ISfinder.blastp | wc -l)
            #echo $seqlength $num_is
            if [ $num_is = 0 ]
            then
                #echo $i $acc $seqlength $num_is $num_arg $seqlength "Plasmid" >> simulated_dist.txt
                echo $i $acc $seqlength $num_is $num_arg $seqlength "Plasmid" >>  $file_in.temp
            else
                mean_dist=$(echo $seqlength $num_is | awk '{print $1/(2*($2+1))}')
                #echo $i $acc $seqlength $num_is $num_arg $mean_dist "Plasmid" >> simulated_dist.txt
                echo $i $acc $seqlength $num_is $num_arg $mean_dist "Plasmid" >> $file_in.temp
            fi
        fi
    done
}

split --number=l/$numthreads -d ARG_accs s
rm -f pids
for i in s??
do
    acc_fun $i &
    echo $! >> pids
done


#Wait for all processes to finish
all_pids=$(cat pids)
for pid in $all_pids; do
    echo $pid
    wait $pid
done

cat s??.temp > simulated_dist.txt
rm s??.temp s??
#13298 simulated_dist.txt



#Run blast against ISfinder db with loci_ARGs.faa as query
#ISfinder db modified with below to more easily get also IS family names
$diamond blastp --db /opt/prokka/db/kingdom/Bacteria/IS_mod --query ARG_loci.faa \
--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp \
--threads 70 --evalue 10E-30 --max-target-seqs 1 --query-cover 90 > protein_IS_all.blastp 
#123545 protein_IS_all.blastp


#Merge ARG and IS coordinates for extracting closest IS elements
#Columns:
#Locus Protein Start End ARO/IS
awk '{print $1,$2}' locus_arg_list | sed 's/:/ /' | awk -F'_' '{print $0,$8,$9}' | awk '{print $1,$2,$3,$4,$5}' > locus_arg_list2
#NZ_CP026125.1_1037050_1064500_3000775 NZ_CP026125.1_pr1052_1049220_1052330 ARO:3000775 1049220 1052330
awk '{print $1,$2}' protein_IS_all.blastp | sed 's/:/ /' | awk -F'_' '{print $0,$8,$9}' | awk '{print $1,$2,$3,$4,$5}' > protein_IS_all.blastp2
#NC_000913.3_382960_408457_3004612 NC_000913.3_pr360_391739_392605 IS3:IS3 391739 392605
cat locus_arg_list2 protein_IS_all.blastp2 | sort -nk1,1 > ARG_IS_list
#300433 ARG_IS_list

locus_fun () {
    inlist=$1
    num=$2
    loci_list=$(cat $inlist | awk '{print $1}')
    for i in $loci_list
    do
        acc=$(echo $i | awk -F'_' '{print $1"_"$2}')
        replicon=$(grep "$acc" simulated_dist.txt | awk '{print $7}')
        sim_cutoff=$(grep "$acc" simulated_dist.txt | awk '{print $6}')
        locus=$(echo $i | awk -F':' '{print $1}')
        protein=$(echo $i | awk -F':' '{print $2}')
        start=$(echo $protein | awk -F'_' '{print $3}')
        end=$(echo $protein | awk -F'_' '{print $4}')
        aro=$(echo $i | awk -F':' '{print $1}' | awk -F'_' '{print $6}')
        is_up=$(grep $locus ARG_IS_list | grep -B100 $protein | grep -v "ARO:" | tail -n1)
        is_down=$(grep $locus ARG_IS_list | grep -A100 $protein | grep -v "ARO:" | head -n1)
        #If there are no IS elements, print a 0 to indicate that no IS elements are found in the given direction
        #Zeroes are used R for calculating ratios per ARO of IS proximity. A ratio of 1 indicates that 100% of genes in an ARO have IS elements in proximity 
        if [ -z "$is_up" ]
        then
            dist_up=0
        else
            dist_up=$(echo $is_up | awk -v arg_start=$start '{print arg_start-$5}')
        fi
        
        if [ -z "$is_down" ]
        then
            dist_down=0
        else
            dist_down=$(echo $is_down | awk -v arg_end=$end '{print $4-arg_end}')
        fi
        
        #For resistance passenger genes (annotated in ISfinder), IS elements have same start/stop as ARGs, resulting in 
        #negative distance(s) for ARGs to IS. Fix these by changing to 1 (0 indicates that there are no IS elements)
        if [ "$dist_down" -lt 0 ]
        then
            dist_down=1
        else
            dist_down=$dist_down
        fi
        if [ "$dist_up" -lt 0 ]
        then
            dist_up=1
        else
            dist_up=$dist_up
        fi
        
        if [ -z "$is_up" ]
        then
            up_name="None"
        else
            up_name=$(echo $is_up | awk '{print $3}')
        fi
        
        if [ -z "$is_down" ]
        then
            down_name="None"
        else
            down_name=$(echo $is_down | awk '{print $3}')
        fi
        #NC_000913_556145_580817_3004039_17_g3	ARO:3004039 311 6090 IS3 IS5 Chr 28302.8 3200.5
        echo $locus "ARO:"$aro $dist_up $dist_down $up_name $down_name $replicon $sim_cutoff | awk '{print $0,($3+$4)/2}' >> ARG_IS_distances2.$num
    done
}

rm -f x?? ARG_IS_distances2.[0-9][0-9] pids
split --number=l/$numthreads -d locus_arg_list
for x in x??
do
    number=$(echo $x | sed 's/x//')
    locus_fun $x $number &
    echo $! >> pids
done

#Wait for all processes to finish
all_pids=$(cat pids)
for pid in $all_pids; do
    echo $pid
    wait $pid
done


cat ARG_IS_distances2.[0-9][0-9] | sort -u > ARG_IS_distances2.txt
#176888 ARG_IS_distances2.txt
rm x?? temp.* ARG_IS_distances2.[0-9][0-9]


#Calculate means of up- and down-distances to IS elements. If one of the distances are 0 (there are no IS elements), print the other one
#If both distances are 0, (there are no IS elements in either direction), print 0
awk '{if ($3 > 0 && $4 > 0) print $0,($3+$4)/2}' ARG_IS_distances2.txt > temp1
awk '{if ($3 == 0 && $4 == 0) print $0,"0"}' ARG_IS_distances2.txt > temp2
awk '{if ($3 == 0 && $4 > 0) print $0,$4}' ARG_IS_distances2.txt > temp3
awk '{if ($3 > 0 && $4 == 0) print $0,$3}' ARG_IS_distances2.txt > temp4 

cat temp[1-4] | sed 's/ /\t/' | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$10}' | sed 's/\.[0-9]*_/_/' > ARG_IS_distances.dat
sed -i 's/None/No_IS/g' ARG_IS_distances.dat
#176888 ARG_IS_distances.dat
rm temp? ARG_IS_distances2.txt


#Filtering on simulated ARG:IS distances is disabled in version 0.3

#Filter and extract unclustered regions based on mean distance against expected random distances
#awk '{if ($8 > $9) print $0}' ARG_IS_distances.dat | sed 's/\.[0-9]*_/_/' | sort > ARG_IS_distances_pass_exp.dat
#sed 's/\.[0-9]*_/_/' ARG_IS_distances.dat | sort > ARG_IS_distances_pass_exp.dat
#awk '{print $1}' ARG_IS_distances_pass_exp.dat | sed -e 's/_[0-9]*_g[0-9]*$//' -e 's/\.[0-9]*_/_/' > passing_regions
#Get more precise region names
#awk '{print $1}' ARG_IS_distances_pass_exp.dat | sed -e 's/\([A0-Z9]*_[A0-Z9]*_[0-9]*_[0-9]*_\)/\1ARO/' -e 's/\(ARO[0-9]*\)_/\1:/' > passing_regions_headers
#sed -e 's/\.[0-9]*:/_/' -e 's/-/_/' -e 's/|/_/' ARG_regions.fasta > ARG_regions_edit.fasta
#Extract regions from simplified headers
#seqtk subseq ARG_regions_edit.fasta passing_regions > passing_regions.fas



#Replace the simple headers with the more precise headers
#fasta_formatter -i passing_regions.fasta -t -o passing_regions.tab
#paste passing_regions_headers passing_regions.tab | awk '{print ">"$1"\n"$3}' > passing_regions.fas
#paste passing_regions_headers passing_regions.tab | awk '{print ">"$1"\n"$3}' | sed -e 's/ARO//' -e 's/:[0-9]*_g[0-9]*//'> passing_regions.fas




#Get mechanism for passing regions for counting number of loci for each major mechanism
aro_mech_fun () {
    inlist=$1
    num=$2
    rm -f templist.$num
    aro_list=$(cat $inlist)
    for aro in $aro_list
    do
        dna_acc=$(grep -w "$aro" /DATA_1/tkn/ProxARG19/CARD_db/aro_index.tsv | sort -u | awk -F'\t' '{print $7}' | sort -u)
        #Get general resistance mechanism
        mech=$(grep -w "$dna_acc" /DATA_1/tkn/ProxARG19/CARD_db/aro_categories_index.tsv | awk -F'\t' '{print $5}' | sort -u | sed 's/ /_/g')
        echo $mech >> templist.$num
    done
}

rm -f x?? templist.?? pids
#awk -F'_' '{print $6}' passing_regions > temparo.list
awk '{print $2}' ARG_IS_distances.dat > temparo.list
split --number=l/$numthreads -d temparo.list
for x in x??
do
    number=$(echo $x | sed 's/x//')
    aro_mech_fun $x $number &
    echo $! >> pids
done

#Wait for all processes to finish
all_pids=$(cat pids)
for pid in $all_pids; do
    echo $pid
    wait $pid
done

cat templist.?? | sort | uniq -c | sed 's/^ *//' > passing_regions_mech
rm -f templist.?? temparo.list x??
#169216 fasta entries passing mean ARG:IS distance > expected by random chance
rm -f passing_regions.tab passing_regions.fasta passing_regions_headers ARG_regions_edit.fasta

#The below is not used in current version

#Get regions that fail the mean distance < expected by random criterium.
#Per ARO, count how many fail and how many pass. Divide them to get a ratio. A high ratio indicates that a given ARO is often placed in a "highly mobilized" context
#awk '{if ($8 < $9) print $0}' ARG_IS_distances.dat | awk '{print $2}' | sort | uniq -c | awk '{print $2,$1}' | sort -n > exp_fail
#awk '{if ($8 > $9) print $0}' ARG_IS_distances.dat | awk '{print $2}' | sort | uniq -c | awk '{print $2,$1}' | sort -n > exp_pass

#join -a 1 -a 2 -j 1 -e '0' -o '2.1 2.2 1.2' exp_fail exp_pass | awk '{print $0,$3/($2+$3)}' | sed '/^0/d' > highly_mobilized_aros
#1154 highly_mobilized_aros

#Cluster similar regions with usearch
sed -e 's/:/_/' -e 's/-/_/' -e 's/|/_/' ARG_regions.fasta > ARG_regions_edit.fasta

mkdir AROs_otus
#Start with the most abundant AROs, since these take longer to cluster
aro_list=$(grep '>' ARG_regions_edit.fasta | awk -F'_' '{print $6}' | sort | uniq -c | sort -nrk1,1 | awk '{print $2}' )

usearch_fun () {
    aro_num2=$1
    grep -A1 "_$aro_num2$" ARG_regions_edit.fasta | grep -v '\-\-' > AROs_otus/$aro_num2.fas
    nohup $usearch -cluster_fast AROs_otus/$aro_num2.fas -centroids AROs_otus/$aro_num2.otus.fa -id 0.99 \
    -uc AROs_otus/$aro_num2.otus.uc -sizeout -minqt 1 -query_cov 0.9 -target_cov 0.9 -strand both -threads 30 -sort length -maxhits 1 &
}

for i in $aro_list
do
    ucount=$(pgrep 'usearch' | wc -l)
    #echo $ucount
    if [ $ucount -lt 80 ]
    then
        #echo $ucount
        usearch_fun $i &
    else
        #echo $ucount
        #echo "Waiting"
        sleep 60
        usearch_fun $i &
    fi
done

#Wait until the last usearch has finished. Takes a long time.


#Extract clustered regions (that pass the simulated cutoffs) to RTUs
grep -h '>' AROs_otus/[0-9][0-9][0-9][0-9][0-9][0-9][0-9].otus.fa | sed -e 's/>//' -e 's/;size=/ size_/' -e 's/;$//' -e 's/\.[0-9]*_/_/' > passing_clusters
#53895 passing_clusters
awk '{print $1}' passing_clusters > cluster_list
split --number=l/$numthreads -d cluster_list

cluster_fun() {
    cluster_list=$1
    num=$2
    clust_list=$(cat $cluster_list)
    rm -f ARG_IS_distances_clusters.$num
    for i in $clust_list
    do
        size=$(grep -w "$i" passing_clusters | sort -u | awk '{print $2}')
        region=$(echo "$i" | sed -e 's/ARO//' -e 's/:/_/')
        #	grep "$i" ARG_IS_distances_pass_exp.dat | awk -v size=$size '{print $0,size}' >> ARG_IS_distances_clusters.$num
        if grep -w -q "$region" ARG_IS_distances.dat
        then
            grep -w "$region" ARG_IS_distances.dat | awk -v size=$size '{print $0,size}' >> ARG_IS_distances_clusters.$num
        else
            echo "$region";
        fi
    done
}

rm -f pids
for i in x??
do
    num=$(echo $i | sed 's/x//')
    cluster_fun $i $num &
    echo $! >> pids
done

all_pids=$(cat pids)
for pid in $all_pids; do
    echo $pid
    wait $pid
done

cat ARG_IS_distances_clusters.?? | sort -u > ARG_IS_distances_clusters.dat
rm ARG_IS_distances_clusters.?? x??
#53895 ARG_IS_distances_clusters.dat

#Extract the resistance mechanism for grouping in R
awk '{print $1}' ARG_IS_distances_clusters.dat > locus_list
rm -f ARG_IS_distances_mech.txt
rm -f z??
split --number=l/$numthreads -d locus_list z

mech_fun () {
    list=$1
    locus_list=$(cat $list)
    for i in $locus_list
    do
        curr_acc=$(grep -w "$i" ARG_IS_distances_clusters.dat | awk '{print $2}' | sort -u)
        dna_acc=$(grep -w "$curr_acc" /DATA_1/tkn/ProxARG19/CARD_db/aro_index.tsv | sort -u | awk -F'\t' '{print $7}' | sort -u)
        #Get general resistance mechanism
        mech=$(grep -w "$dna_acc" /DATA_1/tkn/ProxARG19/CARD_db/aro_categories_index.tsv | awk -F'\t' '{print $5}' | sort -u | sed 's/ /_/g')
        #Get also the specific mechanism (e.g. beta-lactamase)
        mech2=$(grep -w "$dna_acc" /DATA_1/tkn/ProxARG19/CARD_db/aro_categories_index.tsv | awk -F'\t' '{print $3}' | sort -u | sed 's/ /_/g')
        dist_data=$(grep -w "$i" ARG_IS_distances_clusters.dat)
        echo "$dist_data $mech $mech2" >> ARG_IS_distances_mech.txt
    done
}

rm -f pids
for b in z??
do
    mech_fun $b &
    echo $! >> pids
done

all_pids=$(cat pids)
for pid in $all_pids; do
    echo $pid
    wait $pid
done


sed -i 's/ /\t/g' ARG_IS_distances_mech.txt
sed -i 's/None/No_IS/g' ARG_IS_distances_mech.txt
rm z??
#53895 ARG_IS_distances_mech.txt

##### Integron_finder #####

cat AROs_otus/*fa > all_AROs_otus.fasta
sed -e 's/\.[0-9]*//' -e 's/:/_/' -e 's/-/_/' -e 's/ARO:/ARO_/' -e 's/;size=/_size/' -e 's/|/_/' -e 's/;//' all_AROs_otus.fasta > all_AROs_otus.fas

#Run integron_finder on ARG regions
#Unlimited line width on ARG regions
fasta_formatter -i all_AROs_otus.fas -o all_AROs_otus.tab -t
split --number=l/$numthreads -d all_AROs_otus.tab
mkdir integron_finder

int_fun () {
    tab_in=$1
    num=$2
    awk '{print ">"$1"\n"$2}' $tab_in > temp$num.fas
    region_list=$(grep '>' temp$num.fas | sed 's/>//')
    for a in $region_list
    do
        suffix=$(head /dev/urandom | tr -dc A-Za-z0-9 | head -c 13 ; echo '')
        grep -w -A1 "$a" temp$num.fas > temp_int.fa.$suffix
        #integron_finder --cpu 10 --outdir integron_finder/"$a" --linear --prodigal /usr/bin/prodigal temp_int.fa.$suffix
        rm -rf integron_finder/"$a"
        icount=$(pgrep 'integron_finder' | wc -l)
        #echo $icount
        #echo $a
    if [ $icount -lt 80 ]
    then
        #echo $icount
        integron_finder --local_max --cpu 10 --outdir integron_finder/"$a" --linear --prodigal /usr/bin/prodigal temp_int.fa.$suffix
    else
        #echo $icount
        #echo "Waiting"
        sleep 60
        integron_finder --local_max --cpu 10 --outdir integron_finder/"$a" --linear --prodigal /usr/bin/prodigal temp_int.fa.$suffix
    fi
    rm temp_int.fa.$suffix
    done
}


#Integron_finder prints a lot to stdout. Run in screen.

for i in x??
do
    num=$(echo $i | sed 's/x//')
    int_fun $i $num &
done

### Takes quite a while to process all loci
#ll integron_finder/ | wc -l
#49121


#Integron hits
#Since integron_finder performs prodigal gene prediction on input DNA loci, overlapping loci can have multiple ARGs (that belong to other loci). \
Get the ARG coordinates from loci_ARGs.faa that have the format:
#>NZ_CP036368.1_1_22388_3002547:NZ_CP036368.1_pr10_9619_10218
#From integron_finder, the format is: 
#NZ_CP036337_192836_218375_3000165_size1_20 start end

#Adjust by start and stop of gene coordinates by -1 to properly adjust gene coordinates with locus coordinates
find integron_finder/ -name 'temp_int.fa.integrons' -exec grep "protein" {} \; | grep -v 'intI' | awk '{print $3,$4-1,$5-1}' | sort -u > int_passengers
awk -F'_' '{print $0,$3}' int_passengers | awk '{print $1,$4+$2"_"$4+$3}' | sed 's/_size[0-9]*_[0-9]*//' > int_passengers2
#43506 int_passengers

#Reformat fasta headers from loci_ARGs.faa. The objective is to extract info on ARGs in integrons from ARG_IS_distances_mech.txt with the format:
#NC_000907_1665901_1691635_3003953
grep '>' loci_ARGs.faa | sed 's/\.[0-9]*//' > loci_ARGs_mod_headers

integron_passenger_fun () {
    plist=$1
    number=$2
    passenger_list=$(cat $plist | awk '{print $1}' )
    for p in $passenger_list
    do
        #echo $p
        coordinates=$(grep $p int_passengers2 | awk '{print $2}')
        for c in $coordinates
        do
            #echo $c
            #grep $p loci_ARGs_mod_headers
            #grep $p loci_ARGs_mod_headers | grep $c
            grep $p loci_ARGs_mod_headers | grep $c | sort -u | sed 's/>//' >> ARG_integrons_passengers0."$num"
        done
        arg_ints=$(awk -F':' '{print $1}' ARG_integrons_passengers0."$num" | sort -u)
        for a in $arg_ints
        do
            grep -w "$a" ARG_IS_distances_mech.txt >> ARG_integrons_passengers."$num"
        done
    done
}

rm -f x?? ARG_integrons_passengers* pids
split --number=l/$numthreads -d int_passengers2

for x in x??
do
    num=$(echo $x | sed 's/x//')
    integron_passenger_fun $x $num &
    echo $! >> pids
done

all_pids=$(cat pids)
for pid in $all_pids; do
    echo $pid
    wait $pid
done

cat ARG_integrons_passengers.[0-9]* | sort -u > ARG_integrons_ISProx.txt
#3723 ARG_integrons_ISProx.txt

rm x?? ARG_integrons_passengers*

#Get parent acc and then hit accession numbers. Use parent region to extract aro info

awk '{print $1}' passing_clusters > cluster_list
rm -f c??
split --number=l/$numthreads -d cluster_list c

sed -i 's/\.[0-9]*_/_/g' AROs_otus/*otus.uc

#In the function below, there are included some sanity checks for debugging. There shouldn't be a problem now, but the files missing_* can be checked. These should be empty.
tax_fun () {
    parent_file=$1
    num=$2
    parents=$(cat $parent_file)
    rm -f aro_tax.$num
    for parent in $parents
    do
        #echo $parent
        parent2=$(echo $parent | sed -E 's/(_)([0-9]*$)/\|\2/')
        aro=$(echo $parent | grep -o "_[0-9]*$" | sed 's/_//')
        #aro=$(echo $parent | grep -o "ARO[0-9]*")
        parent3=$(echo $parent | sed -e 's/\.[0-9]*//' -e 's/:/_/g' -e 's/-/_/' -e 's/ARO//')
        aro_info=$(grep -m 1 "$parent3" ARG_IS_distances_mech.txt)
        #echo $aro_info
        if [ -z "$aro_info" ]
        then
            echo "NOTHING TO SEE HERE (No aro_info) $parent $aro" >> missing_parents
        else
            accessions=$(grep $parent AROs_otus/$aro.otus.uc | awk '{print $9}' | awk -F'_' '{print $1"_"$2}' | sort -u)
            num_acc=$(grep $parent AROs_otus/$aro.otus.uc | awk '{print $9}' | sed 's/:.*//' | sort -u | wc -l)
            #echo "There are $num_acc accessions under $parent2"
            if [ -z "$accessions" ]
            then
                echo "NOTHING TO SEE HERE (No accessions) $parent $aro" >> missing_accs
            else
                for acc in $accessions
                do
                    if [ -z "$acc" ]
                    then
                        echo "NOTHING TO SEE HERE (No accession) $parent $aro" >> missing_acc
                    else
                        #echo $acc
                        taxid=$(grep -m 1 -w "$acc" /DATA_1/tkn/ProxARG19/nucl_gb.accession2taxid.subset | awk '{print $3}' )
                        if [ -z "$taxid" ]
                        then
                            echo "NOTHING TO SEE HERE (No tax_id) $parent $aro $acc" >> missing_taxid
                        else
                            taxonomy=$(grep -m 1 -w "^$taxid" /DATA_1/tkn/ProxARG19/fullnamelineage_edit.dmp | sed -e 's/ /_/g' | awk -F'|' '{print $3}' | awk -F';_' '{print $1";_"$2";_"$3";_"$4";_"$5";_"$6";_"$7}')
                            #Include species
                            #taxonomy=$(grep -m 1 -w "^$taxid" /DATA_1/tkn/ProxARG19/fullnamelineage_edit.dmp | sed -e 's/ /_/g' | awk -F'|' '{print $3}' | awk -F';_' '{print $1";_"$2";_"$3";_"$4";_"$5";_"$6";_"$7";_"$8}')
                            #echo $taxonomy
                            echo $aro_info $acc $taxonomy >> aro_tax."$num"
                        fi
                    fi
                done
            fi
        fi
    done
}

rm -f pids aro_tax.??
for i in c??
do
    num=$(echo $i | sed 's/c//')
    tax_fun $i $num &
    echo $! >> pids
done

all_pids=$(cat pids)
for pid in $all_pids; do
    echo $pid
    wait $pid
done



sed -e 's/ /\t/g' -e 's/;_/\t/g' aro_tax.[0-9]* | sort -u > aro_tax.Rdat
#176612 aro_tax.Rdat
rm aro_tax.?? c??

#Get a table of CARD protein accession numbers per ARO
aros=$(awk '{print $2}' ARG_IS_distances_mech.txt | sort -u)
rm -f aro_acc.txt
for aro in $aros
do
    acc=$(grep -m1 "$aro" /DATA_1/tkn/ProxARG19/CARD_db/protein_fasta_protein_homolog_model.fasta | awk -F'|' '{print $2}')
    nacc=$(grep -m1 "$aro" /DATA_1/tkn/ProxARG19//CARD_db/nucleotide_fasta_protein_homolog_model.fasta | awk -F'|' '{print $2}')
    echo $aro $nacc $acc >> aro_acc.txt
done

#The following files are imported in R for analyses

#File							Product of
#ARG_IS_distances_mech.txt		this script
#simulated_dist.txt				this script
#aro_tax.Rdat					this script
#aro_acc.txt					this script
#diamond_hits_rgi_cutoff.pass	this script
#diamond_hits_rgi_cutoff.fail	this script
#ARG_integrons_ISProx.txt		this script
#highly_mobilized_aros			this script # Not used in current version
#passing_regions_mech           this script
#16S_IS_distance				16S_control.sh
#plasmid_tax					all_acc_to_tax.sh
#chr_tax						all_acc_to_tax.sh
#refseq_1st_chr.tax				all_acc_to_tax.sh
#post_clust_tax_c				all_acc_to_tax.sh
#protein_homolog.tax			fasta_formatter -i /DATA_1/tkn/ProxARG19/CARD_db/protein_fasta_protein_homolog_model.fasta -t | grep -o "\[.*\]" | awk '{print $1}' | sed -e 's/\[//' -e 's/\]//' | sort | uniq -c | sed -e 's/^ *//' > protein_homolog.tax

#Prepare tar.gz file for transferring to local computer and analyses in RStudio
tar -czvf argmob_data.tar.gz ARG_IS_distances_mech.txt simulated_dist.txt \
aro_tax.Rdat aro_acc.txt diamond_hits_rgi_cutoff.pass diamond_hits_rgi_cutoff.fail \
ARG_integrons_ISProx.txt passing_regions_mech \
16S_IS_distance plasmid_tax chr_tax refseq_1st_chr.tax post_clust_tax_c \
protein_homolog.tax

#Final cleanup
rm -f pids int_passengers* cluster_list all_AROs_otus.tab \
all_AROs_otus.fasta locus_list passing_clusters ARG_IS_list \
ARG_IS_distances2.txt protein_IS_all.blastp2 locus_arg_list2 \
nohup.out passing_regions ARG_IS_distances_pass_exp.dat \
temp* protein_IS_all.blastp 
