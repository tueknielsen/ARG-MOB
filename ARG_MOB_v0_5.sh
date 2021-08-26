#!/usr/bin/env bash

#ARG-MOB version 0.5
#Tue Kjærgaard Nielsen
#August 12 2021

#Disclaimer:
#This bash script performs identification of ARGs and MGEs in RefSeq complete bacterial genomes and preparation for further analyses in R. 
#It is mainly provided as documentation but has been written such that analyses can be replicated. All used databases have been uploaded to GigaDB and are required for replication of analysis. 
#Before running this script, it is required that you change the paths below for applied programs and databases. 

#Many tasks are parallelized via functions that are set to run 100 tasks in parallel. This may not be optimal for all computers/servers and can be changed on line 23 
#Analyses were performed on our server "thoth" with the following specifications: 
#HPE ProLiant DL580 Gen10 Server equipped with four Intel Xeon Gold 6254 processors and 3 TB RAM running Ubuntu 18.04 LTS

#During development, the pipeline was named 'ProxARG' and this name may still occur in the script.

#At many places in the script, there are checkpoints inserted that show the expected number of lines in important (intermediate) files.

#Important note: some ARGs and IS elements are encoded by two or more ORFs. These will be overrepresented in final results, but should not overall affect results much.

#Define the number of threads used for parallel processing
numthreads=100

#Paths to programs. Edit these to match your paths
diamond=/home/tkn/Programs/diamond  #v0.9.25.126
#ncbi-genome-download=/usr/local/bin/ncbi-genome-download    #0.2.11 #Not used in main analysis
fasta_formatter=/usr/bin/fasta_formatter    #0.0.14. Part of the FASTX toolkit.
prodigal=/usr/bin/prodigal  #V2.6.3
samtools=/home/tkn/.local/opt/anaconda3/bin/samtools    #1.10
seqtk=/home/tkn/.local/opt/anaconda3/bin/seqtk  #1.3-r106
usearch=usearch_64  #v11.0.667_i86linux64   We used the paid 64 bit version. The free 32 bit version has not been tested. The free alternative VSEARCH could likely replace USEARCH here. 
integron_finder=/usr/local/bin/integron_finder  #1.5.1
barrnap=/usr/local/bin/barrnap  #0.9

#Databases required for running analysis. Start by unpacking ARGMOB_dbs.tar.gz and then define the variable $database_dir
#tar -xzvf ARGMOB_dbs.tar.gz
database_dir=/DATA_1/tkn/ARGMOB/ARGMOB_dbs


#Below are some notes on how data was preprocessed for the analyses. \
#The actual analysis starts on line 200


#Refseq complete genomes downloaded on Dec. 12 2019 using the command: \
#ncbi-genome-download -l complete -p 20 -F fasta bacteria
#refseq_complete.fasta

#Chromosomal headers from Refseq
#refseq_complete_chr_headers.txt

#Plasmid headers from Refseq
#refseq_complete_plasmid_headers.txt

#Gene predictions on complete refseq genomes.
#all_refseq_complete_prodigal_single_meta.faa

#CARD database downloaded on Dec. 13 2019. Diamond database prepared with: \
#diamond makedb --in /DATA_1/tkn/ProxARG19/CARD_db/protein_fasta_protein_homolog_model.fasta -d /DATA_1/tkn/ProxARG19/CARD_db/protein_fasta_protein_homolog_model
#protein_fasta_protein_homolog_model

#CARD ARO index
#aro_index.tsv

#RGI database
#proteindb.fsa

#NCBI fullnamelineage
#fullnamelineage_edit.dmp

#ISfinder database as included in prokka 
#IS

#NCBI acc to taxid
#nucl_gb.accession2taxid.subset

#Note on 'nucl_gb.accession2taxid.subset': \
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

#Note on 'ncbi_fullname' database: \
#Superphyla are confusing the taxonomy order, since they push the phylum one level down for bacteria with superphyla assigned. \
#Remove supergroups from fullnamelineage.dmp
#sed -e 's/Terrabacteria group; //' -e 's/PVC group; //' \
#-e 's/FCB group; Bacteroidetes\/Chlorobi group; //' \
#-e 's/FCB group; //' -e 's/Cyanobacteria\/Melainabacteria group; //' fullnamelineage.dmp > fullnamelineage_edit.dmp


#Download all gff annoations for extracting 16S control genes (Not part of main analysis)
#ncbi-genome-download -l complete -p 20 -F gff bacteria -o complete_gff -v
#ncbi-genome-download -l complete -p 20 -F faa bacteria -o complete_faa -v


#Genes were predicted from all RefSeq complete genomes with Prodigal. \
#Since Prodigal trains itself on the nucleotide sequence before predicting genes, \
#Prodigal was run separately on each genus in the database. Furthermore, \
#some genes are missed in "meta" mode compared to "single" mode and vice versa. \
#Hence, Prodigal was run twice per genus, once in single and once in meta modes. \
#The two runs were concatenated per genus and redundancies (genes predicted in both modes) \
#were removed. The commands below were performed only once to produce the file all_refseq_complete_prodigal_single_meta.faa


##Prodigal prediction doesn't really work well with a collection of very varying genomes. Do subprocesses according to genera.
#mkdir prodigal_genus
#refseq_genera=$(awk '{print $2}' refseq_complete_headers.txt | sed -e 's/\[//' -e 's/\]//' -e "s/[']//g" | sort | uniq -c | sort -nrk1,1 | awk '{print $2}')
#seqtk_prodigal_fun () {
#   genus=$1
#   echo "Extracting $genus"
#   seqtk subseq /database/ncbi/ncbi_complete_genomes_121219/refseq_complete_w0.fasta prodigal_genus/"$genus".list > prodigal_genus/"$genus".fasta
#   echo "Single Prodigal on $genus"
#   prodigal -p single -a prodigal_genus/"$genus"_single.faa -d prodigal_genus/"$genus"_single.fna -i prodigal_genus/"$genus".fasta
#   echo "Meta Prodigal on $genus"
#   prodigal -p meta -a prodigal_genus/"$genus"_meta.faa -d prodigal_genus/"$genus"_meta.fna -i prodigal_genus/"$genus".fasta
#}
#for genus in $refseq_genera
#    do
#    awk -v gen=$genus '{if ($2 == gen) print $0}' refseq_complete_headers.txt | sed 's/>//' > prodigal_genus/"$genus".list
#    pcount=$(pgrep 'seqtk' | wc -l)
#    if [ $pcount -lt 50 ]
#    then
#        seqtk_prodigal_fun $genus &
#    else
#        sleep 120
#        seqtk_prodigal_fun $genus &
#    fi
#done

##Concatenate single and meta prodigal predictions and remove redundancies (per genus)
#cat prodigal_genus/*single.faa > prodigal_per_genus.faa
#59015148
#cat prodigal_genus/*meta.faa > prodigal_per_genus_meta.faa
#58964013 fasta headers
#Format the fasta headers
#awk -F'#' '{print $1,$2,$3}' prodigal_per_genus.faa | sed -e 's/ \+\([0-9]*\) \+\([0-9]*\)/_\1_\2/' -e 's/__//' | fasta_formatter -t | awk -F'_' '{print $1"_"$2"_"$4"_"$5"\t"$3"\ts"}' > temp1
#awk -F'#' '{print $1,$2,$3}' prodigal_per_genus_meta.faa | sed -e 's/ \+\([0-9]*\) \+\([0-9]*\)/_\1_\2/' -e 's/__//' | fasta_formatter -t | awk -F'_' '{print $1"_"$2"_"$4"_"$5"\t"$3"\tm"}' > temp2

#Merge gene calls from single and meta modes. Remove redundancies. 
#cat temp1 temp2 | sort -u -t$'\t' -k 1,1 -k 2,2 | awk -F$'\t' '{print $1"_pr"$4$3"\t"$2}' | sed 's/ //g' > all_refseq_complete_prodigal_single_meta.faa.tab
#62396260 all_refseq_complete_prodigal_single_meta.faa.tab

#Tab to fasta
#awk -F'\t' '{print ">"$1"\n"$2}' all_refseq_complete_prodigal_single_meta.faa.tab  > all_refseq_complete_prodigal_single_meta.faa
#62396260 fasta headers

#rm temp?


#Genomes in refseq complete genomes have different prefixes. Grep'ing in a file containing all genomes (or their amino acids) takes a long time \
#Subset them according to their prefixes
#grep -o 'NZ_[A-Z][A-Z]' all_refseq_complete_prodigal_single_meta.faa.tab | sort | uniq > NZ_prefixes 
#grep '^NC_' all_refseq_complete_prodigal_single_meta.faa.tab > NC.faa.tab

#NZ_prefixes=$(cat NZ_prefixes)
#for prefix in $NZ_prefixes
#do
#    grep "^$prefix" all_refseq_complete_prodigal_single_meta.faa.tab > "$prefix".faa.tab
#done

#sed 's/\./CP/' NZ_CP.faa.tab > test
#for x in 10000 20000 30000 40000 50000 60000
#do
#    upper=$x
#    lower=$(echo $x | awk '{print $1-10000}')
#    echo $upper $lower
#    awk -v up=$upper -v low=$lower -F'CP' '{if ($2 < up && $2 >= low) print $0}' test > NZ_CP.$x.tab
#done
#sed -i 's/CP/\./2' NZ_CP.[0-9]*.tab

#Split up the CP prefixes, since there are MANY of those - for faster grep'ing later. 
#for x in 5000 10000 15000 20000 25000 30000 35000 40000 45000 50000 55000 60000
#do
#    upper=$x
#    lower=$(echo $x | awk '{print $1-5000}')
#    echo $upper $lower
#    awk -v up=$upper -v low=$lower -F'CP' '{if ($2 < up && $2 >= low) print $0}' test > NZ_CP.$x.tab
#done
#sed -i 's/CP/\./2' NZ_CP.[0-9]*.tab

#rm temp? all_refseq_complete_prodigal_single_meta.faa.tab prodigal_per_genus.faa prodigal_per_genus_meta.faa NZ_CP.faa.tab












###Start of analysis

#Blast all genomes against CARD database. 

$diamond blastp --query $database_dir/single_meta_prodigal/all_refseq_complete_prodigal_single_meta.faa \
--threads $numthreads --db $database_dir/protein_fasta_protein_homolog_model \
--max-target-seqs 1 --evalue 10E-10 --query-cover 80 --subject-cover 80 \
--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp \
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
        model=$(grep "$acc" $database_dir/aro_index.tsv | awk '{print $4"_"$3}')
        #echo $model
        #Get RGI cutoff values
        bitscore_cutoff=$(grep "$model" $database_dir/proteindb.fsa | awk '{print $7}')
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
    echo "Filtering DIAMOND hits $pid"
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
        taxid=$(grep -m 1 -w "$acc" $database_dir/nucl_gb.accession2taxid.subset | awk '{print $3}' )
        taxonomy=$(grep -m 1 -w "^$taxid" $database_dir/fullnamelineage_edit.dmp | sed -e 's/ /_/g' \
        | awk -F'|' '{print $3}' | awk -F';_' '{print $1";_"$2";_"$3";_"$4";_"$5";_"$6";_"$7}')
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
    echo "Finding taxonomy of filtered hits $pid"
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
    echo "Processing hits that pass RGI cutoffs $pid"
    wait $pid
done

cd ..
cat tempdir/tax.?? > rgi_pass.tax
rm -r tempdir diamond_hits_rgi_cutoff
#115268 rgi_pass.tax

paste diamond_hits_rgi_pass rgi_pass.tax > diamond_hits_rgi_pass.Rdat
paste diamond_hits_rgi_fail_ID_pass rgi_fail_ID_pass.tax > diamond_hits_rgi_fail_ID_pass.Rdat

#Hit regions might extend beyond the end, but only until the end is extracted. The end coordinate will still be wrong (larger than the total length). 
#If start coordinate is negative, replace it with 1
halfdist=12170

awk '{print $1}' CARD_refseq_filt.diamond | sed 's/>//' | awk -F'_' -v var="$halfdist"  '{print $1"_"$2" "$3-var" "$4+var}' \
| awk '{if ($2 < 0) print $1":1-"$3; else print $1":"$2"-"$3}' | sed '/\#/d' > ARG_coordinates.txt
awk '{print $2}' CARD_refseq_filt.diamond | awk -F'\|' '{print $3}' | sed -e 's/ARO://' -e '/^$/d' > AROS
awk '{print $1}' CARD_refseq_filt.diamond | awk -F'_' '{print $5}' > prID
wc -l ARG_coordinates.txt AROS prID
#176888 ARG_coordinates.txt
#176888 AROS
#176888 prID

xargs samtools faidx $database_dir/refseq_complete.fasta < ARG_coordinates.txt > temp0
grep -c '>' temp0
#176888

echo "Extracted fasta sequences of ARG regions passing filters"

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
            grep "$acc"  $database_dir/single_meta_prodigal/NC.faa.tab | sed 's/ /_/g' | \
            awk -F'_' -v start=$start_coord -v end=$end_coord -v locus=$i '{if ($3 >= start && $4 <= end) print locus":"$0}' >> temp.$num
        else
            prefix2=$(echo $i | awk -F'_' '{print $2}' | sed 's/[0-9]*\.[0-9]*//')
        if [ $prefix2 == "CP" ]
        then
            number=$(echo $i | awk -F'_' '{print $2}' | sed -e 's/CP//' -e 's/\.[0-9]*//' -e 's/^0*//')
            rnum=$(round5k $number)
            grep "$acc"  $database_dir/single_meta_prodigal/NZ_CP.$rnum.tab | sed 's/ /_/g' | \
            awk -F'_' -v start=$start_coord -v end=$end_coord -v locus=$i '{if ($3 >= start && $4 <= end) print locus":"$0}' >> temp.$num
        else
            grep "$acc"  $database_dir/single_meta_prodigal/NZ_$prefix2.faa.tab | sed 's/ /_/g' | \
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
    echo "Extracting amino acid subsets from ARG loci $pid"
    wait $pid
done

#Takes a while

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
        cacc=$(grep "$acc" $database_dir/refseq_complete_chr_headers.txt)
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
rm -f templist.?? pids CARD_refseq_filt.diamond
for l in l??
do
    number=$(echo $l | sed 's/l//')
    replicon_fun $l $number &
    echo $! >> pids
done


#Wait for all processes to finish
all_pids=$(cat pids)
for pid in $all_pids; do
    echo "Extracting plasmid/chromosome info $pid"
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
$diamond blastp --db $database_dir/IS \
--query $database_dir/single_meta_prodigal/all_refseq_complete_prodigal_single_meta.faa \
--outfmt 6 --threads 70 --evalue 10E-30 --query-cover 90 --max-target-seqs 1 > prodigal_all_ISfinder.blastp 
#843626 prodigal_all_ISfinder.blastp
#790628

grep '>' $database_dir/refseq_complete_plasmid_headers.txt | awk '{print $1}' | sed 's/>//' > all_refseq_plasmid_acc
grep '>' $database_dir/refseq_complete_chr_headers.txt | awk '{print $1}' | sed 's/>//' > all_refseq_chr_acc
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
            seqlength=$(grep -w "$acc" $database_dir/refseq_complete_seqlengths | awk '{print $2}')
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
            seqlength=$(grep -w "$acc" $database_dir/refseq_complete_seqlengths | awk '{print $2}')
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
    echo "Compiling information on ARG loci $pid"
    wait $pid
done

cat s??.temp > simulated_dist.txt
rm s??.temp s??
#13298 simulated_dist.txt


#End of simulated distance calculations (not used in current manuscript)




#Run blast against ISfinder db with loci_ARGs.faa as query
#ISfinder db modified to more easily get also IS family names
$diamond blastp --db $database_dir/IS_mod --query ARG_loci.faa \
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
    echo "Compiling more info on ARG loci $pid"
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


#Filtering on simulated ARG:IS distances is disabled since version 0.3

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
        dna_acc=$(grep -w "$aro" $database_dir/aro_index.tsv | sort -u | awk -F'\t' '{print $7}' | sort -u)
        #Get general resistance mechanism
        mech=$(grep -w "$dna_acc" $database_dir/aro_categories_index.tsv | awk -F'\t' '{print $5}' | sort -u | sed 's/ /_/g')
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
    echo "Extracting mode of action of predicted ARGs $pid"
    wait $pid
done

cat templist.?? | sort | uniq -c | sed 's/^ *//' > passing_regions_mech
rm -f templist.?? temparo.list x??

rm -f passing_regions.tab passing_regions.fasta passing_regions_headers ARG_regions_edit.fasta

#The below is not used in current version

#Get regions that fail the mean distance < expected by random criterium.
#Per ARO, count how many fail and how many pass. Divide them to get a ratio. A high ratio indicates that a given ARO is often placed in a "highly mobilized" context
#awk '{if ($8 < $9) print $0}' ARG_IS_distances.dat | awk '{print $2}' | sort | uniq -c | awk '{print $2,$1}' | sort -n > exp_fail
#awk '{if ($8 > $9) print $0}' ARG_IS_distances.dat | awk '{print $2}' | sort | uniq -c | awk '{print $2,$1}' | sort -n > exp_pass

#join -a 1 -a 2 -j 1 -e '0' -o '2.1 2.2 1.2' exp_fail exp_pass | awk '{print $0,$3/($2+$3)}' | sed '/^0/d' > highly_mobilized_aros
#1154 highly_mobilized_aros

#Cluster similar regions with usearch

usearch_fun () {
    aro_num2=$1
    grep -A1 "_$aro_num2$" ARG_regions_edit.fasta | grep -v '\-\-' > AROs_otus/$aro_num2.fas
    nohup $usearch -cluster_fast AROs_otus/$aro_num2.fas -centroids AROs_otus/$aro_num2.otus.fa -id 0.99 \
    -uc AROs_otus/$aro_num2.otus.uc -sizeout -minqt 1 -query_cov 0.9 -target_cov 0.9 -strand both -threads 30 -sort length -maxhits 1 >/dev/null 2>&1 &
}

sed -e 's/:/_/' -e 's/-/_/' -e 's/|/_/' ARG_regions.fasta > ARG_regions_edit.fasta
mkdir AROs_otus
#Start with the most abundant AROs, since these take longer to cluster
aro_list=$(grep '>' ARG_regions_edit.fasta | awk -F'_' '{print $6}' | sort | uniq -c | sort -nrk1,1 | awk '{print $2}' )

rm -f cluster_pids
for i in $aro_list
do
    ucount=$(pgrep 'usearch' | wc -l)
    #echo $ucount
    maxcount=$(($numthreads/2))
    if [ $ucount -lt $maxcount ]
    then
        #echo $ucount
        usearch_fun $i &
        echo $! >> cluster_pids
    else
        echo "We are running $ucount out of a total $maxcount USEARCH processes. Waiting for some to finish"
        sleep 30
        usearch_fun $i &
        echo $! >> cluster_pids
    fi
done


#Wait for all clustering to finish - takes some time. Only look for usearch processes created by you (by user id)
usid=$(whoami)
while [[ ${?} == 0 ]]      # Repeat until the process has terminated.
do
    pids=$(pgrep -u $usid 'usearch')
    echo "Waiting for USEARCH to finish process IDs:" $pids
    sleep 10s              
    ps -p $pids >/dev/null
done


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
        dna_acc=$(grep -w "$curr_acc" $database_dir/aro_index.tsv | sort -u | awk -F'\t' '{print $7}' | sort -u)
        #Get general resistance mechanism
        mech=$(grep -w "$dna_acc" $database_dir/aro_categories_index.tsv | awk -F'\t' '{print $5}' | sort -u | sed 's/ /_/g')
        #Get also the specific mechanism (e.g. beta-lactamase)
        mech2=$(grep -w "$dna_acc" $database_dir/aro_categories_index.tsv | awk -F'\t' '{print $3}' | sort -u | sed 's/ /_/g')
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
        suffix=$(head /dev/urandom | tr -dc A-Za-z0-9 | head -c 20 ; echo '')
        grep -w -A1 "$a" temp$num.fas > temp_int.fa.$suffix
        rm -rf integron_finder/"$a"
        icount=$(pgrep 'integron_finder' | wc -l)
    if [ $icount -lt 20 ]
    then
        $integron_finder --local_max --cpu 1 --outdir integron_finder/"$a" --linear --prodigal $prodigal temp_int.fa.$suffix >/dev/null 2>&1
    else
        sleep 10
        $integron_finder --local_max --cpu 1 --outdir integron_finder/"$a" --linear --prodigal $prodigal temp_int.fa.$suffix >/dev/null 2>&1
    fi
    rm -f temp_int.fa.$suffix
    done
}


rm -f pids_int
for i in x??
do
    num=$(echo $i | sed 's/x//')
    int_fun $i $num &
    echo $! >> pids_int
done

#Wait for all clustering to finish - takes some time
pids=$(cat pids_int)
while [[ ${?} == 0 ]]      # Repeat until the process has terminated.
do
    active=$(ps -p $pids | grep 'S+' | wc -l)
    active_pids=$(ps -p $pids | grep 'S+' | awk '{print $1}')
    echo "Waiting for integron_finder to finish $active out of $numthreads processes with these IDs:" $active_pids
    sleep 10s
    ps -p $pids >/dev/null
done


#Wait for integron_finder processes to finish
#all_pids=$(cat pids_int)
#for pid in $all_pids; do
#    echo "Waiting for integron_finder processes to finish $pid"
#    wait $pid
#done
#rm pids_int

### Takes quite a while to process all loci
#ll integron_finder/ | wc -l
#53896


#Integron hits
#Since integron_finder performs prodigal gene prediction on input DNA loci, overlapping loci can have multiple ARGs (that belong to other loci). \
#Get the ARG coordinates from loci_ARGs.faa that have the format:
#>NZ_CP036368.1_1_22388_3002547:NZ_CP036368.1_pr10_9619_10218
#From integron_finder, the format is: 
#NZ_CP036337_192836_218375_3000165_size1_20 start end

#Adjust by start and stop of gene coordinates by -1 to properly adjust gene coordinates with locus coordinates
find integron_finder/ -name 'temp_int.fa.integrons' -exec grep "protein" {} \; | grep -v 'intI' | awk '{print $3,$4-1,$5-1}' | sort -u > int_passengers
awk -F'_' '{print $0,$3}' int_passengers | awk '{print $1,$4+$2"_"$4+$3}' | sed 's/_size[0-9]*_[0-9]*//' > int_passengers2
#43504 int_passengers... For some reason this number has varied. If you don't see 43504 lines in file int_passengers, repeat the int_fun loop above. It may be advised to run in screen session. 

rm temp??.fas

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
    echo "Processing integron passenger genes $pid"
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
                        taxid=$(grep -m 1 -w "$acc" $database_dir/nucl_gb.accession2taxid.subset | awk '{print $3}' )
                        if [ -z "$taxid" ]
                        then
                            echo "NOTHING TO SEE HERE (No tax_id) $parent $aro $acc" >> missing_taxid
                        else
                            taxonomy=$(grep -m 1 -w "^$taxid" $database_dir/fullnamelineage_edit.dmp | sed -e 's/ /_/g' | awk -F'|' '{print $3}' | awk -F';_' '{print $1";_"$2";_"$3";_"$4";_"$5";_"$6";_"$7}')
                            #Include species
                            #taxonomy=$(grep -m 1 -w "^$taxid" $database_dir/fullnamelineage_edit.dmp | sed -e 's/ /_/g' | awk -F'|' '{print $3}' | awk -F';_' '{print $1";_"$2";_"$3";_"$4";_"$5";_"$6";_"$7";_"$8}')
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
    echo "Processing taxonomy of ARG loci $pid"
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
    acc=$(grep -m1 "$aro" $database_dir/protein_fasta_protein_homolog_model.fasta | awk -F'|' '{print $2}')
    nacc=$(grep -m1 "$aro" $database_dir/nucleotide_fasta_protein_homolog_model.fasta | awk -F'|' '{print $2}')
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
#16S_IS_distance				16S_control.sh (included in ARGMOB_dbs folder)
#plasmid_tax					all_acc_to_tax.sh (included in ARGMOB_dbs folder)
#chr_tax						all_acc_to_tax.sh (included in ARGMOB_dbs folder)
#refseq_1st_chr.tax				all_acc_to_tax.sh (included in ARGMOB_dbs folder)
#post_clust_tax_c				all_acc_to_tax.sh (included in ARGMOB_dbs folder)
#protein_homolog.tax			(included in ARGMOB_dbs folder). Made with the following command:
#fasta_formatter -i /DATA_1/tkn/ProxARG19/CARD_db/protein_fasta_protein_homolog_model.fasta -t \
#| grep -o "\[.*\]" | awk '{print $1}' | sed -e 's/\[//' -e 's/\]//' | sort | uniq -c | sed -e 's/^ *//' > protein_homolog.tax 

cp $database_dir/16S_IS_distance .
cp $database_dir/plasmid_tax .
cp $database_dir/chr_tax .
cp $database_dir/refseq_1st_chr.tax .
cp $database_dir/post_clust_tax_c .
cp $database_dir/protein_homolog.tax .
cp $database_dir/tnp_registry_supp2.csv .


#Prepare tar.gz file for transferring to local computer and analyses in R/RStudio
tar -czvf argmob_data.tar.gz ARG_IS_distances_mech.txt simulated_dist.txt \
aro_tax.Rdat aro_acc.txt diamond_hits_rgi_cutoff.pass diamond_hits_rgi_cutoff.fail \
ARG_integrons_ISProx.txt passing_regions_mech \
16S_IS_distance plasmid_tax chr_tax refseq_1st_chr.tax post_clust_tax_c \
protein_homolog.tax tnp_registry_supp2.csv rgi_pass.tax rgi_fail_ID_pass.tax

#The file argmob_data.tar.gz can be moved to another computer for importing into RStudio. 
#The provided Rmarkdown script runs further analyses, statistics and makes plots. 


###Important: in the start of the Rmarkdown script, there is a list of required R packages. Make sure to have those installed before running the script.
#If Knitting the Rmarkdown script, it is highly recommended to make a html output file, in order to get interactive tables and plots. 

md5sum argmob_data.tar.gz
#89285bcb2b9b036bb803339ced36d9aa  argmob_data.tar.gz


#Compress some stuff
tar -czvf arg_loci_with_proteins.tar.gz ARG_regions.fasta ARG_loci.faa
tar -czvf integron_finder.tar.gz integron_finder/
tar -czvf AROs_otus.tar.gz AROs_otus/
md5sum integron_finder.tar.gz
#c651c34ee7d342ba7c457b03cc3e1ed9  integron_finder.tar.gz
rm -r integron_finder
md5sum arg_loci_with_proteins.tar.gz
#9274dfa3b448fb4c2ec02f71ed40606f  arg_loci_with_proteins.tar.gz
md5sum AROs_otus.tar.gz
#979b9a34d9aac88208eac800ac1d2bb9  AROs_otus.tar.gz
rm -r AROs_otus

#Final cleanup - WAIT WITH THIS IF YOU NEED TO TROUBLESHOOT/INVESTIGATE FURTHER
rm -f pids int_passengers* cluster_list all_AROs_otus.tab \
all_AROs_otus.fasta locus_list passing_clusters ARG_IS_list \
ARG_IS_distances2.txt protein_IS_all.blastp2 locus_arg_list2 \
nohup.out passing_regions ARG_IS_distances_pass_exp.dat \
temp* protein_IS_all.blastp diamond_hits_rgi_cutoff \
diamond_hits_rgi_fail_ID_pass diamond_hits_rgi_pass diamond_hits_rgi_pass.Rdat \
diamond_hits_rgi_fail_ID_pass.Rdat ARG_coordinates.txt \
AROS ARG_regions.fasta ARG_loci.faa ARG_list \
locus_arg_list replicon_acc locus_arg_list_repl loci_ARGs.faa \
ARG_accs prodigal_all_ISfinder.blastp all_refseq_plasmid_acc \
all_refseq_chr_acc ARG_IS_distances.dat passing_regions_mech \
ARG_regions_edit.fasta ARG_IS_distances_clusters.dat \
all_AROs_otus.fas aro_acc.txt loci_ARGs_mod_headers
