grep '>' /database/ncbi/ncbi_complete_genomes_121219/refseq_complete_plasmid_headers.txt | awk '{print $1}' | sed 's/>//' > all_refseq_plasmid_acc
grep '>' /database/ncbi/ncbi_complete_genomes_121219/refseq_complete_chr_headers.txt | awk '{print $1}' | sed 's/>//' > all_refseq_chr_acc
plasmids=$(cat all_refseq_plasmid_acc)
rm -f plasmid_tax
split --number=l/100 -d all_refseq_plasmid_acc p

rm -f plasmid_tax
plasmid_fun () {
	list_in=$1
	plist=$(cat $list_in)
	for i in $plist
	do
	taxid=$(grep -m 1 -w "$i" /DATA_1/tkn/ProxARG19/nucl_gb.accession2taxid.subset | awk '{print $3}' )
	taxonomy=$(grep -m 1 -w "^$taxid" /DATA_1/tkn/ProxARG19/fullnamelineage_edit.dmp | sed -e 's/ /_/g' | awk -F'|' '{print $3}' | awk -F';_' '{print $1";_"$2";_"$3";_"$4";_"$5";_"$6";_"$7}')
	echo $i $taxonomy >> plasmid_tax
	done
}
for p in p??
do
plasmid_fun $p &
echo $! >> pids
done

all_pids=$(cat pids)
for pid in $all_pids; do
echo $pid
wait $pid
done

rm -f chr_tax

chr_fun () {
	list_in=$1
	clist=$(cat $list_in)
	for i in $clist
	do
	taxid=$(grep -m 1 -w "$i" /DATA_1/tkn/ProxARG19/nucl_gb.accession2taxid.subset | awk '{print $3}' )
	taxonomy=$(grep -m 1 -w "^$taxid" /DATA_1/tkn/ProxARG19/fullnamelineage_edit.dmp | sed -e 's/ /_/g' | awk -F'|' '{print $3}' | awk -F';_' '{print $1";_"$2";_"$3";_"$4";_"$5";_"$6";_"$7}')
	echo $i $taxonomy >> chr_tax
done
}
split --number=l/100 -d all_refseq_chr_acc c

rm -f pids
for c in c??
do
chr_fun $c &
echo $! >> pids
done

all_pids=$(cat pids)
for pid in $all_pids; do
echo $pid
wait $pid
done


rm p?? c??

sed -i -e 's/ /\t/g' -e 's/;_/\t/g' chr_tax
sed -i -e 's/ /\t/g' -e 's/;_/\t/g' plasmid_tax



#Tax not divided into chr and plasmids. Some strains have MANY plasmids and some may have multiple chromosomes. 
#accessions=$(ls /database/ncbi/ncbi_complete_genomes_121219/refseq/bacteria)
#rm -f all_refseq_acc.tax
#for i in $accessions
#do
#ucount=$(pgrep -u tkn 'perl' | wc -l)
#echo $ucount
#if [ $ucount -lt 10 ]
#then
#echo $ucount
#./esearch_multi.sh $i &
#else
#echo $ucount
#echo "Waiting"
#sleep 10
#fi
#done


acc_function () {
	acc_list=$1
	number=$2
	accessions=$(cat $acc_list)
	rm -f refseq_1st_chr.$number
	for a in $accessions
	do
	#echo $a
	access=$(fasta_formatter -t -i /database/ncbi/ncbi_complete_genomes_121219/refseq/bacteria/"$a"/*fna | awk -v acc=$i '{ print length($NF)"\t"$1"\t"}' | sort -nr | head -n1 | awk '{print $2}')
	taxid=$(grep -m 1 -w "$access" /DATA_1/tkn/ProxARG19/nucl_gb.accession2taxid.subset | awk '{print $3}' )
	taxonomy=$(grep -m 1 -w "^$taxid" /DATA_1/tkn/ProxARG19/fullnamelineage_edit.dmp | sed -e 's/ /_/g' | awk -F'|' '{print $3}' | awk -F';_' '{print $1";_"$2";_"$3";_"$4";_"$5";_"$6";_"$7}')
	echo $access $taxid $taxonomy >> refseq_1st_chr.$number
	done
}

rm -f pids
ls /database/ncbi/ncbi_complete_genomes_121219/refseq/bacteria > all_accessions
split --number=l/100 -d all_accessions a
for acc_list in a??
do
num=$(echo $acc_list | sed 's/a//')
acc_function $acc_list $num &
echo $! >> pids
done

all_pids=$(cat pids)
for pid in $all_pids; do
echo $pid
wait $pid
done


cat refseq_1st_chr.?? > refseq_1st_chr
rm refseq_1st_chr.??
sed -e 's/ /\t/g' -e 's/;_/\t/g' refseq_1st_chr > refseq_1st_chr.tax
rm a?? 


#Get tax after clustering
otu_acc=$(grep '>' all_AROs_otus.fasta | sed -e 's/>//' | awk -F '_' '{print $1"_"$2}')

grep '>' all_AROs_otus.fasta | sed -e 's/>//' | awk -F '_' '{print $1"_"$2}' > otu_acc
split --number=l/100 -d otu_acc o

rm -f post_clust_tax_*

otu_fun () {
	list_in=$1
	olist=$(cat $list_in)
	for i in $olist
	do
	acc_p=$(grep -m 1 "$i" all_refseq_plasmid_acc)
	if [ -z "$acc_p" ]
	then
	acc=$(grep -m 1 "$i" all_refseq_chr_acc)
	taxid=$(grep -m 1 -w "$acc" /DATA_1/tkn/ProxARG19/nucl_gb.accession2taxid.subset | awk '{print $3}' )
	taxonomy=$(grep -m 1 -w "^$taxid" /DATA_1/tkn/ProxARG19/fullnamelineage_edit.dmp | sed -e 's/ /_/g' | awk -F'|' '{print $3}' | awk -F';_' '{print $1";_"$2";_"$3";_"$4";_"$5";_"$6";_"$7}')
	echo $acc $taxonomy >> post_clust_tax_c
	else
	acc=$acc_p
	taxid=$(grep -m 1 -w "$acc" /DATA_1/tkn/ProxARG19/nucl_gb.accession2taxid.subset | awk '{print $3}' )
	taxonomy=$(grep -m 1 -w "^$taxid" /DATA_1/tkn/ProxARG19/fullnamelineage_edit.dmp | sed -e 's/ /_/g' | awk -F'|' '{print $3}' | awk -F';_' '{print $1";_"$2";_"$3";_"$4";_"$5";_"$6";_"$7}')
	echo $acc $taxonomy >> post_clust_tax_p
	fi
	#echo $acc
done
}
rm -f pids
for o in o??
do
otu_fun $o &
echo $! >> pids
done

all_pids=$(cat pids)
for pid in $all_pids; do
echo $pid
wait $pid
done


sed -i -e 's/ /\t/g' -e 's/;_/\t/g' post_clust*
sed -i 's/\tdelta\/epsilon_subdivisions//' refseq_1st_chr.tax
rm o??

#Get taxonomy of CARD
fasta_formatter -i /DATA_1/tkn/ProxARG19/CARD_db/protein_fasta_protein_homolog_model.fasta -t | grep -o "\[.*\]" | awk '{print $1}' | sed -e 's/\[//' -e 's/\]//' | sort | uniq -c | sed -e 's/^ *//' > protein_homolog.tax
