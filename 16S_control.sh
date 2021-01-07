#Barrnap was run on all complete refseq genomes to find 16S rRNA genes
#/database/ncbi/ncbi_complete_genomes_121219
#barrnap /database/ncbi/ncbi_complete_genomes_121219/refseq_complete.fasta --threads 50 > refseq_complete_barrnap.txt

#Get 16S sequences
grep '16S_rRNA' refseq_complete_barrnap.txt > 16S_barrnap.txt
#Extract 100kbp up-and downstream of 16S genes
awk '{print $1,$4-100000,$5+100000}' 16S_barrnap.txt | awk '{if ($2 < 0) print $1,1,$3;else print $0}' | sed -e 's/ /:/' -e 's/ /-/' > 16S_coordinates
xargs samtools faidx /database/ncbi/ncbi_complete_genomes_121219/refseq_complete.fasta < 16S_coordinates > 16S_regions.fasta
fasta_formatter -i 16S_regions.fasta -o 16S_regions.tab -t

rm x??
split -l 400 16S_regions.tab

#Do prodigal gene prediction on extracted regions
rm -f x??.faa
for i in x??
do
awk '{print ">"$1"\n"$2}' $i > $i.2
nohup prodigal -a $i.faa -i $i.2 -p meta -q &
done

cat x??.faa > 16S_regions_all.faa

#Find IS elements in extracted regions
nohup diamond blastp --db /opt/prokka/db/kingdom/Bacteria/IS --query 16S_regions_all.faa --outfmt 6 \
--threads 70 --evalue 10E-30 --max-target-seqs 1 --query-cover 90 > 16S_protein_IS.blastp &

rm x??
rm x??.faa
rm x??.2
#Sort the blastn results by query sequence, then bitscore, e-value and finally % ID
awk '{print $1"\t"$0}' 16S_protein_IS.blastp | sed -e 's/_\([0-9]*\)\t/\t\1\t/' | sort -u -k1,1 > sortedFile

#Keep only best hit per query
sort -u -k1,1 sortedFile | tail +4 | grep -v 'diamond' > uniqFile
awk '{print $3}' uniqFile > IS_prots.list

#Get a list of proteins with IS blastp matches 
seqtk subseq 16S_regions_all.faa IS_prots.list > IS_prots.faa

#Get only the coordinates of the IS elements. Protein coordinates are relative to the subregion. Adjust these using the start-end coordinates of the region. 
fasta_formatter -i IS_prots.faa -t | awk '{print $1,$1,$3,$5}' | sed -e 's/[A-Z]*_[A0-Z9]*.[0-9]*://' -e 's/-[0-9]*_[0-9]*//' | awk '{print $2,$3+$1,$4+$1}' | sort -nk1,1 > IS_coordinates
#get the IS annotations
sort -nk1,1 uniqFile | awk '{print $1,$2}' > IS_annotations
#Merge annotations and coordinates
paste IS_coordinates IS_annotations > IS_coord_annot

#Get the original 16S coordinates for each region 
awk '{print $1,$4-100000,$5+100000,$4,$5}' 16S_barrnap.txt | awk '{if ($2 < 0) print $1,1,$3,$4,$5;else print $0}' |  sed -e 's/ /:/' -e 's/ /-/' > 16S_regions_orig_coord

split --number=l/100 -d 16S_regions_orig_coord
for split_region in x{00..99}
do
number=$(echo $split_region | sed 's/x//')
regions=$(awk '{print $1}' $split_region)
rm -f 16S_IS_coords.$number
nohup ./processing_IS_coords.sh $number &
done

cat 16S_IS_coords.?? > 16S_IS_coords
rm 16S_IS_coords.??
rm x??

#rm -f 16S_IS_coords
#for i in $regions
#do
#orig_coord=$(grep "$i" 16S_regions_orig_coord | awk '{print $2,$3}')
#grep "$i" IS_coord_annot | awk -v orig="$orig_coord" '{print $1,$2,$3,$5,orig}' >> 16S_IS_coords
#done

#Get the up-and down-distance to closest IS elements
awk '{if ($3 < $5) print $0,$5-$3; else print $0,$2-$6}' 16S_IS_coords > 16S_IS_distance
#Header=
#region	IS_start	IS_end	IS_element	16S_start	16S_end	distance(up or down)

regions_with_IS=$(awk '{print $1}' 16S_IS_coords | sed 's/_[0-9]*$//')

#How many 16S regions with IS elements
awk '{print $1}' 16S_IS_coords | sed 's/_[0-9]*$//' | sort -u | wc -l
#40969
#How many 16S regions in total?
wc -l 16S_regions_orig_coord
#80143
#16S regions without IS:
#80143 - 40969 = 39174








#Identify rpoB genes. Downloading all gff annotations. 

