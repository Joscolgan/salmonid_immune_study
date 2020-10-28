
## Extract single pike matches to putative brown trout immune ohnologs:
tail -n +2 pike_comp_sub_one2many_immune_ohno_strutta.txt | cut -f 1,5 > pike_comp_sub_one2many_immune_ohno_strutta_gene_list.txt 

## Create a tmp folder:
mkdir dir
cd tmp

## Copy file to current directory:
cp -p ../pike_comp_sub_one2many_immune_ohno_strutta_gene_list.txt .

## Split file every two files to achieve brown trout ohnolog gene name and single pike ortholog: 
split -l 2 pike_comp_sub_one2many_immune_ohno_strutta_gene_list.txt 

## Delete original file to clean directory:
rm pike_comp_sub_one2many_immune_ohno_strutta_gene_list.txt 

## Create a single unique list of genes of interest:
for name in *;
do
sed 's/\t/\n/g' "$name" | sort | uniq > "$name"_sorted.txt;
done

## Create a single file of longest protein isoform for pike and trout genes of interest:
cd ../../data/
mkdir pike_longest_protein_isoform
mkdir trout_longest_protein_isoform

cd pike_longest_protein_isoform/
ln -s ~/projects/2020-04-11_salmonid_immune_project/2020-04-11_orthofinder_of_immune_genes/input/primary_transcripts/E_lucius.fa .

cd ../trout_longest_protein_isoform/
ln -s ~/projects/2020-04-11_salmonid_immune_project/2020-04-11_orthofinder_of_immune_genes/input/primary_transcripts/S_trutta.fa .

## Go back to results directory:
cd ../../results/tmp

for name in *sorted.txt;
do
new_name="$(echo "$name" | cut -d '.' -f 1 )";
while read line;
do
grep "$line" ../../data/pike_longest_protein_isoform/E_lucius.fa | sed 's/>//g' >> "$new_name"_protein_list.tmp;
grep "$line" ../../data/trout_longest_protein_isoform/S_trutta.fa | sed 's/>//g' >> "$new_name"_protein_list.tmp;
done < "$name"; done 

cat ../../data/pike_longest_protein_isoform/E_lucius.fa ../../data/trout_longest_protein_isoform/S_trutta.fa > pike_trout_longest_protein_isoform_combined.fa

for name in *tmp;
do
new_name="$(echo "$name" | cut -d '.' -f 1 )";
seqtk subseq ./pike_trout_longest_protein_isoform_combined.fa "$name" -l 80 >> "$new_name".fa;
done

## If gene list does not contain three genes, remove:
grep -c '>' *fa | sed 's/:/\t/g' | awk '$2!=3' | cut -f 1 | tail -n +2 > remove_list.txt
while read line; do rm "$line"; done < remove_list.txt

## Align using mafft:
for name in x*fa; do new_name="$(echo "$name" | cut -d '.' -f 1 )"; mafft --auto "$name" > "$new_name"_output; done

## Calculate evolutionary distance:
for name in *output; do new_name="$(echo "$name" | cut -d '.' -f 1 )"; ~/src/EMBOSS-6.6.0/emboss/distmat "$name" -protmethod 2 -outfile "$new_name".phy; done

for name in *output; do tmp="$(grep '>' "$name" | tr '\n' '\t' | sed 's/>//g')"; echo "$tmp" >> ohnologs_proitein_match.txt; done
 
grep 'SEL' *.phy  | cut -f 1,3,4,6 > combined_distmat.txt

paste combined_distmat.txt ohnologs_proitein_match.txt > combined_distmat_with_protein_names.txt

awk -F $'\t' '{ if ($2 > $3) print $2,$3,$5,$6,$7;}' combined_distmat_with_protein_names.txt >> ancestral_copy.txt
awk -F $'\t' '{ if ($2 < $3) print $3,$2,$5,$7,$6;}' combined_distmat_with_protein_names.txt >> ancestral_copy.txt 
awk -F $'\t' '{ if ($2 == $3) print $3,$2,$5,$7,$6;}' combined_distmat_with_protein_names.txt >> ancestral_copy.txt 
 
for name in *list.fa; do seqtk comp "$name" > "$name"_comp.txt; done

for name in *list.fa_comp.txt; do new="$(cut -f 2 "$name" | tr '\n' '\t')"; echo "$new" >> total_protein_comp.txt; done

paste ancestral_copy.txt total_protein_comp.txt > ancestral_copy_plus_protein_length.txt
mv ancestral_copy_plus_protein_length.txt ../
