#!/bin/sh
#############################################################
#
# Author: Joe Colgan
#
# Program:
#
# Purpose:
#
#############################################################

## Extract coding sequence for longest transcript for putative diverged gene:
while read line;
do
pike="$(echo "$line" | cut -f 1 | cut -d '.' -f 1 )";
trout="$(echo "$line" | cut -f 2 | cut -d '.' -f 1 )";
echo "$pike";
echo "$trout";
zgrep "$pike" ../data/cds/*.gz | cut -d ':' -f 2- | sed 's/>//g' > "$trout"_tmp;
longest_transcript="$(seqtk subseq <(zcat ../data/cds/*.gz) "$trout"_tmp | seqtk comp | sort -k2,2nr | cut -f 1 | head -n 1 )";
echo "$longest_transcript";
grep "$longest_transcript" <(zcat ../data/cds/*.gz) | sed 's/>//g' >> "$trout"_longest_seq_div.txt;
echo "$trout";
zgrep "$trout" ../data/cds/* | cut -d ':' -f 2- | sed 's/>//g' > "$trout"_tmp;
longest_transcript="$(seqtk subseq <(zcat ../data/cds/*.gz) "$trout"_tmp | seqtk comp | sort -k2,2nr | cut -f 1 | head -n 1 )";
grep "$longest_transcript" <(zcat ../data/cds/*.gz) | sed 's/>//g' >> "$trout"_longest_seq_div.txt;
seqtk subseq <(zcat ../data/cds/*.gz) "$trout"_longest_seq_div.txt -l 80  >> "$trout"_orthologs_div.fasta;
echo "Complete for" "$trout";
done < diverged_genes.txt

## Extract coding sequence for longest transcript for putative conserved gene:
while read line;
do
pike="$(echo "$line" | cut -f 1 | cut -d '.' -f 1 )";
trout="$(echo "$line" | cut -f 2 | cut -d '.' -f 1 )";
echo "$pike";
echo "$trout";
zgrep "$pike" ../data/cds/*.gz | cut -d ':' -f 2- | sed 's/>//g' > "$trout"_tmp;
longest_transcript="$(seqtk subseq <(zcat ../data/cds/*.gz) "$trout"_tmp | seqtk comp | sort -k2,2nr | cut -f 1 | head -n 1 )";
echo "$longest_transcript";
grep "$longest_transcript" <(zcat ../data/cds/*.gz) | sed 's/>//g' >> "$trout"_longest_seq_con.txt;
echo "$trout";
zgrep "$trout" ../data/cds/* | cut -d ':' -f 2- | sed 's/>//g' > "$trout"_tmp;
longest_transcript="$(seqtk subseq <(zcat ../data/cds/*.gz) "$trout"_tmp | seqtk comp | sort -k2,2nr | cut -f 1 | head -n 1 )";
grep "$longest_transcript" <(zcat ../data/cds/*.gz) | sed 's/>//g' >> "$trout"_longest_seq_con.txt;
seqtk subseq <(zcat ../data/cds/*.gz) "$trout"_longest_seq_con.txt -l 80  >> "$trout"_orthologs_con.fasta;
echo "Complete for" "$trout";
done < converged_genes.txt

## Extract protein sequence for longest transcript for putative conserved gene:
while read line;
do
pike="$(echo "$line" | cut -f 1 | cut -d '.' -f 1 )";
trout="$(echo "$line" | cut -f 2 | cut -d '.' -f 1 )";
echo "$pike";
echo "$trout";
zgrep "$pike" ../data/cds/*.gz | cut -d ':' -f 2- | sed 's/>//g' > "$trout"_tmp;
longest_transcript="$(seqtk subseq <(zcat ../data/pep/*.gz) "$trout"_tmp | seqtk comp | sort -k2,2nr | cut -f 1 | head -n 1 )";
echo "$longest_transcript";
grep "$longest_transcript" <(zcat ../data/pep/*.gz) | sed 's/>//g' >> "$trout"_longest_seq_div_pep.txt;
echo "$trout";
zgrep "$trout" ../data/pep/* | cut -d ':' -f 2- | sed 's/>//g' > "$trout"_tmp;
longest_transcript="$(seqtk subseq <(zcat ../data/pep/*.gz) "$trout"_tmp | seqtk comp | sort -k2,2nr | cut -f 1 | head -n 1 )";
grep "$longest_transcript" <(zcat ../data/pep/*.gz) | sed 's/>//g' >> "$trout"_longest_seq_div_pep.txt;
seqtk subseq <(zcat ../data/pep/*.gz) "$trout"_longest_seq_div_pep.txt -l 80  >> "$trout"_orthologs_div.faa;
echo "Complete for" "$trout";
done < diverged_genes.txt

## Extract protein sequence for longest transcript for putative conserved gene:
while read line;
do
pike="$(echo "$line" | cut -f 1 | cut -d '.' -f 1 )";
trout="$(echo "$line" | cut -f 2 | cut -d '.' -f 1 )";
echo "$pike";
echo "$trout";
zgrep "$pike" ../data/cds/*.gz | cut -d ':' -f 2- | sed 's/>//g' > "$trout"_tmp;
longest_transcript="$(seqtk subseq <(zcat ../data/pep/*.gz) "$trout"_tmp | seqtk comp | sort -k2,2nr | cut -f 1 | head -n 1 )";
echo "$longest_transcript";
grep "$longest_transcript" <(zcat ../data/pep/*.gz) | sed 's/>//g' >> "$trout"_longest_seq_con_pep.txt;
echo "$trout";
zgrep "$trout" ../data/pep/* | cut -d ':' -f 2- | sed 's/>//g' > "$trout"_tmp;
longest_transcript="$(seqtk subseq <(zcat ../data/pep/*.gz) "$trout"_tmp | seqtk comp | sort -k2,2nr | cut -f 1 | head -n 1 )";
grep "$longest_transcript" <(zcat ../data/pep/*.gz) | sed 's/>//g' >> "$trout"_longest_seq_con_pep.txt;
seqtk subseq <(zcat ../data/pep/*.gz) "$trout"_longest_seq_con_pep.txt -l 80  >> "$trout"_orthologs_con.faa;
echo "Complete for" "$trout";
done < converged_genes.txt

for name in *tmp;
do
new_name="$(echo "$name" | cut -d '_' -f 1 )";
mkdir "$new_name";
cd "$new_name";
mv ../"$new_name"_orthologs*.faa .;
mv ../"$new_name"_orthologs*.fasta .;
cp -p ../../dnds/codeml.ctl .;
cd ../;
done

## Align protein sequences first and then generate codon-based alignments: 
for name in ENSSTUG000000*/;
do
new_name="$(echo "$name" | cut -d '/' -f 1)";
echo "$new_name";
cd "$new_name";
clustalo -i "$new_name"_orthologs*.faa -o "$new_name"_aln.faa --force; pal2nal.pl "$new_name"_aln.faa "$new_name"_orthologs*.fasta -output paml -nogap > "$new_name".pal2nal;
tmp="$(echo "$new_name".pal2nal)";
sed -i "s/cluster_1.pal2nal/$tmp/g" codeml.ctl; codeml;
python2.7 ../../dnds/parse_codeml_output.py codeml.txt >> "$new_name"_result.txt; cd ../;
done

## Concatenate output:
while read line;
do
div_copy="$(echo "$line" | awk '{ print $4 }' | cut -d '.' -f 1 )";
con_copy="$(echo "$line" | awk '{ print $5 }' | cut -d '.' -f 1 )";
paste <(cat "$div_copy"/*result.txt | tail -n 1) <(cat "$con_copy"/*result.txt | tail -n 1) >> combined_div_con_dnds.txt;
done < ../data/ancestral_copy_plus_protein_length.txt 
