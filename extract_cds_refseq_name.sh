

for name in *refseq.txt;
do new_name="$(echo "$name" | cut -d '.' -f 1 )";
echo "$new_name";
while read line;
do
zgrep "$line" ~/projects/2020-04-11_salmonid_immune_project/2020-04-11_orthofinder_of_immune_genes/data/2020-04-11_refseq_datasets/*/*gff* | \
awk '$3=="CDS"' | head -n 1 | cut -f 9 | cut -d ';' -f 2 | \
sed 's/Parent=rna-//g' >> "$new_name"_transcript_names.txt;
done < "$name";
done
