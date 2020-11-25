
## Extract gene name and distance matrix values for each brown trout ohnolog:
for name in *distmat;
do
new_name="$(echo "$name" | cut -d '.' -f 1)";
test_1="$(grep 'ENSSTUG' "$name" | rev | cut -f 1 | rev | cut -d ' ' -f 2 | head -n 1)";
gene_1="$(grep 'ENSSTUG' "$name" | rev | cut -f 1 | rev | cut -d ' ' -f 1 | head -n 1)";
test_2="$(grep 'ENSSTUG' "$name" | rev | cut -f 1 | rev | cut -d ' ' -f 2 | tail -n 1)";
gene_2="$(grep 'ENSSTUG' "$name" | rev | cut -f 1 | rev | cut -d ' ' -f 1 | tail -n 1)";
echo "$gene_1" "$gene_2" >> "$new_name"_distmat_trout.txt;
grep 'ENSE' "$name" | cut -f "$test_1","$test_2" >> "$new_name"_distmat_trout.txt;
done
