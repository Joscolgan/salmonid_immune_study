while read line;
do
pike="$(echo "$line" | cut -f 1 )";
trout="$(echo "$line" | cut -f 2 )";
salmon="$(echo "$line" | cut -f 3)";
rainbow="$(echo "$line" | cut -f 4)";
## Pull out sequence headers for pike:
echo "$pike"
zgrep "$pike" ../data/cds/*.gz | cut -d ':' -f 2- | sed 's/>//g' > "$trout"_tmp;
longest_transcript="$(seqtk subseq <(zcat ../data/cds/*.gz) "$trout"_tmp | seqtk comp | sort -k2,2nr | cut -f 1 | head -n 1 )";
echo "$longest_transcript"
grep "$longest_transcript" <(zcat ../data/cds/*.gz) | sed 's/>//g' >> "$trout"_longest_seq_div.txt;
## Pull out sequence headers for trout:
echo "$trout"
zgrep "$trout" ../data/cds/* | cut -d ':' -f 2- | sed 's/>//g' > "$trout"_tmp;
longest_transcript="$(seqtk subseq <(zcat ../data/cds/*.gz) "$trout"_tmp | seqtk comp | sort -k2,2nr | cut -f 1 | head -n 1 )";
grep "$longest_transcript" <(zcat ../data/cds/*.gz) | sed 's/>//g' >> "$trout"_longest_seq_div.txt;
## Pull out sequence headers for salmon:
echo "$salmon"
zgrep "$salmon" ../data/cds/* | cut -d ':' -f 2- | sed 's/>//g' > "$trout"_tmp;
longest_transcript="$(seqtk subseq <(zcat ../data/cds/*.gz) "$trout"_tmp  | seqtk comp | sort -k2,2nr | cut -f 1 | head -n 1 )";
grep "$longest_transcript" <(zcat ../data/cds/*.gz) | sed 's/>//g' >> "$trout"_longest_seq_div.txt;
## Pull out sequence headers for rainbow trout:
echo "$rainbow"
zgrep "$rainbow" ../data/cds/* | cut -d ':' -f 2- | sed 's/>//g' > "$trout"_tmp;
longest_transcript="$(seqtk subseq <(zcat ../data/cds/*.gz) "$trout"_tmp | seqtk comp | sort -k2,2nr | cut -f 1 | head -n 1 )";
grep "$longest_transcript" <(zcat ../data/cds/*.gz) | sed 's/>//g' >> "$trout"_longest_seq_div.txt;
seqtk subseq <(zcat ../data/cds/*.gz) "$trout"_longest_seq_div.txt -l 80  >> "$trout"_orthologs_div.fasta;
echo "Complete for" "$trout";
done < test.txt
