#remove the read number added earlier --yy5
awk 'BEGIN{OFS="\t"}{if($1 ~ /^@/) {print} else {$2=and($2,compl(2048)); print substr($1,2),$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}}'