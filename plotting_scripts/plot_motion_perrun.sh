#! /bin/bash

subjIDs=("RR" "MM" "PP" "MK" "AB" "AD" "LA" "AE" "TP" "NM" "AF"  \
	"AG" "AH" "AI" "GG" "SL" "UV" "PQ" "KQ" "LN" "RT" "PT"  \
	"PL" "NS")

cd /projectnb/somerslab/tom/projects/spacetime_network/data/unpacked_data_nii

for ss in ${subjIDs[@]}
do
	plot-twf-sess -mc -s $ss -d ./ -fsd bold
done
