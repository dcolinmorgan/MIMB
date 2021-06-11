#!/bin/bash
#bash MIMB/slop_mega.sh

cut -f1,2,3,4,5 uniq_mega_BS_Ch0.txt|uniq >mega_motif.txt

../../../../udd/rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a mega_motif.txt -b hg38_tss_coord.txt > mega_motif_refseq.txt

../../../../udd/rekrg/Tools/bedtools2/bin/bedtools slop -i mega_motif_refseq.txt -g ../../../../udd/rekrg/Tools/bedtools2/genomes/human.hg38.genome -r 100 -l 100  > slop_mega_motif_refseq.txt

 # cut -f1,2,4,5,9 slop_mega_motif_refseq.txt |sort -k1,1 -k2,2n| uniq > slop_mega_motif_refseq01.txt

../../../../udd/rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a slop_mega_motif_refseq.txt -b GM12878_BS_Ch01.txt > slop_mega_BS_Ch0.txt

cut -f1,2,3,4,5,9,13,14,15,16 slop_mega_BS_Ch0.txt |sed -e "s/(//" -e "s/)//" > slop_mega_BS_Ch1.txt 

source /proj/relibs/relib00/conda/bin/activate
source activate mypy3 ## 
chmod +x MIMB/take_uniq.py

python MIMB/take_uniq.py
