##domain scale auc bestThres sens spec
sort -nk3 ROCminimumError.txt | awk '($4!="is"){print $1" "$2" "$4}' > Error_motif_scale_threshold.txt
sort -nk3 ROCmaxCorr.txt | uniq | awk '($4!="threshold"&& $4>0.6){print $1" "$2" "$4}' > Corr_Domain_scale_threshold.txt

join <(sort Corr_Domain_scale_threshold.txt) <(sort Error_motif_scale_threshold.txt ) | awk '{print $1" "$2" "$3" "$5}' >ErrorCorr.joined

##corrAccum	#poi togliere a mano i "+" e sostituirle con " "
cp ../../results/candidateGutfreundroc_domainsWithScale0.1.txt .
cp ../../results/candidateGutfreundroc_domainsWithScale0.05.txt .

#maxcorr	#poi togliere a mano i "+" e sostituirle con " "
cp ../../results/candidateLysateBil0MaxCorr_domainsWithScale0.1.txt .
cp ../../results/candidateLysateBil0MaxCorr_domainsWithScale0.05.txt .


grep -f candidateLysateBil0MaxCorr_domainsWithScale0.05.txt ErrorCorr.joined > FinalDomainScaleMaxCorr_CandidateLysateBil00.05.txt
grep -f candidateLysateBil0MaxCorr_domainsWithScale0.1.txt ErrorCorr.joined > FinalDomainScaleMaxCorr_CandidateLysateBil00.1.txt

grep -f candidateLysateBil0MaxCorr_domainsWithScale0.05.txt ErrorCorr.joined > FinalDomainScaleAccumCorr_CandidateLysateBil00.05.txt
grep -f candidateLysateBil0MaxCorr_domainsWithScale0.1.txt ErrorCorr.joined > FinalDomainScaleAccumCorr_CandidateLysateBil00.1.txt

#candidateGutfreundMaxCorr_domainsWithScale0.05.txt
#candidateGutfreundroc_domainsWithScale0.05.txt
#candidateLysateBil0MaxCorr_domainsWithScale0.05.txt

#candidateGutfreundMaxCorr_domainsWithScale0.1.txt
#candidateGutfreundroc_domainsWithScale0.1.txt
#candidateLysateBil0MaxCorr_domainsWithScale0.1.txt

