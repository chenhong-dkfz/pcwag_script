 
cohs="PBCA-DE PACA-AU PACA-CA PAEN-AU SARC-US OV-AU OV-US BRCA-EU BRCA-UK PRAD-UK PRAD-US EOPC-DE" #applied cohorts
els="enhancer promoter"

for coh in ${cohs};do
  for el in ${els};do
    python  /home/hongc/CaSINo/code/analysis_scripts/fe2cgcmap.py ${el} ${coh}
  done
done