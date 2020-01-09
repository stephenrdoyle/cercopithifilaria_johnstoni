
# Repeat model and mask of filarial genomes

```shell
cd /nfs/users/nfs_s/sd21/lustre118_link/kirsty
```
### get genomes from WBP
```
mkdir genomes
cd genomes

wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/dirofilaria_immitis/PRJEB1797/dirofilaria_immitis.PRJEB1797.WBPS14.genomic.fa.gz -O DM.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/brugia_timori/PRJEB4663/brugia_timori.PRJEB4663.WBPS14.genomic.fa.gz -O BT.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/brugia_pahangi/PRJEB497/brugia_pahangi.PRJEB497.WBPS14.genomic.fa.gz -O BP.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/brugia_malayi/PRJNA10729/brugia_malayi.PRJNA10729.WBPS14.genomic.fa.gz -O BM.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/loa_loa/PRJNA246086/loa_loa.PRJNA246086.WBPS14.genomic.fa.gz -O LL.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/litomosoides_sigmodontis/PRJEB3075/litomosoides_sigmodontis.PRJEB3075.WBPS14.genomic.fa.gz -O LS.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/onchocerca_ochengi/PRJEB1204/onchocerca_ochengi.PRJEB1204.WBPS14.genomic.fa.gz -O OO.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/onchocerca_volvulus/PRJEB513/onchocerca_volvulus.PRJEB513.WBPS14.genomic.fa.gz -O OV.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/onchocerca_flexuosa/PRJNA230512/onchocerca_flexuosa.PRJNA230512.WBPS14.genomic.fa.gz -O OF.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/wuchereria_bancrofti/PRJNA275548/wuchereria_bancrofti.PRJNA275548.WBPS14.genomic.fa.gz -O WB.fa.gz


for i in *.gz; do gunzip ${i}; done

#--- kirsty sent C.johnstoni via WeTransfer - scp'ed into this dir and called it "CJ.fa"

cd ..
```


### make databases
```
for i in $( cd genomes/ ; ls -1 | sed 's/.fa//g' ); do \
mkdir ${i}_RM_OUT ; cd ${i}_RM_OUT ; bsub.py 1 01_RM_builddb "/nfs/users/nfs_s/sd21/lustre118_link/software/REPEATMASKER/RepeatModeler-open-1.0.11/BuildDatabase -name ${i} ../genomes/${i}.fa"  ; cd .. ;
done
```


### run modeller
```
for i in $( cd genomes/ ; ls -1 | sed 's/.fa//g' ); do \
cd ${i}_RM_OUT ; bsub.py --threads 20 10 02_RM_model  "/nfs/users/nfs_s/sd21/lustre118_link/software/REPEATMASKER/RepeatModeler-open-1.0.11/RepeatModeler -pa 20 -engine ncbi -database ${i}"; cd .. ;
done
```


### run masker
```
for i in $( cd genomes/ ; ls -1 | sed 's/.fa//g' ); do \
cd ${i}_RM_OUT ; bsub.py --threads 7 10 03_RM_mask  "/nfs/users/nfs_s/sd21/lustre118_link/software/REPEATMASKER/RepeatMasker/RepeatMasker -e ncbi -pa 7 -s -dir ./ -small -gff -lib RM_*/consensi.fa.classified ../genomes/${i}.fa"; cd .. ;
done
```

### collate output
```
> summary.data
for i in $( cd genomes/ ; ls -1 | sed 's/.fa//g' ); do \
echo ${i} >> summary.data;
cat ${i}_RM_OUT/${i}.fa.tbl | grep -e ^"base" -e ^"SINE" -e ^"LINE" -e ^"LTR" -e ^"DNA" -e ^"Unclassified" -e ^"Total" -e ^"Simple" -e ^"Low" -e ^"Small" -e ^"Satellites" >> summary.data;
done
```
