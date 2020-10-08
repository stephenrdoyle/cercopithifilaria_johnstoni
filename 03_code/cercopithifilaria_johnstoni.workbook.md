# Workspace for the Cercopithifilaria johnstoni genome project

## Authors: 
Kirsty McCann and Stephen Doyle


# Genome Scope
```bash
# trim the reads first - the full length reads are quite error prone

module load trimmomatic/0.39--1

run_trimmomatic CJ Cj3-500-700_S1_L001_R1_001.fastq.gz Cj3-500-700_S1_L001_R2_001.fastq.gz

# trimmomatic settings: ILLUMINACLIP:/nfs/users/nfs_s/sd21/databases/trimmomatic_Illumina-adapters.fa:2:30:10 CROP:150 SLIDINGWINDOW:10:20 MINLEN:100

# un gizip the raw data
for i in *gz; do
zcat ${i} > ${i%.gz}; done

# run jellyfish
jellyfish count -C -m 21 -s 1000000000 -t 10 *.fastq -o reads.jf

jellyfish histo -t 10 reads.jf > reads.histo

# use the "reads.histo" as input to genomescope

# output of run: http://genomescope.org/analysis.php?code=8JHWvmV1sukF1nGshxjY


```

# genome improvement

- output of spades (as presented in Kirsty's thesis)
/nfs/users/nfs_s/sd21/lustre118_link/cercopithifilaria_johnstoni/scaffolds.fasta

- strategy
     - redundans (without reads) - remove haplotypes
     - blobtools - remove contaminants
     - opera - rescaffold
     - redundans (with reads) - scaffold and gapfill



## redundans
```
/nfs/users/nfs_s/sd21/lustre118_link/software/GENOME_IMPROVEMENT/redundans/redundans.py --fasta scaffolds.fasta --outdir REDUNDANS_OUT

```

## blobtools
```bash

# get taxdump
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -


# blast genome

#--- fix reference, and split it into 100 reads per file - this will help speed up the blast
fastaq to_fasta -l0 scaffolds.reduced.fa scaffolds.reduced.l0.fa

split --lines=200 scaffolds.reduced.l0.fa

#--- run blast, with a separate blast job per 100 reads
for i in x*; do \
bsub.py 10 --queue long --threads 16 blast /lustre/scratch118/infgen/team133/ea10/miniconda/bin/blastn \
-db /data/blastdb/Supported/NT/nt \
-query ${i} \
-outfmt \"6 qseqid staxids bitscore std\" \
-max_target_seqs 10 \
-max_hsps 1 \
-evalue 1e-25 \
-num_threads 16 \
-out ${i}_blast.out;
done

cat *blast.out > blast.out
rm x*


# generate coverage data
#--- use minimap to map reads to generate a bam for blobtools coverage stats

minimap2 -ax sr -t 16 scaffolds.reduced.fa Cj3-500-700_S1_L001_R1_001.fastq.gz Cj3-500-700_S1_L001_R2_001.fastq.gz | samtools sort -@16 -O BAM -o cj_mapped.bam -



# run blobtools
/nfs/users/nfs_s/sd21/lustre118_link/software/ASSEMBLY_QC/blobtools/blobtools create -i scaffolds.reduced.fa -b cj_mapped.bam -t blast.out -o CJ_blobtools_out

/nfs/users/nfs_s/sd21/lustre118_link/software/ASSEMBLY_QC/blobtools/blobtools plot -i CJ_blobtools_out.blobDB.json

/nfs/users/nfs_s/sd21/lustre118_link/software/ASSEMBLY_QC/blobtools/blobtools view -i CJ_blobtools_out.blobDB.json


# filter blobtools output
#--- only keep blast hits to nematoda and no-hit
grep "Nematoda\|no-hit" CJ_blobtools_out.blobDB.table.txt | cut -f1 > blobtools_filter_byhit.list

#--- only keep sequences with the read depth >10
awk '{if($3>10) print $1}' CJ_blobtools_out.cj_mapped.bam.cov > blobtools_filter_bycov_gt10.list

#--- only keep sequences that pass both filters
cat blobtools_filter_byhit.list blobtools_filter_bycov_gt10.list | sort | uniq -c | awk '{if($1==2) print $2}' > blobtools_filter_bycov_gt10_byhit.list


# remake fasta using filtered list of sequences
while read i; do samtools faidx scaffolds.reduced.fa ${i} >> scaffolds.reduced.bt_filter.fa; done < blobtools_filter_bycov_gt10_byhit.list

```



## rescaffold with opera
```bash
# run bowtie to generate a mapping file
bowtie2-build scaffolds.reduced.bt_filter.fa contigs

bowtie2 -k 5 -x contigs -X 1000 --rf -p 32 -1 Cj3-500-700_S1_L001_R1_001.fastq.gz -2 Cj3-500-700_S1_L001_R2_001.fastq.gz -S cj_bowtie2.sam

# run opera
/nfs/users/nfs_s/sd21/lustre118_link/software/GENOME_SCAFFOLDING/OPERA-LG_v2.0.4/bin/OPERA-LG scaffolds.reduced.bt_filter.fa cj_bowtie2.sam opera_out
```



## rescaffold the opera output
```
# run redundans, this time with reads, to scaffold and gapfill
/nfs/users/nfs_s/sd21/lustre118_link/software/GENOME_IMPROVEMENT/redundans/redundans.py -f scaffoldSeq.fasta -i ../Cj3-500-700_S1_L001_R1_001.fastq.gz ../Cj3-500-700_S1_L001_R2_001.fastq.gz -o rescaffold -t 16


```




## check genome completeness using BUSCO and CEGMA
```bash
cd /nfs/users/nfs_s/sd21/lustre118_link/cercopithifilaria_johnstoni/BUSCO

bsub.py --queue long --threads 30 20 busco5 ~sd21/bash_scripts/run_busco_nematode.sh ORIGINAL original.fa
bsub.py --queue long --threads 30 20 busco6 ~sd21/bash_scripts/run_busco_nematode.sh REDUNDANS redundans.fa
bsub.py --queue long --threads 30 20 busco7 ~sd21/bash_scripts/run_busco_nematode.sh BLOBTOOLS blobtools_filter.fa
bsub.py --queue long --threads 30 20 busco10 ~sd21/bash_scripts/run_busco_nematode.sh OPERA opera.fa
bsub.py --queue long --threads 30 20 busco11 ~sd21/bash_scripts/run_busco_nematode.sh OPERA_RESCAFFOLD opera_rescaffold.fa

```



## AUGUSTUS
```bash
/software/pathogen/external/apps/usr/bin/splitMfasta.pl HAEM_V3.3.chr.fa --outputpath=./
for f in *.split.*; do NAME=`grep ">" $f`; mv $f ${NAME#>}.fa; done


/software/pathogen/external/apps/usr/bin/summarizeACGTcontent.pl HAEM_V3.3.chr.fa > basesummary.out




--species=/nfs/users/nfs_s/sd21/software/augustus-3.2.1/config/species/BUSCO_OPERA_RESCAFFOLD_busco3.02.nematoda_3730906830

grep "bases" basesummary.out | awk -v PWD=$PWD -v HINTS=merged_hints.gff '{print PWD"/"$3".fa",PWD"/"HINTS,"1",$1}' OFS="\t" > sequences.list

/software/pathogen/external/apps/usr/bin/createAugustusJoblist.pl --sequences=sequences.list --wrap="#" --overlap=50000 --chunksize=3000000 --outputdir=augustus_split_out --joblist=jobs.lst --jobprefix=augsplit --command \
"/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/augustus_v3.3/bin/augustus --species=HAEM_V3.3.chr2 --strand=both --genemodel=partial --protein=on --introns=on --start=on --stop=on --cds=on --codingseq=on --UTR=off --nc=off --gff3=on --alternatives-from-evidence=on --extrinsicCfgFile=/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/augustus_v3.3/config/species/HAEM_V3.3.chr2/extrinsic.HAEM_V3.3.chr2.cfg --AUGUSTUS_CONFIG_PATH=/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/augustus_v3.3/config/"



mkdir augustus_split_out


for i in augsplit*; do echo -e "bsub.py 4 augsplit_log ./${i}" >> run_augsplit; done        #Max memory used was 1.2Gb
chmod a+x run_augsplit
bsub.py --queue yesterday 1 run_splitter ./run_augsplit



cat augustus_split_out/*gff | /software/pathogen/external/apps/usr/bin/join_aug_pred.pl > augustus_merge.final.gff
grep "AUGUSTUS" augustus_merge.final.gff > augustus_merge.final_2.gff


#update to augustus
/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/augustus_v3.3/bin/augustus
--AUGUSTUS_CONFIG_PATH=/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/augustus_v3.3/config
--extrinsicCfgFile=/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/augustus_v3.3/config/species/HAEM_V3.3.chr2/extrinsic.HAEM_V3.3.chr2.cfg

```


#--------------------------------------------------------------------------------------------------------------

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


#---------------------------------------------------------------------------------------------------------------------



# Assembly statistics for filarial nematodes + C. johnstoni
```
cd miniconda3/bin/

#---New C. johnstoni genome assembly used
#---nematode genome versions recorded

assembly-stats /home/kirstmac/Documents/Files/cercopithifilaria_johnstoni-master/02_data/cjohnstoni_genome_200917.fasta
assembly-stats /home/kirstmac/Documents/whole-genomes/acanthocheilonema_viteae.PRJEB1697.WBPS13.genomic.fa
assembly-stats /home/kirstmac/Documents/whole-genomes/brugia_malayi.PRJNA10729.WBPS13.genomic.fa
assembly-stats /home/kirstmac/Documents/whole-genomes/brugia_pahangi.PRJEB497.WBPS13.genomic.fa
assembly-stats /home/kirstmac/Documents/whole-genomes/brugia_timori.PRJEB4663.WBPS13.genomic.fa
assembly-stats /home/kirstmac/Documents/whole-genomes/dirofilaria_immitis.PRJEB1797.WBPS13.genomic.fa
assembly-stats /home/kirstmac/Documents/whole-genomes/litomosoides_sigmodontis.PRJEB3075.WBPS13.genomic.fa
assembly-stats /home/kirstmac/Documents/whole-genomes/loa_loa.PRJNA246086.WBPS13.genomic.fa
assembly-stats /home/kirstmac/Documents/whole-genomes/onchocerca_flexuosa.PRJEB512.WBPS13.genomic.fa
assembly-stats /home/kirstmac/Documents/whole-genomes/onchocerca_ochengi.PRJEB1204.WBPS13.genomic.fa
assembly-stats /home/kirstmac/Documents/whole-genomes/onchocerca_volvulus.PRJEB513.WBPS14.genomic.fa
assembly-stats /home/kirstmac/Documents/whole-genomes/wuchereria_bancrofti.PRJEB536.WBPS13.genomic.fa

```


#---------------------------------------------------------------------------------------------------------------------

# Filarial nematode BUSCO analyses whole genomes

```
cd /Documents/Programs/busco-master

python scripts/run_BUSCO.py -i /home/kirstmac/Documents/whole-genomes/brugia_pahangi.PRJEB497.WBPS13.genomic.fa -o B_pahangi_Busco -l nematoda_odb9/ -m geno
python scripts/run_BUSCO.py -i /home/kirstmac/Documents/whole-genomes/brugia_timori.PRJEB4663.WBPS13.genomic.fa -o B_timori_Busco -l nematoda_odb9/ -m geno
python scripts/run_BUSCO.py -i /home/kirstmac/Documents/whole-genomes/acanthocheilonema_viteae.PRJEB1697.WBPS13.genomic.fa -o A_viteae_BUSCO -l nematoda_odb9/ -m geno
python scripts/run_BUSCO.py -i /home/kirstmac/Documents/whole-genomes/brugia_malayi.PRJNA10729.WBPS13.genomic.fa -o B_malayi_BUSCO -l nematoda_odb9/ -m geno
python scripts/run_BUSCO.py -i /home/kirstmac/Documents/whole-genomes/dirofilaria_immitis.PRJEB1797.WBPS13.genomic.fa -o D_immitis_BUSCO -l nematoda_odb9/ -m geno
python scripts/run_BUSCO.py -i /home/kirstmac/Documents/whole-genomes/litomosoides_sigmodontis.PRJEB3075.WBPS13.genomic.fa -o L_sigmodontis_BUSCO -l nematoda_odb9/ -m geno
python scripts/run_BUSCO.py -i /home/kirstmac/Documents/whole-genomes/loa_loa.PRJNA246086.WBPS13.genomic.fa -o L_loa_BUSCO -l nematoda_odb9/ -m geno
python scripts/run_BUSCO.py -i /home/kirstmac/Documents/whole-genomes/onchocerca_flexuosa.PRJEB512.WBPS13.genomic.fa -o O_flexuosa_BUSCO -l nematoda_odb9/ -m geno
python scripts/run_BUSCO.py -i /home/kirstmac/Documents/whole-genomes/onchocerca_ochengi.PRJEB1204.WBPS13.genomic.fa -o O_ochengi_BUSCO -l nematoda_odb9/ -m geno
python scripts/run_BUSCO.py -i /home/kirstmac/Documents/whole-genomes/onchocerca_volvulus.PRJEB513.WBPS14.genomic.fa -o O_volvulus_BUSCO -l nematoda_odb9/ -m geno
python scripts/run_BUSCO.py -i /home/kirstmac/Documents/whole-genomes/wuchereria_bancrofti.PRJEB536.WBPS13.genomic.fa -o W_bancrofti_BUSCO -l nematoda_odb9/ -m geno


#---just noticed I have used a different version of Wb, Ls and Of genomes for all my analyses than what you have Steve. Will amend these and adjust notes when I have.

```


#---------------------------------------------------------------------------------------------------------------------



# BUSCO for new CJ proteome and filarial nematode proteomes
```

cd /Documents/Programs/busco-master

#---using the new CJ genome annotation that was trained using BUSCO

python scripts/run_BUSCO.py -i /home/kirstmac/Documents/Files/cercopithifilaria_johnstoni-master/200920_CJ_BUSCO.fa -o CJ_BUSCOtrained-assembly -l nematoda_odb9/ -m prot

#---then running the same analysis on the remaining filarial nematodes to compare CJ against
#---note that all the nematode protein files are from the same genome versions used from assembly statistics and all other analyses

python scripts/run_BUSCO.py -i /home/kirstmac/Documents/Files/Proteins/A_viteae.fa -o Av-proteins -l nematoda_odb9/ -m prot
python scripts/run_BUSCO.py -i /home/kirstmac/Documents/Files/Proteins/B_malayi.fa -o Bm-proteins -l nematoda_odb9/ -m prot
python scripts/run_BUSCO.py -i /home/kirstmac/Documents/Files/Proteins/B_pahangi.fa -o Bp-proteins -l nematoda_odb9/ -m prot
python scripts/run_BUSCO.py -i /home/kirstmac/Documents/Files/Proteins/B_timori.fa -o Bt-proteins -l nematoda_odb9/ -m prot
python scripts/run_BUSCO.py -i /home/kirstmac/Documents/Files/Proteins/D_immitis.fa -o Di-proteins -l nematoda_odb9/ -m prot
python scripts/run_BUSCO.py -i /home/kirstmac/Documents/Files/Proteins/L_loa.fa -o Ll-proteins -l nematoda_odb9/ -m prot
python scripts/run_BUSCO.py -i /home/kirstmac/Documents/Files/Proteins/L_sigmodontis.fa -o Ls-proteins -l nematoda_odb9/ -m prot
python scripts/run_BUSCO.py -i /home/kirstmac/Documents/Files/Proteins/O_flexuosa.fasta -o Of-proteins -l nematoda_odb9/ -m prot
python scripts/run_BUSCO.py -i /home/kirstmac/Documents/Files/Proteins/O_ochengi.fasta -o Oc-proteins -l nematoda_odb9/ -m prot
python scripts/run_BUSCO.py -i /home/kirstmac/Documents/Files/Proteins/O_volvulus.fasta -o Ov-proteins -l nematoda_odb9/ -m prot
python scripts/run_BUSCO.py -i /home/kirstmac/Documents/Files/Proteins/W_bancrofti.fa -o Wb-proteins -l nematoda_odb9/ -m prot

```


#---------------------------------------------------------------------------------------------------------------------

# 12S-COI phylogeny
```

module load raxml-gcc/8.0.19

echo "Starting at $(date)"
raxmlHPC -s 20-09-30_12S-COI_alignment_reducedoutgroups_gb.phy -n 20-09-30_12SCOI_reducedoutgroups_best -m GTRCAT -p 6 -T 2 -# 20
raxmlHPC -s 20-09-30_12S-COI_alignment_reducedoutgroups_gb.phy -n 20-09-30_12SCOI_reducedoutgroups.bootall -m GTRCAT -p 6 -b 6 -T 2 -# 1000

#---combine bootstrap with best 12S-COI tree

>raxmlHPC.exe -f b -t D:\PhDanalyses2\C_johnstoni_analyses\Cercopithifilaria_genes\201001_reduced_gblocks_tree\RAxML_result.20-09-30_12SCOI_reducedoutgroups_best.RUN.6 -z D:\PhDanalyses2\C_johnstoni_analyses\Cercopithifilaria_genes\201001_reduced_gblocks_tree\RAxML_bootstrap.20-09-30_12SCOI_reducedoutgroups.bootall -m GTRCAT -n CJ_12S-COI_reduced.bipart

```


#---------------------------------------------------------------------------------------------------------------------

# Circos plot
```
#---Running analysis on the linux system (ASUS) - the HPC not set up properly
#---run PROMER

cd Documents/Programs/MUMmer3.23/
./promer --mum --prefix CJv2_OV O_volvulus_chromosomes.fa /home/kirstmac/Documents/Files/cercopithifilaria_johnstoni-master/02_data/cjohnstoni_genome_200917.fasta
	
#---run SHOW COORDS
./show-coords -dTlro -L 100 CJv2_OV.delta > CJv2_OV.coords

#---use NUCMER.2.CIRCOS to create the circos files
perl Nucmer.2.circos.pl --promer --ref_order=relw --parse_contig_names --coord-file /home/kirstmac/Documents/Files/cercopithifilaria_johnstoni-master/Updated-promer-CJv2_OV/CJv2_OV.coords --label_size=20 --min_chr_len=10000 --min_hit_len=1000 --ribbons --no_query_labels --min_ID=80 --query_order=rel --min_query_chr_len=50000 --flipquery

#---amend the housekeeping.conf for analysis to run
#---files that need to be in the working directory:
	circos.conf
		#---edited this by removing the etc/ from etc/housekeeping.conf
	karyotype.txt
	links.txt
	housekeeping.conf
		#---editing the
			max_ticks            = 5000
			max_ideograms        = 550
			max_links            = 200000
			max_points_per_track = 300000
	
module load circos/0.67_5
circos

```
#---------------------------------------------------------------------------------------------------------------------

# GC vs Coverage plot
```bash

#---mapping reads to the updated CJ genome
#---index the genome
bwa index /home/kirstmac/Documents/Files/cercopithifilaria_johnstoni-master/02_data/cjohnstoni_genome_200917.fasta
bwa mem /home/kirstmac/Documents/Files/cercopithifilaria_johnstoni-master/02_data/cjohnstoni_genome_200917.fasta /home/kirstmac/Documents/Files/cj_paired_R1.fastq /home/kirstmac/Documents/Files/cj_paired_R2.fastq > CJraw_mapped_201004.sam

#---sort the file for use in bedtools
samtools sort CJraw_mapped_201004.sam > CJraw_sorted_201005.bam

#---index the genome
samtools faidx /home/kirstmac/Documents/Files/cercopithifilaria_johnstoni-master/02_data/cjohnstoni_genome_200917.fasta

#---move to new code for full scaffolds if required
---------------------------------------------------------------------------------------------------------------------------------------------------------
#---cut the first two columns from the index file - contig name and length - for bedtools to make windows
cut -f1,2 cjohnstoni_genome_200917.fasta.fai > CJ.genome

#---make windows - choosing 10000
bedtools makewindows -g CJ.genome -w 10000 > CJ_genome_10k.bed

#---Run bedtools coverage to determine the coverage of CJ reads
#---SUPER IMPORTANT - version bedtools-gcc/2.20.1 - works when running the following command. Version 2.26.0 does not work. Most likely because the commands have switched around but not clear why they don't produce the same output.
    
bedtools coverage -b CJ_genome_10K.bed -abam CJraw_sorted_201005.bam > coverage_10k_201007_V3.txt

#---GC content using bedtools nuc tools

bedtools nuc -fi /home/kirstmac/Documents/Files/cercopithifilaria_johnstoni-master/02_data/cjohnstoni_genome_200917.fasta -bed /home/kirstmac/Documents/Files/cercopithifilaria_johnstoni-master/CJ_genome_10K.bed > GCcontent_201007.txt

#---Combine the GC analysis with the coverage analysis
Paste GCcontent_201007.txt coverage_10K_201007_V3.txt > GC_COV_201007.txt
---------------------------------------------------------------------------------------------------------------------------------------------------------

#---windows of 10000 did not work needed to have the entire scaffolds
#---new code for full scaffolds

samtools dict cjohnstoni_genome_200917.fasta | cut -f2,3 | awk -F'[\t:]' '$1=="SN" {print $2,"1",$4}' OFS="\t" > CJgenome_new.bed

bedtools coverage -b CJgenome_new.bed -abam CJraw_sorted_201005.bam > coverage_new_201007.txt

bedtools nuc -fi cjohnstoni_genome_200917.fasta -bed CJgenome_new.bed > GCcontent_201007_new.txt


# extract promer hits linking Cj and Ov sequences
show-coords -lTH CJv2_OV.delta | awk '{print $15,$14}' OFS="\t" | sort | uniq -c | sort -k2,2 -k1,1nr | awk '!seen[$2]++ {print $2,$3}' OFS="\t" > cj_ov_promer.hits

# make a new file to output to
>cj_ov_promer.complete.hits

# make of list of Cj scaffolds in the reference, check to see if they are in the promer output - if yes, print the Cj-Ov promer hit, if no, print the Cj name and that there was no hit
samtools dict scaffolds.reduced.fa |\
	cut -f2,3 |\
	awk -F'[\t:]' '$1=="SN" {print $2}' OFS="\t" |\
	sort | uniq |\
	while read -r NAME; do if grep -q ${NAME} cj_ov_promer.hits; then grep ${NAME} cj_ov_promer.hits >> cj_ov_promer.complete.hits; else echo -e ${NAME}\\t"no_hit" >> cj_ov_promer.complete.hits; fi; done
	
	
paste GCcontent_201007_new.txt coverage_new_201007.txt cj_ov_promer.complete.hits > NEW_GC_COV_201007.txt
```

```R
#---R script

library(ggplot2)

A <- read.csv(file="D://PhDanalyses2//C_johnstoni_analyses//whole_genome//20-10-07_updated_GCcov//NEW_GC_COV_201007_OVchromdata_MATCHED.csv")

ggplot(A, aes(x=GC, y=COV/150)) + geom_point() + ylim(0,80)

SS <- ggplot() + geom_point(data = A,
                            aes(x = GC, 
                                y = COV/150,
                                color = chromosome), size = 1.5) + 
  ylim(0,1500) +
  xlim(0.5,0.9) +
  labs(y = "Coverage (read depth)", x = "GC content (%)") + 
  labs(color = "O_volvulus chromosomes")

```
SS 
![Chromosome Figure](../04_analysis/GCcov_plot_20-10-08.png)

```
#---plot generated with the log scale, log10

LOG <- ggplot() + geom_point(data = A,
                            aes(x = GC, 
                                y = log10(COV/150),
                                color = chromosome), size = 1.5) + 
  ylim(-2,4) +
  xlim(0.5,0.95) +
  labs(y = "Coverage log10 (read depth)", x = "GC content (%)") + 
  labs(color = "O_volvulus chromosomes")

```
LOG
![Chromosome Figure](../04_analysis/GCcov_plot_log10_201008.png)

```
#---seperate chromosome plots

CHROM <- ggplot() + geom_point(data = A,
                             aes(x = GC, 
                                 y = log10(COV/150),
                                 color = chromosome), size = 1.5) + 
  facet_grid(chromosome ~.) +
  ylim(-2,4) +
  xlim(0.5,0.95) +
  labs(y = "Coverage log10 (read depth)", x = "GC content (%)") + 
  labs(color = "O_volvulus chromosomes")

```
CHROM
![Chromosome Figure](../04_analysis/seperate_plot_chrom_20-10-08.png)

