#HOMEWORK 5 BIOL7210

#make directory for assignment
mkdir -pv /home/knatarajan37/ex5

# use environment from read cleaning and genome assembly assignment
conda activate ex3

# retrieve files
fasterq-dump \
 SRR3214715 SRR3215024 SRR3215107 \
 --outdir ./raw_data \
 --split-files \
 --skip-technical

cd sra
mamba install -c bioconda fastp

# run fastp for qaqc
# For SRR3214715
fastp \
 -i SRR3214715_1.fastq \
 -I SRR3214715_2.fastq \
 -o ~/ex5/sra/cleaned_SRR3214715_1.fq.gz\
 -O ~/ex5/sra/cleaned_SRR3214175_2.fq.gz \
 --html SRR3214715_fastp_report.html
# For SRR3215024
fastp \
 -i SRR3215024_1.fastq \
 -I SRR3215024_2.fastq \
 -o ~/ex5/sra/cleaned_SRR3215024_1.fq.gz\
 -O ~/ex5/sra/cleaned_SRR3215024_2.fq.gz \
 --html SRR3215024_fastp_report.html
# For SRR3215107
fastp \
 -i SRR3215107_1.fastq \
 -I SRR3215107_2.fastq \
 -o ~/ex5/sra/cleaned_SRR3215107_1.fq.gz\
 -O ~/ex5/sra/cleaned_SRR3215107_2.fq.gz \
 --html SRR3215107_fastp_report.html

# run skesa for assembly

# For SRR3214715
spades.py --pe1-1 ~/ex5/sra/cleaned_SRR3214715_1.fq.gz --pe1-2 ~/ex5/sra/cleaned_SRR3214715_2.fq.gz -o ~/ex5/trim/SRR3214715_assembly_k17 -k 17 --memory 30 --isolate --threads 4

# For SRR3215024
spades.py --pe1-1 ~/ex5/sra/cleaned_SRR3215024_1.fq.gz --pe1-2 ~/ex5/sra/cleaned_SRR3215024_2.fq.gz -o ~/ex5/trim/SRR3215024_assembly_k17 -k 17 --memory 30 --isolate --threads 4

# For SRR3215107
spades.py --pe1-1 ~/ex5/sra/cleaned_SRR3215107_1.fq.gz --pe1-2 ~/ex5/sra/cleaned_SRR3215107_2.fq.gz -o ~/ex5/trim/SRR3215107_assembly_k17 -k 17 --memory 30 --isolate --threads 4



#check all assemblies are approximately the same size
 ls -alh SRR3215107_assembly_k17/contigs.fasta SRR3215024_assembly_k17/contigs.fasta SRR3214715_assembly_k17/contigs.fasta

# run MLST with docker
 mlst ~/ex5/trim/SRR3214715_assembly_k17/contigs.fasta > mlst.tsv
 mlst ~/ex5/trim/SRR3215024_assembly_k17/contigs.fasta > mlst.tsv
 mlst ~/ex5/trim/SRR3214715_assembly_k17/contigs.fasta > mlst2.tsv
 mlst ~/ex5/trim/SRR3215107_assembly_k17/contigs.fasta > mlst3.tsv
 

# concatenate all MLST summary files
cat *mlst.tsv > mlst_noheaders.tsv ### I renamed this in the end

conda deactivate
conda activate fastani

# download heliobactor type strain genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/017/821/535/GCF_017821535.1_ASM1782153v1/GCF_017821535.1_ASM1782153v1_genomic.fna.gz
gunzip *.fna
mv GCF_017821535.1_ASM1782153v1_genomic.fna reference.fna


# run fastani
fastANI \
  --query  ~/ex5/trim/SRR3214715_assembly_k17/contigs.fasta \
  --ref reference.fna \
  --fragLen 1000 \
  --output SRR3214715_fastANI_Output.tsv

fastANI \
  --query ~/ex5/trim/SRR3215024_assembly_k17/contigs.fasta \
  --ref reference.fna \
  --fragLen 1000 \
  --output SRR3215024_fastANI_Output.tsv

  fastANI \
  --query /~/ex5/trim/SRR3215107_assembly_k17/contigs.fasta \
  --ref reference.fna \
  --fragLen 1000 \
  --output SRR3215107_fastANI_Output.tsv

# combine tsv files
cat *_fastANI_Output.tsv > FastANI_Output.tsv

# calculate alignment percent and length
awk \
  '{alignment_percent = $4/$5*100} \
   {alignment_length = $4*3000} \
   {print $0 "\t" alignment_percent "\t" alignment_length}' \
  FastANI_Output.tsv \
  > FastANI_Output_With_Alignment.tsv
sed \
  "1i Query\tReference\t%ANI\tNum_Fragments_Mapped\tTotal_Query_Fragments\t%Query_Aligned\tBasepairs_Query_Aligned" \
  FastANI_Output_With_Alignment.tsv \
  > FastANI_Output_With_Alignment_With_Header.tsv
column -ts $'\t' FastANI_Output_With_Alignment_With_Header.tsv | less -S

# replace query names with accession IDs
awk 'BEGIN {FS=OFS="\t"} NR==1 {print $0} NR>1 {match($1, /SRR[0-9]+/); print substr($1, RSTART, RLENGTH), $2, $3, $4, $5, $6, $7}' FastANI_Output_With_Alignment_With_Header.tsv > FastANI_Output_With_Alignment_With_Header_modified.tsv

conda deactivate
cd /home/amckee/biol7210/homework/hw5/checkm
conda activate checkm

# quality assessment
# checkm

# get checkm database
wget https://zenodo.org/records/7401545/files/checkm_data_2015_01_16.tar.gz
tar zxvf checkm_data_2015_01_16.tar.gz
echo 'export CHECKM_DATA_PATH=/home/knatarajan37/ex5/checkm/db' >> ~/.bashrc
source ~/.bashrc
echo "${CHECKM_DATA_PATH}"

checkm taxon_list | grep Campylobacter
checkm taxon_set species "Campylobacter jejuni" Cj.markers
checkm \
  analyze \
  Cj.markers \
  ~ex5/checkm \
  analyze_output
checkm \
  qa \
  -f checkm.tax.qa.out \
  -o 1 \
  Cj.markers \
  analyze_output
sed 's/ \+ /\t/g' checkm.tax.qa.out > checkm.tax.qa.out.tsv
cut -f 2- checkm.tax.qa.out.tsv > tmp.tab && mv tmp.tab checkm.tax.qa.out.tsv
sed -i '1d; 3d; $d' checkm.tax.qa.out.tsv

# extract the first 
awk 'BEGIN {FS=OFS="\t"} NR==1 {print $0} NR>1 {sub(/_filtered_assembly$/, "", $1); print}' checkm.tax.qa.out.tsv > quality.tsv
