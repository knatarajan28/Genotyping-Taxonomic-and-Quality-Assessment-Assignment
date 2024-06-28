# Genotyping-Taxonomic-and-Quality-Assessment-Assignment
- Three Helicobacter pylori isolates were sequenced as part of an outbreak analysis. SRA accessions are on NCBI for each: SRR3214715 SRR3215024 SRR3215107
- Use all (3) SRA accessions above and fetch the Illumina sequence data from NCBI. Perform fastANI against the species type strain. You will have to find and fetch the type strain genome.
- Fetch all read sets with sra-tools as performed previously with fasterq-dump
- Quick read clean with fastp or what you're most comfortable with
- Quick assembly with skesa or what you're most comfortable with
- Filter out low coverage and short contigs
- Verify filesizes look similar with ls -alh *.fna in your output directory containing all assemblies. If they're not similar in filesizes, refine trim and assembly parameters. All 3 are highly related and in the outbreak.
- Genotype all 3 assemblies with MLST. For just 1 assembly of the 3, estimate its completeness and contamination levels.
