# LazypipePlus
Bergthaler Lab's development version of Lazypipe metagenomic pipeline

+ This version of Lazypipe contains rule for generating coverage plots based on the Virosaurus collection of reference genomes: https://viralzone.expasy.org/8676
+ The rule coverage_plots does bwa-mem alignment to the virosaurus reference fasta and filters out all reads with mapping quality <20. It produces two types of output: 
++ A) it retains multi-mapping reads (files with "allmap"
++ B) removes all multi-mapping reds (files with "uniqmap")
+ Coverage files are produced as static png plots and interactive html files.
