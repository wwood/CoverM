library(data.table)

contigs = fread('cargo run --release -- contig -b v80.public.B06.V06.A.pos.lim.R1stq_filtered.bam',header=F)
genomes = fread('cargo run --release -- genome -s . -b v80.public.B06.V06.A.pos.lim.R1stq_filtered.bam --min-covered-fraction 0',header=F)
genomes_coverage_hist = fread('cargo run --release -- genome -s . -b v80.public.B06.V06.A.pos.lim.R1stq_filtered.bam --min-covered-fraction 0 -m coverage_histogram', header=F)

setnames(contigs, c('sample','contig','coverage'))
setnames(genomes, c('sample','genome','mean_coverage'))
setnames(genomes_coverage_hist, c('sample','genome','coverage','count'))

stopifnot(
  nrow(merge(genomes_coverage_hist[, .(calculated_mean_coverage=sum(coverage*count)/sum(count)), by='genome'], genomes[mean_coverage > 0], all=T, by='genome')[abs(calculated_mean_coverage - mean_coverage)>0.0001]) == 0)

## Are the correct lengths of the genomes found?
lengths = fread('samtools view -H v80.public.B06.V06.A.pos.lim.R1stq_filtered.bam |grep -v \'^@PG\'',header=F)
l2 = lengths[2:nrow(lengths), .(contig = gsub('SN:','',V2), contig_length=as.numeric(gsub('LN:','',V3)))]
l2[, genome := gsub('\\..*','',contig)]
l2[, .(genome_length=sum(contig_length)), by=genome]
stopifnot(nrow(merge(l2[, .(genome_length=sum(contig_length)), by=genome], genomes_coverage_hist[,.(coverage_hist_len=sum(count)), by=genome], by='genome',all=T)[!is.na(coverage_hist_len)][genome_length != coverage_hist_len]) == 0)

## Are there the correct number of contigs in the contig output?
stopifnot(nrow(contigs) == nrow(l2))
