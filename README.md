# bc_umi_gen
Generate random FASTQs with designated barcode and unique molecular identifier structure to test sequencing data processing pipelines.

Single cell sequencing and especially custom library designs heavily utilize whitelisted barcodes (BC), sometimes more than one and split by spacer sequences, and random unique molecular identifiers (UMI), which correct for amplification bias. Here we aim to provide a tool that can flexibly generate synthetic FASTQ files, for testing and presenting small examples.

Currently the tool handles:
1. split barcodes from whitelists
2. UMI of requested length
3. T stretch of requested length
4. mappable sequence from multiple genes in fasta
5. paired-end reads
6. designate general read quality

To include in next version:
1. setting up library design as JSON config
2. collapsible UMI and read pairs
3. custom mismatch rate
4. speed improvements
5. barcodes that need correcting
