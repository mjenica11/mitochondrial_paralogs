library(Rsubread)

simReads(# transcript database and wanted abundance
		 transcript.file="/scratch/mjpete11/subread_simulations/fastas/SLC25A4.fa", 
		 expression.levels=10, 
		 # the name of the output
		 output.prefix="1000_reads_SLC25A4", 
		 # options on the output
		 library.size=1000, read.length=100, truth.in.read.names=TRUE, 
		 # simulating sequencing errors
		 simulate.sequencing.error=TRUE, quality.reference=NULL,
		 # parameters for generating paired-end reads
		 paired.end=TRUE, 
		 fragment.length.min=100,
		 fragment.length.max=500,
		 fragment.length.mean=180,
		 fragment.length.sd=40)
