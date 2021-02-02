We developed Amplicon Sequence Analysis Pipeline 2 (ASAP 2) to analyze marker gene amplicon sequence data automatically and consistently. It was designed to be us
er-friendly and time-saving. Users just need to organize their fastq (demultiplexed or multiplexed, single-end or pair-end) or fasta (demultiplexed joined paired-
end or single-end) and metadata properly. Multiple projects can be put together and the pipeline will automatically merge the data. People can even put their own
data and data from collaborators or NCBI together to compare. Also, if you have intermediate files (feature table and representative sequence) from collaborators
and you don't want to repeat the previous steps, just start from them.

######### Dependency #########
Please install the following packages.

Python>=3.6
QIIME2>=2020.6
pandas>=0.25.3
biopython>=1.77
vegan>=2.5-6 # R package

If you install QIIME 2 in a virtual environment, remember to activate it before running ASAP 2.

######### Prepare data files and run #########
Please prepare your sequence file and metadata according to the manual. Check once more the folder structure and naming according to the example data. Then put al
l the project data in one folder (e.g. input/) and run the bin/asap2.py outside the folder.

IMPORTANT: note that only projects using the same region (the same primer set, or regions with >80% overlap if you don't mind the bias resulted from primer amplif
ication) can be put together and merged, or else, the phylogenetic tree will tell a quite different story.

