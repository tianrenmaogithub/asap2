#! python

import os
import argparse
home=os.path.dirname(os.path.realpath(__file__))
exec(open(home+'/../lib/functions.py','r').read())

parser=argparse.ArgumentParser(prog='ASAP', description='This program wraps all the process of QIIME 2 including  importing data, demultiplexing, denoise, dereplicate, ASV, classification, ploting, diversity analysis (alpha, beta) etc. It accepts many data type, including multiplexed and demultiplexed, single-end and paired-end, barcode inside and barcode outside, read data and feature data. For a better view of log, please use \'python -u\' to run this script.') 
parser.add_argument('-v','--version', action='version', version='%(prog)s 2.2')
parser.add_argument('-i','--inputdir',dest='indir',required=True,help='Directory of the input data. See the manual for preparation of input data.')
parser.add_argument('-c',dest='clas',required=True,help='Classifier model for taxonomic classification. Download classifier models from https://docs.qiime2.org/2020.11/data-resources/ and put them in '+home+'/../lib/classifier/ if you haven\'t done it.',choices=os.listdir(home+'/../lib/classifier/'))
parser.add_argument('-q','--qual',dest='qual',required=False,default=20,type=int,help='Quality score cutoff to trim the demultiplexed sequences, 10-40 [default: 20].')
parser.add_argument('-g','--gap',dest='gap',required=False,default=0.2,type=float,help='Maximum percent of gap in a column in the alignment filtering for phylogenetic tree construction, 0-1 [default: 0.2].')
parser.add_argument('-r','--resampling-depth',dest='resa',required=False,default=0,type=int,help='Resampling depth. Make it 0 if you want the pipeline to decide it [default: 0].')
parser.add_argument('-s','--step',dest='step',required=False,default=30,type=int,help='Step number for alpha rarefaction curve, 10-100 [default: 30].')
parser.add_argument('-t','--iteration',dest='iter',required=False,default=30,type=int,help='Iteration number for alpha rarefaction curve, 10-100 [default: 30].')
parser.add_argument('-n','--confidence',dest='conf',required=False,default=0.7,type=float,help='Confidence cutoff of taxonomic classification, 0.5-1 [default: 0.7].')
parser.add_argument('-p','--cpu',dest='proc',required=False,default=10,type=int,help='Processors to use in certain steps such as dada2 and taxonomic classification, 1-10 [default: 10].')
parser.add_argument('-o','--outdir',dest='outdir',required=False,default='asap2_out',help='Directory for outputs [Default: asap2_out]')
args=parser.parse_args()


# check parameters
if args.qual>40 or args.qual<10:
    raise ValueError('-q Quality score should be 10 - 40')
if args.gap>1 or args.gap<0:
    raise ValueError('-g Gap percentage should be 0 - 1')
if args.step>100 or args.step<10:
    raise ValueError('-s Step should be 10 - 100')
if args.iter>100 or args.iter<10:
    raise ValueError('-t Iteration should be 10 - 100')
if args.proc>10 or args.proc<1:
    raise ValueError('-p CPU should be 1 - 10')
if args.conf>1 or args.conf<0.5:
    raise ValueError('-n Confidence should be 0.5 - 1')

# check input data format
print('\n---Checking the input directory for data format---\n')
checkDir(args.indir)
call('rm -rf '+args.outdir+'; mkdir '+args.outdir,args.outdir+' created.','Error: fail to create '+args.outdir+'/.')
print('\nOutput will be put in '+args.outdir+'\n')


# check metadata format


# import data
print('\n---Importing data from '+args.indir+'/---\n')
call('mkdir '+args.outdir+'/imported',args.outdir+'/imported/ created.','Error: fail to create '+args.outdir+'/imported/.')
importDir(args.indir,args.outdir+'/imported')


# demultiplex
print('\n---Demultiplexing data---\n')
call('mkdir '+args.outdir+'/demultiplexed',args.outdir+'/demultiplexed/ created.','Error: fail to create '+args.outdir+'/demultiplexed/.')

for i in os.listdir(args.outdir+'/imported'):
    project=re.sub(r'\.qza$','',i)
    if i[:4]=='fqMu':
        demultiplex(args.outdir+'/imported/'+i,args.indir+'/'+project+'/metadata.tsv',args.outdir+'/demultiplexed')
    elif i[:4]=='fqDe':
        call('cp '+args.outdir+'/imported/'+i+' '+args.outdir+'/demultiplexed/'+project+'-demux.qza','Already demultiplexed data: \n'+args.outdir+'/imported/'+i+' copied to '+args.outdir+'/demultiplexed/'+project+'-demux.qza')


# summarize and visualize demultiplexed qza
print('\n---Summarising demultiplexed qza---\n')
call('mkdir '+args.outdir+'/demuxSumViz',args.outdir+'/demuxSumViz/ created.','Error: fail to create '+args.outdir+'/demuxSumViz/.')
demuxSumViz(args.outdir+'/demultiplexed',args.outdir+'/demuxSumViz')
call('mkdir '+args.outdir+'/readSummary',args.outdir+'/readSummary created.','Error: fail to create '+args.outdir+'/readSummary.')
readProfileDemux(args.outdir+'/demuxSumViz',args.outdir+'/readSummary/read_summary_demux.tsv')


# quality trimming, denoise, chrimeric, abundance
print('\n---Quality trimming, denoise, chrimeric, abundance---\n')
call('mkdir '+args.outdir+'/feature',args.outdir+'/feature/ created.','Error: fail to create '+args.outdir+'/feature/.')
for i in os.listdir(args.outdir+'/demultiplexed'):
    if i[-9:]=='demux.qza':
        denoise(args.outdir+'/demultiplexed/'+i,args.outdir+'/feature',q=args.qual,p=args.proc)
readProfileDenoise(args.outdir+'/feature',args.outdir+'/readSummary/read_summary_denoise.tsv')


# convert input feature tables to qza
print('\n---Converting input feature table data to qza---\n')
convertInFt(args.indir,args.outdir+'/feature')


# merge feature tables and sequences of multiple projects
print('\n---Merging multiple projects---\n')
call('mkdir '+args.outdir+'/featureMerged',args.outdir+'/featureMerged created.','Error: fail to create '+args.outdir+'/featureMerged.')
mergeProjects(args.outdir+'/feature',args.outdir+'/featureMerged')
getFileFromZ(args.outdir+'/featureMerged/merged-table.qzv','sample-frequency-detail.csv',args.outdir+'/readSummary/read_summary_merged.tsv')


# merge metadata
print('\n---Merging multiple metadata---\n')
mergeMetadata(args.indir,args.outdir+'/featureMerged/merged-metadata.tsv')


# phylogenetic tree construction
print('\n---Constructing phylogenetic tree---\n')
call('mkdir '+args.outdir+'/phylogeny',args.outdir+'/phylogeny created.','Error: fail to create '+args.outdir+'/phylogeny.')
call('qiime phylogeny align-to-tree-mafft-fasttree --i-sequences '+args.outdir+'/featureMerged/merged-rep-seqs.qza --o-alignment '+args.outdir+'/phylogeny/aligned-rep-seqs.qza --o-masked-alignment '+args.outdir+'/phylogeny/masked-aligned-rep-seqs.qza --o-tree '+args.outdir+'/phylogeny/unrooted-tree.qza --o-rooted-tree '+args.outdir+'/phylogeny/rooted-tree.qza --p-mask-max-gap-frequency '+str(args.gap)+' --p-n-threads '+str(args.proc),'Phylogenetic tree saved.','Phylogenetic tree construction failed.')

# alpha beta diversity
print('\n---Alpha and beta diversity and statistical tests---\n')
sampleDepth=resampleDepth(args.outdir+'/readSummary/read_summary_merged.tsv',args.resa)
call('qiime diversity core-metrics-phylogenetic --i-phylogeny '+args.outdir+'/phylogeny/rooted-tree.qza --i-table '+args.outdir+'/featureMerged/merged-table.qza --p-sampling-depth '+str(sampleDepth)+' --m-metadata-file '+args.outdir+'/featureMerged/merged-metadata.tsv --output-dir '+args.outdir+'/diversity','Diversity analysis done successfully.','Diversity analysis failed.')
call('mkdir '+args.outdir+'/diversity/alpha '+args.outdir+'/diversity/beta; mv '+args.outdir+'/diversity/*_vector.qza '+args.outdir+'/diversity/alpha; mv  '+args.outdir+'/diversity/*matrix.qza '+args.outdir+'/diversity/*results.qza '+args.outdir+'/diversity/*emperor.qzv '+args.outdir+'/diversity/beta','Alpha and Beta diversity results moved into alpha/ and beta/.','Failed to move Alpha and Beta diversity results in alpha/ and beta/.')
alphaGroupSig(args.outdir+'/diversity/alpha/',args.outdir+'/featureMerged/merged-metadata.tsv')
alphaCorrelation(args.outdir+'/diversity/alpha/',args.outdir+'/featureMerged/merged-metadata.tsv')


# alpha rarefaction
print('\n---Alpha rarefaction---\n')
r=pd.read_csv(args.outdir+'/readSummary/read_summary_merged.tsv',index_col=0)
md=int(np.median(r.iloc[:,0]))

call('qiime diversity alpha-rarefaction --i-table '+args.outdir+'/featureMerged/merged-table.qza --i-phylogeny '+args.outdir+'/phylogeny/rooted-tree.qza --p-max-depth '+str(md)+' --o-visualization '+args.outdir+'/diversity/alpha/alpha-rarefaction.qzv --p-steps '+str(args.step)+' --p-iterations '+str(args.iter),'Alpha rarefaction done successfully.','Alpha rarefaction failed.')
call('qiime tools export --input-path '+args.outdir+'/diversity/rarefied_table.qza --output-path '+args.outdir+'/featureMerged','Rarefied feature table exported to '+args.outdir+'/featureMerged/feature-table.biom','Exporting rarefied feature table failed.')
call('biom convert -i '+args.outdir+'/featureMerged/feature-table.biom -o '+args.outdir+'/featureMerged/feature-table.tsv --to-tsv','Converting to tsv done successfully.','Converting to tsv failed.')
call('mv '+args.outdir+'/featureMerged/feature-table.biom '+args.outdir+'/featureMerged/merged-table-rarefied.biom; mv '+args.outdir+'/featureMerged/feature-table.tsv '+args.outdir+'/featureMerged/merged-table-rarefied.tsv','Renaming to merged-table-rarefied.biom and merged-table-rarefied.tsv done successfully.','Renaming failed.')


# taxonomic classification
print('\n---Taxonomic classification and visualization---\n')
call('mkdir '+args.outdir+'/taxonomy',args.outdir+'/taxonomy created.','Error: fail to create '+args.outdir+'/taxonomy.')
classification(args.outdir+'/featureMerged/merged-rep-seqs.qza',args.outdir+'/featureMerged/merged-table.qza',home+'/../lib/classifier/'+args.clas,args.outdir+'/featureMerged/merged-metadata.tsv',args.outdir+'/taxonomy',args.conf,args.proc) # changed from 'diversity/rarefied_table.qza' to 'featureMerged/merged-table.qza' in version 2.2
getFileFromZ(args.outdir+'/taxonomy/classification.qza','taxonomy.tsv',args.outdir+'/taxonomy/taxonomy.tsv')


# ordination analysis: CCA, RDA and variable selection 
print('\n---Ordination analysis: CCA, RDA and variable selection---\n')
call('mkdir '+args.outdir+'/ordination',args.outdir+'/ordination created.','Error: fail to create '+args.outdir+'/ordination.')
try:
    ordination(args.outdir+'/taxonomy/taxa-bar-plots.qzv',args.outdir+'/featureMerged/merged-metadata.tsv',args.outdir+'/ordination')
except:
    pass

print('\nAll the processes have been finished successfully. Enjoy it!')
