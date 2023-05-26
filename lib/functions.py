#! python

import os
import sys
import re
import subprocess as subp
from collections import Counter
import pandas as pd
import numpy as np
from Bio import Seq
from Bio import SeqIO
import zipfile
import gzip
from io import StringIO
import io

def call(cmd='',out='',err=''):
    if subp.call('set -ex; '+cmd,shell=True)==0:
        if out:
            print('\n'+out+'\n')
    else:
        if err:
            print('\n'+err+'\n')
        sys.exit(1)

# check file presence

def checkFiles(path,files):
    #s1=set([i for i in os.listdir(path) if os.path.isfile(path+'/'+i)])
    s1=set(os.listdir(path))
    s1 = {i for i in s1 if i[0] != '.'}
    s2=set(files)
    if s1==s2:
        print('All files required were found in '+path)
    else:
        if s1-s2:
            print('Error: '+str(s1-s2)+' are not required in '+path+'\nPlease read the manual carefully for data preparation. The test data is a good reference to use.')
        if s2-s1:
            print('Error: '+str(s2-s1)+' are not found in '+path+'\nPlease read the manual carefully for data preparation. The test data is a good reference to use.')
        sys.exit(1)

def checkfqDePe(indir):
    sample=[]
    for i in [f for f in os.listdir(indir) if f[0] != '.']:
        if not re.match(r'.*_.*_.*_R[12]_.*\.fastq.gz',i):
            print('Error: '+indir+'/'+i+' is in wrong format. It should be .*_.*_.*_R[12]_.*\.fastq.gz')
            sys.exit(1)
        sample.append(re.split(r'_R[12]',i)[0])
    bad=[i for i in Counter(sample).keys() if Counter(sample)[i]!=2]
    if bad:
        print('Error: Each sample should have R1 and R2 files. '+str(bad)+' have only one or more than two read files.')
    print('All the formats are correct in '+indir+'.')

def checkfqDeSe(indir):
    for i in [f for f in os.listdir(indir) if f[0] != '.']:
        if not re.match(r'.*_.*_.*_R1_.*\.fastq.gz',i):
            print('Error: '+indir+'/'+i+' is in wrong format. It should be .*_.*_.*_R1_.*\.fastq.gz')
            sys.exit(1)
    print('All the formats are correct in '+indir+'.')

def checkFtDir(indir):
    f=[f for f in os.listdir(indir) if f[0] != '.']
    if not len(f)==3:
        print(str(f)+' found.')
        print('Error: Only three files: feature-table.tsv/biom rep-seqs.fasta and metadata.tsv are required in '+indir)
        sys.exit(1)
    if not ('metadata.tsv' in f and 'rep-seqs.fasta' in f):
        print('Error: metadata.tsv and rep-seqs.fasta are required in '+indir)
        sys.exit(1)
    if not ('feature-table.tsv' in f or 'feature-table.biom' in f):
        print('Error: either feature-table.tsv or feature-table.biom is required in '+indir)
        sys.exit(1)
    print('All the formats are correct in '+indir+'.')

def checkDir(indir):
    for i in [f for f in os.listdir(indir) if f[0] != '.']:
        print('Checking directory '+i)
        if i.split('-')[0]=='fqMuBiPe':
            checkFiles(indir+'/'+i,['seqData','metadata.tsv'])
            checkFiles(indir+'/'+i+'/seqData',['forward.fastq.gz','reverse.fastq.gz'])
        elif i.split('-')[0]=='fqMuBiSe':
            checkFiles(indir+'/'+i,['seqData','metadata.tsv'])
            checkFiles(indir+'/'+i+'/seqData',['sequences.fastq.gz'])
        elif i.split('-')[0]=='fqMuBoPe':
            checkFiles(indir+'/'+i,['seqData','metadata.tsv'])
            checkFiles(indir+'/'+i+'/seqData',['barcodes.fastq.gz','forward.fastq.gz','reverse.fastq.gz'])
        elif i.split('-')[0]=='fqMuBoSe':
            checkFiles(indir+'/'+i,['seqData','metadata.tsv'])
            checkFiles(indir+'/'+i+'/seqData',['barcodes.fastq.gz','sequences.fastq.gz'])
        elif i.split('-')[0]=='fqDePe':
            checkFiles(indir+'/'+i,['seqData','metadata.tsv'])
            checkfqDePe(indir+'/'+i+'/seqData')
        elif i.split('-')[0]=='fqDeSe':
            checkFiles(indir+'/'+i,['seqData','metadata.tsv'])
            checkfqDeSe(indir+'/'+i+'/seqData')
        elif i.split('-')[0]=='featureTable':
            checkFtDir(indir+'/'+i)
        else:
            print('Error: '+indir+'/'+i+" is not in the format of one of [featureTable-project_1','fqMuBiPe-project_emp1','fqMuBoSe-project_emp1','fqDePe-project_cas1','fqMuBiSe-project_emp1','fqDeSe-project_cas1','fqMuBoPe-project_emp1']\nYou should put one or multiple folders (e.g. fqMuBiPe-project_emp1) in another folder (e.g. input) and zip it.")
            sys.exit(1)

# import data of different formats

def importData(indir,outdir):
    dirname=os.path.basename(re.sub(r'/+$','',indir))
    fmt=dirname.split('-')[0]

    if fmt=='fqDePe':
        call('qiime tools import --type SampleData[PairedEndSequencesWithQuality] --input-path '+indir+'/seqData --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path '+outdir+'/'+dirname+'.qza', indir+' imported successfully.','Importing '+indir+' failed.')
    elif fmt=='fqDeSe':
        call('qiime tools import --type SampleData[SequencesWithQuality] --input-path '+indir+'/seqData --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path '+outdir+'/'+dirname+'.qza', indir+' imported successfully.','Importing '+indir+' failed.')
    elif fmt=='fqMuBoPe':
        call('qiime tools import --type EMPPairedEndSequences --input-path '+indir+'/seqData --output-path '+outdir+'/'+dirname+'.qza', indir+' imported successfully.','Importing '+indir+' failed.')
    elif fmt=='fqMuBiPe':
        call('qiime tools import --type MultiplexedPairedEndBarcodeInSequence --input-path '+indir+'/seqData --output-path '+outdir+'/'+dirname+'.qza', indir+' imported successfully.','Importing '+indir+' failed.')
    elif fmt=='fqMuBoSe':
        call('qiime tools import --type EMPSingleEndSequences --input-path '+indir+'/seqData --output-path '+outdir+'/'+dirname+'.qza', indir+' imported successfully.','Importing '+indir+' failed.')
    elif fmt=='fqMuBiSe':
        call('qiime tools import --type MultiplexedSingleEndBarcodeInSequence --input-path '+indir+'/seqData/sequences.fastq.gz --output-path '+outdir+'/'+dirname+'.qza', indir+' imported successfully.','Importing '+indir+' failed.')

def importDir(indir,outdir):
    for i in os.listdir(indir):
        print('Importing data '+i)
        importData(indir+'/'+i,outdir)


# demultiplex fqMuBoSe and fqMuBoPe

def grepBarcode(barcode,barcodeFile):
    num1=subp.check_output('zcat '+barcodeFile+' | head -1000000 | grep '+barcode+' | wc -l',shell=True,universal_newlines=True)
    num1=int(num1.strip())

    barcode=str(Seq.Seq(barcode).reverse_complement())
    num2=subp.check_output('zcat '+barcodeFile+' | head -1000000 | grep '+barcode+'| wc -l',shell=True,universal_newlines=True)
    num2=int(num2.strip())

    ratio=num1/(num2+1)
    
    if ratio>10:
        return('AI')
    elif num2>0 and ratio<0.1:
        return('RC')
    else:
        return('ND')

def checkBarcodeOrientation(metadataFile,barcodeFile):
    print('Checking orientation of 10 barcodes ...')
    df=pd.read_csv(metadataFile,sep='\t',header=0,index_col=0,comment='#')
    barcodes=list(df['barcode-sequence'][1:11])
    l=[grepBarcode(i,barcodeFile) for i in barcodes]
    nd=np.array(barcodes)[np.array(l)=='ND']
    asis=np.array(barcodes)[np.array(l)=='AI']
    revc=np.array(barcodes)[np.array(l)=='RC']

    if len(nd)>0:
        print(str(nd)+' are not found in the barcode file.')

    asisNum=Counter(l)['AI']
    revcNum=Counter(l)['RC']
    ratio=asisNum/(revcNum+1)
    print('Barcodes that are the same in the barcode file and matadata file.\nBarcodes: '+' '.join(asis)+'\n#: '+str(asisNum)+'\n')
    print('Barcodes that are the reverse complement in the barcode file and matadata file.\nBarcodes: '+' '.join(revc)+'\n#: '+str(revcNum)+'\n')

    if ratio>2:
        print('The orientation of barcodes in the barcode file and the metadata file are the same.')
        return('AI')
    elif ratio<0.5:
        print('The orientation of barcodes in the barcode file and the metadata file are reverse complement.')
        return('RC')
    else:
        print('Cannot determine if the barcodes in the barcode file and the metadata file are the same or reverse complement.')
        return('ND')

def getReadLen(inQza):
    z = zipfile.ZipFile(inQza)
    f = [i for i in z.filelist if re.search('fastq.gz$',i.filename.split('/')[-1])][0]
    l = [i[1].decode().strip() for i in list(enumerate(gzip.GzipFile(fileobj=io.BytesIO(z.open(f).read()))))[:100]]
    return(round(np.mean([len(l[i]) for i in range(1, 100, 4)])))

def demultiplex(inQza,metadata,outdir):
    print('\nDemultiplexing '+inQza)
    readLen = round(getReadLen(inQza) * 0.9)

    # fqMuBoSe
    if os.path.basename(inQza).split('-')[0]=='fqMuBoSe':
        call('qiime demux emp-single --i-seqs '+inQza+' --m-barcodes-file '+metadata+' --m-barcodes-column barcode-sequence --o-per-sample-sequences '+outdir+'/'+re.sub(r'\.qza$','',os.path.basename(inQza))+'-demux.qza --o-error-correction-details '+outdir+'/'+re.sub(r'\.qza$','',os.path.basename(inQza))+'-demux-details.qza','Demultiplexing '+inQza+' done successfully.','Demultiplexing '+inQza+' failed.')
    
    # fqMuBoPe
    elif os.path.basename(inQza).split('-')[0]=='fqMuBoPe':
        print('Checking barcode direction in metadata file and barcode file.')
        barcodeFile=os.path.dirname(os.path.realpath(metadata))+'/seqData/barcodes.fastq.gz'
        # check if the barcodes in metadata file and barcode file are the same or reverse completmentary
        orien=checkBarcodeOrientation(metadata,barcodeFile)
        if orien=='ND':
            print('Error: Cannot determine direction of barcodes in barcode file and metadata file. Exiting.')
            sys.exit(1)
        elif orien=='AI':
            call('qiime demux emp-paired --i-seqs '+inQza+' --m-barcodes-file '+metadata+' --m-barcodes-column barcode-sequence --o-per-sample-sequences '+outdir+'/'+re.sub(r'\.qza$','',os.path.basename(inQza))+'-demux.qza --o-error-correction-details '+outdir+'/'+re.sub(r'\.qza$','',os.path.basename(inQza))+'-demux-details.qza','Demultiplexing '+inQza+' done successfully.','Demultiplexing '+inQza+' failed.')
        elif orien=='RC':
            call('qiime demux emp-paired --i-seqs '+inQza+' --m-barcodes-file '+metadata+' --m-barcodes-column barcode-sequence --p-rev-comp-mapping-barcodes --o-per-sample-sequences '+outdir+'/'+re.sub(r'\.qza$','',os.path.basename(inQza))+'-demux.qza --o-error-correction-details '+outdir+'/'+re.sub(r'\.qza$','',os.path.basename(inQza))+'-demux-details.qza','Demultiplexing '+inQza+' done successfully.','Demultiplexing '+inQza+' failed.')

    # fqMuBiSe
    elif os.path.basename(inQza).split('-')[0]=='fqMuBiSe':
        call('qiime cutadapt demux-single --i-seqs '+inQza+' --m-barcodes-file '+metadata+' --m-barcodes-column barcode-sequence --o-per-sample-sequences '+outdir+'/'+re.sub(r'\.qza$','',os.path.basename(inQza))+'-demux.qza --o-untrimmed-sequences '+outdir+'/'+re.sub(r'\.qza$','',os.path.basename(inQza))+'-untrimmed.qza --p-minimum-length '+str(readLen),'Cutadapt and demultiplexing '+inQza+' done successfully.','Cutadapt and demultiplexing '+inQza+' failed.')

    # fqMuBiPe
    elif os.path.basename(inQza).split('-')[0]=='fqMuBiPe':
        #check if metadata has reverse read barcoding: barcode-sequence2
        df=pd.read_csv(metadata,sep='\t',header=0,index_col=0,comment='#')
        if 'barcode-sequence2' in df.columns:
            call('qiime cutadapt demux-paired --i-seqs '+inQza+' --m-forward-barcodes-file '+metadata+' --m-forward-barcodes-column barcode-sequence --m-reverse-barcodes-file '+metadata+' --m-reverse-barcodes-column barcode-sequence2 --o-per-sample-sequences '+outdir+'/'+re.sub(r'\.qza$','',os.path.basename(inQza))+'-demux.qza --o-untrimmed-sequences '+outdir+'/'+re.sub(r'\.qza$','',os.path.basename(inQza))+'-untrimmed.qza --p-minimum-length '+str(readLen),'Cutadapt and demultiplexing '+inQza+' done successfully.','Cutadapt and demultiplexing '+inQza+' failed.')
        else:
            call('qiime cutadapt demux-paired --i-seqs '+inQza+' --m-forward-barcodes-file '+metadata+' --m-forward-barcodes-column barcode-sequence --o-per-sample-sequences '+outdir+'/'+re.sub(r'\.qza$','',os.path.basename(inQza))+'-demux.qza --o-untrimmed-sequences '+outdir+'/'+re.sub(r'\.qza$','',os.path.basename(inQza))+'-untrimmed.qza --p-minimum-length '+str(readLen),'Cutadapt and demultiplexing '+inQza+' done successfully.','Cutadapt and demultiplexing '+inQza+' failed.')

def demuxSumViz(indir,outdir):
    for i in os.listdir(indir):
        if i[-10:]=='-demux.qza':
            print('Summarising '+i)
            call('qiime demux summarize --i-data '+indir+'/'+i+' --o-visualization '+outdir+'/'+i[:-10]+'.qzv','Summarising '+indir+'/'+i+' done successfully.','Summarising '+indir+'/'+i+' failed.')

def getRegion(l):
    l1=list(l)
    l1.append(False)
    l2=[]
    n=0
    for i in l1:
        if i:
            n+=1
            l2.append(1)
        else:
            l2.append(n)
            n=0
    l2=np.array(l2)
    return([l2.argmax()-l2.max(),l2.argmax()])

def getGoodRegion(inQzv,o='forward',q=20):
    print('Parsing '+inQzv)
    z=zipfile.ZipFile(inQzv)
    f=[i for i in z.filelist if re.search(o+'.*summaries.tsv',i.filename.split('/')[-1])]
    df1=pd.read_csv(z.open(f[0]),sep='\t',header=0,index_col=0,comment='#')
    med=df1.loc['50%',:]
    med1=med.rolling(5).mean()
    med1[:4]=med[:4]
    med1=med1>=q
    r=getRegion(med1)
    print('Will trim sequences at '+str(r[0])+' bp and '+str(r[1])+' bp and remove other regions with moving average of median quality score < '+str(q)+'\nPlease ensure there will be at least 12 bp overlap after trimming for paired-end sequence data')
    return(r)

def multiIndexDataframe(df,name,axis=0):
    df1=df.copy()
    if axis==0:
        row=list(zip([name]*df1.shape[0],df1.index))
        ind=pd.MultiIndex.from_tuples(row)
        df1.index=ind
        return(df1)
    else:
        col=list(zip([name]*df1.shape[1],df1.columns))
        ind=pd.MultiIndex.from_tuples(col)
        df1.columns=ind
        return(df1)

def readProfileDemux(indir,out):
    print('Summarising read number after demultiplexing of all projects based on '+indir)
    df=pd.DataFrame()
    for i in os.listdir(indir):
        if not i.split('.')[-1]=='qzv':
            continue
        project=i[:-4]
        z=zipfile.ZipFile(indir+'/'+i)
        f=[i for i in z.filelist if re.search(r'counts.tsv$',i.filename.split('/')[-1])]
        df1=pd.read_csv(z.open(f[0]),sep='\t',header=0,index_col=0,comment='#')
        r=df1.iloc[:,0]
        print('\n'+project+':\n'+'# sample: '+str(len(r))+'\nMax read: '+str(np.max(r))+'\nMin read: '+str(np.min(r))+'\nMedian read: '+str(np.median(r))+'\nMean read: '+str(round(np.mean(r))))
        df1= multiIndexDataframe(df1,project)
        df=pd.concat([df,df1],sort=False)

    df.to_csv(out,sep='\t')
    print('\nResults saved to '+out)

def denoise(inQza,outdir,q=20,p=10):
    print('Processing '+inQza)
    if os.path.basename(inQza).split('-')[0][-2:]=='Se':
        project=os.path.basename(inQza)[:-10]
        qzv=os.path.dirname(inQza)+'/../demuxSumViz/'+project+'.qzv'
        trim=getGoodRegion(qzv,o='forward',q=q)
        call('qiime dada2 denoise-single \
                --i-demultiplexed-seqs '+inQza+'\
                --p-trim-left '+str(trim[0])+' \
                --p-trunc-len '+str(trim[1])+' \
                --p-trunc-q '+str(q)+' \
                --p-n-threads '+str(p)+' \
                --o-representative-sequences '+outdir+'/'+project+'-rep-seqs.qza \
                --o-table '+outdir+'/'+project+'-table.qza \
                --o-denoising-stats '+outdir+'/'+project+'-stats.qza'\
                ,'Abundance results saved in '+outdir+'.','Abundance summarization failed.')
    elif os.path.basename(inQza).split('-')[0][-2:]=='Pe':
        project=os.path.basename(inQza)[:-10]
        qzv=os.path.dirname(inQza)+'/../demuxSumViz/'+project+'.qzv'
        trimF=getGoodRegion(qzv,o='forward',q=q)
        trimR=getGoodRegion(qzv,o='reverse',q=q)
        call('qiime dada2 denoise-paired \
                --i-demultiplexed-seqs '+inQza+'\
                --p-trim-left-f '+str(trimF[0])+' \
                --p-trunc-len-f '+str(trimF[1])+' \
                --p-trim-left-r '+str(trimR[0])+' \
                --p-trunc-len-r '+str(trimR[1])+' \
                --p-trunc-q '+str(q)+' \
                --p-n-threads '+str(p)+' \
                --o-representative-sequences '+outdir+'/'+project+'-rep-seqs.qza \
                --o-table '+outdir+'/'+project+'-table.qza \
                --o-denoising-stats '+outdir+'/'+project+'-stats.qza'\
                ,'Abundance results saved in '+outdir+'.','Abundance summarization failed. Is there an overlap of forward and reverse reads after trimming?')
    else:
        print('Error: Unknown demultiplexed data input.')
        sys.exit(1)

def readProfileDenoise(indir,out):
    print('Summarising read number after denoising of all projects based on '+indir)
    df=pd.DataFrame()
    for i in os.listdir(indir):
        if not i.split('-')[-1]=='stats.qza':
            continue
        project=i[:-10]
        call('qiime metadata tabulate --m-input-file '+indir+'/'+i+' --o-visualization '+indir+'/'+project+'-stats.qzv','Tabulating '+indir+'/'+i+' done successfully.','Tabulating '+indir+'/'+i+' failed.')
        z=zipfile.ZipFile(indir+'/'+project+'-stats.qzv')
        f=[i for i in z.filelist if re.search(r'metadata.tsv$',i.filename.split('/')[-1])]
        df1=pd.read_csv(z.open(f[0]),sep='\t',header=0,index_col=0,comment='#')
        df1= multiIndexDataframe(df1,project)
        df=pd.concat([df,df1],sort=False)

    df.to_csv(out,sep='\t')
    print('Results saved to '+out)

def convertInFt(indir,outdir):
    print('Converting input feature table to qza')
    for i in os.listdir(indir):
        if not i.split('-')[0]=='featureTable':
            continue
        print('processing '+indir+'/'+i)

        for j in os.listdir(indir+'/'+i):
            if j=='feature-table.tsv' or j=='feature-table.biom':
                call('biom convert -i '+indir+'/'+i+'/'+j+' -o '+outdir+'/'+i+'-feature-table.biom --to-hdf5','Converted to format biom version 2.1 hdf5 as '+outdir+'/'+i+'-feature-table.biom.','Converting to format biom version 2.1 hdf5 failed.')
                call('qiime tools import --input-path '+outdir+'/'+i+'-feature-table.biom --type FeatureTable[Frequency] --input-format BIOMV210Format --output-path '+outdir+'/'+i+'-table.qza; rm -rf '+outdir+'/'+i+'-feature-table.biom','Biom imported as qza','Importing biom as qza failed.')
            elif j=='rep-seqs.fasta':
                call('qiime tools import --input-path '+indir+'/'+i+'/'+j+' --output-path '+outdir+'/'+i+'-rep-seqs.qza --type FeatureData[Sequence]','Representative sequence imported as qza as '+outdir+'/'+i+'-rep-seqs.qza.','Importing representative sequence as qza failed.')

def join_pe(indir, outdir):
    call('rm -rf '+outdir+'; mkdir '+outdir, outdir+' created', outdir+' not created')
    for i in os.listdir(indir):
        if i[-9:] == 'demux.qza' and i.split('-')[0][-2:] == 'Pe':
            call('qiime vsearch join-pairs --i-demultiplexed-seqs %s/%s --o-joined-sequences %s/%s' % (indir, i, outdir, i), i+' paired end reads joined', i+' paired end reads joining failed')
        elif i[-9:] == 'demux.qza' and i.split('-')[0][-2:] == 'Se':
            call('cp %s/%s %s/%s' % (indir, i, outdir, i))

def quality_filter(indir, score, outdir, statsdir):
    call('rm -rf '+outdir+'; mkdir '+outdir, outdir+' created', outdir+' not created')
    for i in os.listdir(indir):
        if i[-9:] != 'demux.qza':
            continue
        call('qiime quality-filter q-score --i-demux %s --p-min-quality %s --o-filtered-sequences %s --o-filter-stats %s' % (indir+'/'+i, score, outdir+'/'+i, statsdir+'/'+i+'_filtered'), 'Quality filtering done', 'Quality filtering failed')

def otu_clustering(indir, identity, outdir, cpu):
    call('rm -rf '+outdir+'; mkdir '+outdir, outdir+' created', outdir+' not created')
    tempdir = outdir+'/temp'
    call('mkdir '+tempdir, tempdir + ' created', tempdir+' not created')

    def read_demux_qza(infile, outfile, tempdir): # extract qza, convert fastq files of each sample into fasta with sample IDs, merge to a single file seqs.fna
        z = zipfile.ZipFile(infile, 'r')
        dirname = z.namelist()[0].split('/')[0]
        z.extractall(tempdir)
        df = pd.read_csv(tempdir+'/'+dirname+'/data/MANIFEST', sep=',', header=0, index_col=0, comment='#')
        outfile = open(outfile, 'a')
        for i in df.index:
            fastq = tempdir+'/'+dirname+'/data/'+df.loc[i, 'filename']
            n = 0
            for rec in SeqIO.parse(gzip.open(fastq, 'rt'), 'fastq'):
                rec.id = str(i)+'_'+str(n)
                SeqIO.write(rec, outfile, 'fasta')
                n += 1
        outfile.close()
        print(infile+' written into fasta')
    
    # process Pe and Se demux.qza and get fasta file seqs.fna
    for i in os.listdir(indir):
        if i[-9:] != 'demux.qza':
            continue
        read_demux_qza(indir+'/'+i, outdir+'/seqs.fna', tempdir)
    os.system('rm -rf '+tempdir)

    # dereplicate
    def remove_new_line(infile):
        f = ''
        n = {}
        for i in open(infile, 'r'):
            if i[0] == '>':
                sample = i.split(' ')[0]
                sample = sample[1:]
                sample = '_'.join(sample.split('_')[:-1])
                if n.get(sample):
                    n[sample] += 1
                else:
                    n[sample] = 1
                f += '\n>'+sample+'_'+str(n[sample])+'\n'
            else:
                f += i[:-1]
        f = f[1:]
        f += '\n'

        outfile = open(infile, 'w')
        outfile.write(f)
        outfile.close()

    remove_new_line(outdir+'/seqs.fna')
    call('qiime tools import --input-path %s/seqs.fna --output-path %s/seqs.qza --type \'SampleData[Sequences]\'' % (outdir, outdir), 'Importing seqs.fna done', 'Importing seqs.fna failed')
    
    call('qiime vsearch dereplicate-sequences --i-sequences %s/seqs.qza --o-dereplicated-table %s/table.qza --o-dereplicated-sequences %s/rep-seqs.qza' % (outdir, outdir, outdir), 'Dereplicate done', 'Dereplicate failed')

    # OTU clustering
    call('qiime vsearch cluster-features-de-novo --i-table %s/table.qza --i-sequences %s/rep-seqs.qza --p-perc-identity %s --o-clustered-table %s/table-dn.qza --o-clustered-sequences %s/rep-seqs-dn.qza --p-threads %s' % (outdir, outdir, identity, outdir, outdir, cpu), 'OTU clustering done', 'OTU clustering failed')

def otu_read_summary(infile, outfile):
    indir = os.path.dirname(infile)
    z = zipfile.ZipFile(infile, 'r')
    dirname = z.namelist()[0].split('/')[0]
    z.extract(dirname+'/data/feature-table.biom', indir)
    os.system('mv '+indir+'/'+dirname+'/data/feature-table.biom '+indir)
    os.system('rm -rf '+indir+'/'+dirname)
    call('biom convert -i %s/feature-table.biom -o %s/feature-table.txt --to-tsv' % (indir, indir))
    df = pd.read_csv(indir+'/feature-table.txt', sep='\t', header=1, index_col=0)
    df.sum().to_csv(outfile, sep=',')

def mergeProjects(indir,outdir):
    table=[i for i in os.listdir(indir) if i[-9:]=='table.qza']
    seqs=[i for i in os.listdir(indir) if i[-12:]=='rep-seqs.qza']

    tableOpt=['--i-tables '+indir+'/'+i for i in table]
    seqsOpt=['--i-data '+indir+'/'+i for i in seqs]
    
    print('Merging '+str(table)+'to '+outdir)
    call('qiime feature-table merge '+' '.join(tableOpt)+' --o-merged-table '+outdir+'/merged-table.qza --p-overlap-method sum',str(table)+' from '+indir+'/ were merged to '+outdir+'/merged-table.qza','Merging multiple table.qza failed.')
    call('qiime feature-table summarize --i-table '+outdir+'/merged-table.qza --o-visualization '+outdir+'/merged-table.qzv','Summarizing merged table done successfully.','Summarizing merged table failed.')
    call('qiime tools export --input-path '+outdir+'/merged-table.qza --output-path '+outdir,'Exporting merged table done successfully.','Exporting merged table failed.')
    call('biom convert -i '+outdir+'/feature-table.biom -o '+outdir+'/feature-table.tsv --to-tsv','Converting biom to tsv done successfully.','Converting biom to tsv failed.')
    call('mv '+outdir+'/feature-table.biom '+outdir+'/merged-table.biom; mv '+outdir+'/feature-table.tsv '+outdir+'/merged-table.tsv','Renaming to merged-table.biom and merged-table.tsv done successfully.','Renaming to merged-table.biom and merged-table.tsv failed.')

    print('Merging '+str(seqs)+'to '+outdir)
    call('qiime feature-table merge-seqs '+' '.join(seqsOpt)+' --o-merged-data '+outdir+'/merged-rep-seqs.qza',str(seqs)+' from '+indir+'/ were merged to '+outdir+'/merged-rep-seqs.qza','Merging multiple rep-seq.qza failed.')
    call('qiime feature-table tabulate-seqs --i-data '+outdir+'/merged-rep-seqs.qza --o-visualization '+outdir+'/merged-rep-seqs.qzv','Summarizing merged rep-seqs done successfully.','Summarizing merged rep-seqs failed.')
    call('qiime tools export --input-path '+outdir+'/merged-rep-seqs.qza --output-path '+outdir,'Exporting merged rep-seqs done successfully.','Exporting rep-seqs failed.')
    call('mv '+outdir+'/dna-sequences.fasta '+outdir+'/merged-rep-seqs.fasta','Renaming to merged-table.biom and merged-table.tsv done successfully.','Renaming to merged-table.biom and merged-table.tsv failed.')

def getFileFromZ(zipFile,target,out):
    z=zipfile.ZipFile(zipFile)
    f=[i for i in z.filelist if i.filename.split('/')[-1]==target]
    open(out,'w').write(z.read(f[0]).decode())

def resampleDepth(infile,resa):
    # This function return resampling depth yielding the maximum observation
    print('Estimating the optimal resampling depth yielding the maximum observation (# samples and # features)')
    r=[round(float(i.strip().split(',')[1])) for i in open(infile,'r').readlines()]
    r.pop(0)
    r.sort()
    print('# sample: '+str(len(r))+'\nMax read: '+str(np.max(r))+'\nMin read: '+str(np.min(r))+'\nMedian read: '+str(np.median(r))+'\nMean read: '+str(round(np.mean(r))))
    if resa != 0:
        print('Using user specified resampling depth '+str(resa)+'.\n# sample left: '+str((np.array(r) > resa).sum())+'\n')
        return(resa)

    s=list(range(1,len(r)+1))
    s=s[::-1]
    h=round(len(r)/2)
    r1=r[:h]
    s1=s[:h]
    t1=np.array(r1)*np.array(s1)
    print('Resampling depth: '+str(r1[t1.argmax()])+'\n# sample left: '+str(s1[t1.argmax()])+'\n')
    return(r1[t1.argmax()])

def updateDataFrame(dfList):
    df=pd.concat(dfList,sort=False)
    df=df.loc[~df.index.duplicated(),:]
    for i in dfList:
        df.update(i)
    df.sort_index(inplace=True)
    df.sort_index(axis=1,inplace=True)
    return(df)

def mergeMetadata(indir,out):
    metadataFile=[indir+'/'+i+'/metadata.tsv' for i in os.listdir(indir) if i[0] != '.']
    df=[pd.read_csv(i,sep='\t',header=0,index_col=0,comment='#') for i in metadataFile]
    df=updateDataFrame(df)
    df.to_csv(out,sep='\t')
    print('Merged metadata saved to '+out)

def alphaGroupSig(indir,metadataFile):
    df=pd.read_csv(metadataFile,sep='\t',header=0,index_col=0)
    if len([i for i in df.columns if len(df[i].unique())>1])>0:
        print('Conducting group significance test using categorial factors for alpha diversity results in '+indir)
        qza=[i for i in os.listdir(indir) if i[-10:]=='vector.qza']
        for i in qza:
            try:
                call('qiime diversity alpha-group-significance --i-alpha-diversity '+indir+'/'+i+' --m-metadata-file '+metadataFile+' --o-visualization '+indir+'/'+i[:10]+'-groupSignificance.qzv',i+' done successfully.',i+' failed')
            except:
                print('Will not conduct group significance test using categorial factors for alpha diversity results in '+indir+' because no qualified columns in the metadata are found.')
    else:
        print('Will not conduct group significance test using categorial factors for alpha diversity results in '+indir+' because no qualified columns in the metadata are found.')

def alphaCorrelation(indir,metadataFile,cor='spearman'):
    df=pd.read_csv(metadataFile,sep='\t',header=0,index_col=0)
    if len([i for i in df.columns if len(df[i].unique())>1])>0:
        print('Conducting correlation test using continuous factors for alpha diversity results in '+indir)
        qza=[i for i in os.listdir(indir) if i[-10:]=='vector.qza']
        for i in qza:
            try:
                call('qiime diversity alpha-correlation --i-alpha-diversity '+indir+'/'+i+' --m-metadata-file '+metadataFile+' --p-intersect-ids --p-method '+cor+' --o-visualization '+indir+'/'+i[:10]+'-correlation.qzv',i+' done successfully.',i+' failed')
            except:
                print('Will not conduct group significance test using categorial factors for alpha diversity results in '+indir+' because no qualified columns in the metadata are found.')
    else:
        print('Will not conduct correlation test using continuous factors for alpha diversity results in '+indir+' because no qualified columns in the metadata are found.')

def classification(seq,table,classifier,metadata,outdir,conf,cpu):
    print('Conducting taxonomic classification')
    call('qiime feature-classifier classify-sklearn --i-classifier '+classifier+' --i-reads '+seq+' --o-classification '+outdir+'/classification.qza'+' --p-confidence '+str(conf)+' --p-n-jobs '+str(cpu),'Classification done successfully.','Classification failed.')
    print('Converting to classification.qzv')
    call('qiime metadata tabulate --m-input-file '+outdir+'/classification.qza --o-visualization '+outdir+'/classification.qzv','Conversion done successfully.','Conversion failed.')
    getFileFromZ(outdir+'/classification.qza','taxonomy.tsv',outdir+'/taxonomy.tsv')
    call('qiime taxa barplot --i-table '+table+' --i-taxonomy '+outdir+'/classification.qza --m-metadata-file '+metadata+' --o-visualization '+outdir+'/taxa-bar-plots.qzv','Bar plot generated.','Bar plot failed.')

def renameTax(tax):
    tax=re.sub(r';__','',tax)
    tax=tax.split(';')[-1]
    return(tax)

def splitTax(infile,outfile):
    df1=pd.read_csv(infile,header=0,index_col=0)
    tax=[i for i in df1.columns if re.search(r'^[a-z]__',i)]
    df1=df1.loc[:,tax]
    df1.columns=[renameTax(i) for i in df1.columns]
    df1.to_csv(outfile,sep='\t')

def getTax(inQzv,outdir):
    z=zipfile.ZipFile(inQzv)
    f=[i for i in z.filelist if re.search(r'level-\d+.csv',i.filename.split('/')[-1])]
    for i in f:
        splitTax(z.open(i),outdir+'/tax-L'+i.filename[-5]+'.tsv')
        print(outdir+'/tax-L'+i.filename[-5]+'.tsv saved.')

def ordination(inQzv,metadataFile,outdir):
    home=os.path.dirname(os.path.realpath(__file__))

    print('Extracting taxonomy data.')
    getTax(inQzv,outdir)

    print('Conducting ordination analysis.')
    tax=[i for i in os.listdir(outdir) if re.search(r'tax-L\d+.tsv$',i)]
    tax.remove('tax-L1.tsv')
    for i in tax:
        call('Rscript '+home+'/../lib/ordination.r '+outdir+'/'+i+' '+metadataFile+' '+outdir+'/RDA-L'+i[5]+'.txt '+outdir+'/RDA-L'+i[5]+'.pdf '+outdir+'/CCA-L'+i[5]+'.txt '+outdir+'/CCA-L'+i[5]+'.pdf '+outdir+'/log.txt','Ordination analysis RDA, CCA and ordistep tests finished for '+i,'Ordination analysis RDA, CCA and ordistep tests failed for '+i)

