import subprocess
from pathlib import Path
import pandas as pd
from collections import defaultdict
import os
import pysam
import re
import gffutils

def run_STAR(fastq_fn, num_threads=4, num_reads_to_map=-1):
    ''' Take a single fastq and run STAR to align and get gene counts, and put the output results in its own subdirectory. There cant be any dots (.) in the file name because of how the name is split.'''
    
    fastq_fn = Path(fastq_fn)
    
    prefix = fastq_fn.name.split('.')[0]

    output_dir = Path(fastq_fn.parent / 'STAR_output_files')
    output_dir.mkdir(exist_ok=True)
    
    full_prefix = output_dir / (prefix + '_') 
    
    STAR_command = [
        'STAR',
        '--runThreadN', str(num_threads),
        '--runMode', 'alignReads',
        '--outSAMtype', 'BAM', 'SortedByCoordinate',
        '--readFilesCommand', 'zcat',
        '--genomeDir', '/nvme/indices/refdata-cellranger-mm10-1.2.0/star',
        '--outFileNamePrefix', str(full_prefix),
        '--readFilesIn', str(fastq_fn),
        '--quantMode', 'GeneCounts',
        '--readMapNumber', str(num_reads_to_map),
    ]
    
    print(fastq_fn)
    subprocess.check_call(STAR_command)
    
    output_bam = str(full_prefix) + 'Aligned.sortedByCoord.out.bam'
    
    pysam.index(output_bam)
    

def run_kallisto(fastq_fn):
    '''Takes a single fastq and runs kallisto to generate abundance files which contain tpms for transcripts (ENMUST), output results in it's own subdirectory. There cant be any dots (.) in the file name because of how the name is split.'''

    fastq_fn = Path(fastq_fn) #creating a path object of the fastq to manipulate
    prefix = fastq_fn.name.split('.')[0]     
    current_dir = fastq_fn.parent

    kallisto_command = [
        'kallisto',
        'quant',
        '--index=/home/lgoodrich/RNAseq/kallisto/mm10_transcripts.idx',
        '--output-dir={}/kallisto_'.format(str(current_dir)) + str(prefix), 
        str(fastq_fn),
        '--single', '-l', '200', '-s', '30' 
        #there isn't an option to only do a few thousand reads
    ]

    print(fastq_fn)
    subprocess.check_call(kallisto_command)

def transcripts_to_genes(input_tsv, new_csv_file_name='collapsed_tpms'):
    '''Starting with the abundance tsv file from Kallisto as input, collapse transcripts into genes (ENMUST to ENMUSG). Save as a new csv in the same directory with the file name specified.
    
    Currently is just outputting the index, not the whole dataframe'''
    
    input_path = Path(input_tsv)
    output_dir = Path(input_path.parent)

    abundance_df = pd.read_csv(input_tsv, sep='\t', index_col=0)

    tpm_df = abundance_df['tpm']

    transcript_to_gene_df = pd.read_csv('/home/lgoodrich/igv/Mus_musculus.GRCm38.90_transcript_ID_to_gene_ID.csv', index_col=0)

    truncated_ID = [i.split('.')[0] for i in tpm_df.index]
    tpm_df.index = truncated_ID

    gene_to_transcripts = defaultdict(list)

    for row in transcript_to_gene_df.itertuples():
        transcript = row[0]
        gene = row[1]
        #check if the transcript is present in the tpm_df
        #if it is, add it, if not, continue
        if transcript in tpm_df.index:
            gene_to_transcripts[gene].append(transcript)

    gene_to_summed_tpms = {}

    for gene, transcripts in gene_to_transcripts.items():
        gene_to_summed_tpms[gene] = tpm_df.loc[transcripts].sum()

    collapsed_df = pd.DataFrame(gene_to_summed_tpms, index=['tpms']).T

    collapsed_df.to_csv(str(output_dir) + '/{}.csv'.format(new_csv_file_name))

def coord_to_counts_pysam(list_coords, bam_fn):
    '''Given a set of coordinates in a list with the format [chr, start, end] and an indexed bam file, return the counts as tuples, in the form (name, counts)'''
    
    bam_fn_path = Path(bam_fn)
    
    bam_name = bam_fn_path.name.split('.')[0]
    
    samfile = pysam.AlignmentFile(bam_fn, 'rb')
    
    counts = samfile.count(list_coords[0], list_coords[1], list_coords[2])
    
    return bam_name, counts

def coords_for_pysam(db, ensembl_ID):
    '''Returns coordinates in the right order to be read by pysam: [chromosome, start, stop].'''
    
    coord = []
    
    attributes = ['chrom', 'start', 'stop']
    for attribute in attributes:
        coord.append(rnaseq_experiment.lookup_attributes_from_gffutils(db, ensembl_ID, attribute))
    
    return coord

def make_dict_counts(bam_files, list_coords, name_of_feature):
    '''Given a list of bam files, and coordinates for a feature of interest, return a single dictionary of bam_fn: count.'''
    
    name_counts = {}
    for bam_file in sorted(bam_files):
        bam_name, counts = coord_to_counts_pysam(list_coords, bam_file)
        name_counts[bam_name] = counts
    
    return name_counts
    

def make_df_coords_counts(bam_files, list_coords, name_of_feature):
    '''Given a list of bam files, and coordinates for a feature of interest, return a pandas dataframe with a single column named after the feature, the index as the bam file names, and values as counts.
    '''
    
    master_dict = {}
    
    name_counts = {}
    for bam_file in sorted(bam_files):
        bam_name, counts = coord_to_counts_pysam(list_coords, bam_file)
        name_counts[bam_name] = counts
    
    master_dict[name_of_feature] = name_counts
    
    df = pd.DataFrame.from_dict(master_dict, orient='index').T
        
    return df

def sample_name_from_STAR_log(string):
    '''Function to get sample name from kallisto output directories. Unique to time course data that has 'day' in the file name
    '''
    match_object = re.search('(.*_day_\w)', string) 

    if match_object is not None:        
        sample_name = match_object.group(0) 

    else: 
        raise ValueError

    return (sample_name)   

def sample_info(string):
    '''add description'''
    sorted_mo = re.search('(.*)_(sorted_day_\w)', string)
    mo = re.search('(.*)_(day_\w)', string) #separates string into two capture groups

    if sorted_mo is not None:            
        sample_name = sorted_mo.group(1)
        sample_day = sorted_mo.group(2)      
       
    elif mo is not None:        
        sample_name = mo.group(1) #everything before "day" is the sample name
        sample_day = mo.group(2)
        
    else: 
        raise ValueError
    
    return (sample_name, sample_day)   

def ref_counts_from_idxstats(string):
    '''Given output of samtools idxstats for alignment turn the ref:mapped_reads into a dictionary'''
    ref_counts = {}
    
    list_of_counts = string.splitlines()
    for line in list_of_counts:
        ref, ref_length, num_mapped_reads, _ = line.split()
        ref_counts[ref] = num_mapped_reads
        
    return(ref_counts)

def make_sample_ref_counts_dictionary(bam_fns):
    
    nested_dict = {}
    
    for bam_fn in bam_fns:
        bam_path = Path(bam_fn)
        file_name = bam_path.name
        sample_name = sample_name_from_STAR_log(str(file_name)) #key of dict
        data = pysam.idxstats(bam_fn)
        nested_dict[sample_name] = ref_counts_from_idxstats(data) #value of dict

    return nested_dict

def STAR_metrics(metric_fns):
    '''metric file format: *Log.final.out'''
    
    master_dict = {}

    for fn in metric_fns:
        head, tail = os.path.split(fn)
        sample_name = sample_name_from_STAR_log(tail)

        master_dict[sample_name] = {}

        metrics = []

        with open(fn) as file:
            for line in file:
                metrics.append(line.strip())

            metrics_values = {}

            for line in metrics:

                num_input_reads_match_object = re.search('Number of input reads \|\t(\d*)', line)
                num_unique_mapped_match_object = re.search('Uniquely mapped reads number \|\t(\d*)', line)
                percent_match_object = re.search('Uniquely mapped reads % \|\t(\d*.\d*)', line)

                if num_input_reads_match_object is not None: 
                    input_reads = int((num_input_reads_match_object).group(1))
                    metrics_values['num_input_reads'] = input_reads         

                elif num_unique_mapped_match_object is not None: #if there is a match for this metric
                    num_unique_mapped_reads = int((num_unique_mapped_match_object).group(1))
                    metrics_values['num_unique_mapped_reads'] = num_unique_mapped_reads #append it to our dictionary

                elif percent_match_object is not None: 
                    percent_mapped = float((percent_match_object).group(1))
                    metrics_values['percent_unique_mapped_reads'] = percent_mapped 

            master_dict[sample_name] = metrics_values

    return master_dict  

def lookup_attributes_from_gffutils(db, name, attr, special=False):
    '''Must first call the gffutils database with .FeatureDB
    
    Special attributes are name, description, biotype.'''
    
    #db = gffutils.FeatureDB('/home/lgoodrich/igv/mm10_gff.db')
    
    full_name = 'gene:' + name
    try:
        feature = db[full_name]
        if special:
            return feature.attributes[attr][0]
        else:
            return getattr(feature, attr)
    except:
        return None
    
def add_attributes_to_df(gff_file, df):
    '''Uses gffutils to lookup and add attributes for a particular ENSEMBL_ID. ID must be the index. Adds the following attributes: name, description, biotype, chromosome, strand, start, stop.'''
    
    db = gffutils.FeatureDB(gff_file)
    
    data_columns = df.columns
        
    df['name'] = [lookup_attributes_from_gffutils(db, name, 'Name', special=True) for name in df.index] 
    df['description'] = [lookup_attributes_from_gffutils(db, name, 'description', special=True) for name in df.index] 
    df['biotype'] = [lookup_attributes_from_gffutils(db, name, 'biotype', special=True) for name in df.index] 
    df['chromosome'] = [lookup_attributes_from_gffutils(db, name, 'chrom') for name in df.index] 
    df['strand'] = [lookup_attributes_from_gffutils(db, name, 'strand') for name in df.index] 
    df['start'] = [lookup_attributes_from_gffutils(db, name, 'start') for name in df.index] 
    df['stop'] = [lookup_attributes_from_gffutils(db, name, 'stop') for name in df.index] 
    attribute_cols =['name', 'description', 'biotype', 'chromosome', 'strand', 'start','stop']
    
    
    
    full_df = df[attribute_cols + list(data_columns)]
    
    return full_df

def genes_over_TPM_cutoff(df, cutoff):
    '''Given a pandas dataframe with column names ending in _tpm, return a list of genes (indices) over a given cutoff. 1 is a good starting number.'''
    samples = [c for c in df.columns if '_tpm' in c]
    just_tpms = df[samples]
    means = just_tpms.mean(axis=1)
    above_cutoff = means[means > cutoff]
    return list(above_cutoff.index)
