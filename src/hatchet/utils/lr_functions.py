import subprocess
import numpy as np
import pandas as pd
from hatchet.utils.ArgParsing import parse_combine_counts_args, parse_genotype_snps_arguments
from hatchet.utils.Supporting import log, logArgs
from hatchet.utils.combine_counts import adaptive_bins_arm, read_snps, read_total_and_thresholds


def genotype_snps(args):
    """
    This functions creates extracts each chromosome from the genotype file,
    and saves it to a separate file in the output directory.
    """
    log(
        msg=(
            '# Parsing the input arguments, checking the consistency of given files, and extracting required ',
            'information\n',
        ),
        level='STEP',
    )
    args = parse_genotype_snps_arguments(args)
    logArgs(args, 80)

    for chrom in args['chromosomes']:
        log(
            msg=('Extracting chromosome {} from genotype file...\n'.format(chrom),),
            level='STEP',
        )

        # call bcftools to extract chromosome
        output_dir = args['outputsnps']
        output_vcf = f'{output_dir}/{chrom}.vcf.gz'
        mindp = args['mincov']
        maxdp = args['maxcov']
        command = [
            args['bcftools'],
            'view',
            '-r',
            chrom,
            '-i',
            f'FILTER="PASS" & strlen(REF)=1 & strlen(ALT)=1 & SUM(FORMAT/AD)>={mindp} & SUM(FORMAT/AD)<={maxdp} & N_ALT=1 & (FORMAT/GT[0]=="1|0" | FORMAT/GT[0]=="0|1")',
            args['snps'],
            '-Oz',
            '-o',
            output_vcf,
        ]
        process = subprocess.Popen(command)
        process.wait()  # Wait for the process to complete


def phase_snps(args):
    """
    This functions creates extracts each phased genotype from genotpes
    file and saves it to phased.vcf.gz in the output directory.
    """
    log(
        msg=(
            '# Parsing the input arguments, checking the consistency of given files, and extracting required ',
            'information\n',
        ),
        level='STEP',
    )
    args = parse_genotype_snps_arguments(args)
    logArgs(args, 80)

    # call bcftools to extract 0|1 and 1|0 snps
    output_dir = args['outputsnps']
    output_vcf = f'{output_dir}/phased.vcf.gz'
    mindp = args['mincov']
    maxdp = args['maxcov']
    command = [
        args['bcftools'],
        'view',
        '-i',
        f'FILTER="PASS" & strlen(REF)=1 & strlen(ALT)=1 & SUM(FORMAT/AD)>={mindp} & SUM(FORMAT/AD)<={maxdp} & N_ALT=1 & (FORMAT/GT[0]=="1|0" | FORMAT/GT[0]=="0|1")',
        args['snps'],
        '-Oz',
        '-o',
        output_vcf,
    ]
    process = subprocess.Popen(command)
    process.wait()  # Wait for the process to complete


def get_b_count(df):
    # Select REF if FLIP == 1, otherwise select ALT, then sum the selected values
    sel = df.apply(lambda row: row['REF'] if row['FLIP'] == 1 else row['ALT'], axis=1)
    total_sum = sel.sum()
    return total_sum


def get_haplostring(df):
    return ''.join(list(df['FLIP'].astype(int).astype(str)))


def combine_counts(args, haplotype_file):
    log(msg='# Parsing and checking input arguments\n', level='STEP')
    args = parse_combine_counts_args(args)
    logArgs(args, 80)

    baffile = args['baffile']
    threads = args['processes']
    chromosomes = args['chromosomes']
    outfile = args['outfile']
    all_names = args['sample_names']
    msr = args['min_snp_reads']
    mtr = args['min_total_reads']
    use_chr = args['use_chr']
    phase = args['phase']
    blocksize = args['blocksize']
    max_snps_per_block = args['max_snps_per_block']
    test_alpha = args['test_alpha']
    multisample = args['multisample']
    nonormalFlag = args['nonormalFlag']
    ponfile = args['ponfile']
    referencefasta = args['referencefasta']
    XX = args['XX']
    rd_array = args['array']

    def read_ranges_from_file(filename):
        '''
        Reads haplotype blocks from gtf file
        '''
        df = pd.read_csv(filename, sep='\t', header=None, comment='#',
                        names=['CHR', 'source', 'feature', 'START', 'END', 'score', 'strand', 'frame', 'attribute'])
        ranges = df[['CHR', 'START', 'END']]
        return ranges
    
    # Read ranges from the haplotype file
    haplotype_blocks = read_ranges_from_file(haplotype_file)

    rows = []
    for ch in chromosomes:
        positions, snp_counts, snpsv = read_snps(baffile, ch, all_names, phasefile=phase)
        total_counts, complete_thresholds = read_total_and_thresholds(ch, rd_array, False)
        blocks = haplotype_blocks[haplotype_blocks.CHR == ch]

        for index, row in blocks.iterrows():
            block_snps_idx = np.where((positions >= row.START) & (positions <= row.END))[0]
            block_snp_pos = positions[block_snps_idx]
            block_snp_counts = snp_counts[block_snps_idx]
            if len(block_snp_counts) == 0:
                continue
            block_thr_idx = np.where((complete_thresholds >= block_snp_pos[0]) & (complete_thresholds <= block_snp_pos[-1]))[0]
            if len(block_snps_idx) > len(block_thr_idx) + 1 or len(block_thr_idx) == 0:
                # centromere loci
                continue
            log(msg=f'snpcount {len(block_snps_idx)} thrcount {len(block_thr_idx)}\n', level='STEP')
            
            # thresholds must begin and end outside the boundaries:
            block_thr_idx = np.insert(block_thr_idx, 0, block_thr_idx[0] - 1)
            block_thr_idx = np.append(block_thr_idx, block_thr_idx[-1] + 1)

            block_thr_pos = complete_thresholds[block_thr_idx]
            block_thr_counts = total_counts[block_thr_idx]

            bins = adaptive_bins_arm(
                snp_thresholds=block_thr_pos,
                total_counts=block_thr_counts,
                snp_positions=block_snp_pos,
                snp_counts=block_snp_counts,
                chromosome=ch,
                min_snp_reads=msr,
                min_total_reads=1,
                nonormalFlag=nonormalFlag,
                use_averages_rd=True
            )

            starts = bins[0]
            ends = bins[1]

            dfs = [snpsv[(snpsv.POS >= starts[i]) & (snpsv.POS <= ends[i])] for i in range(len(starts))]
            dfs = [df.dropna(subset=['FLIP']) for df in dfs]

            for i, df in enumerate(dfs):
                if len(df) == 0:
                    continue
                for sample, df2 in df.groupby('SAMPLE'):
                    start = starts[i]
                    end = ends[i]
                    total = df2.TOTAL.sum()
                    bcount = get_b_count(df2)
                    bcount = min(bcount, total - bcount)
                    haplostring = get_haplostring(df2)
                    log(msg=f'{ch} {sample} {start} {end} {total} {bcount} {bcount/total}\n', level='STEP')
                    rows.append(
                        [
                            ch,
                            "unit",
                            start,
                            end,
                            sample,
                            bins[3][i][0],
                            bins[2][i][1],
                            bins[2][i][0],
                            len(df2),
                            bcount,
                            total,
                            "",
                            "",
                            "",
                            "",
                            bcount / total,
                        ]
                    )
            # snp_dfs = snpsv[(snpsv.POS >= row.START) & (snpsv.POS <= row.END)]
            # if len(snp_dfs) == 0:
            #     continue

        big_bb = pd.DataFrame(
            rows,
            columns=[
                'CHR',
                'UNIT',
                'START',
                'END',
                'SAMPLE',
                'RD',
                'TOTAL_READS',
                'NORMAL_READS',
                'SNPS',
                'BCOUNT',
                'TOTAL_SNP_READS',
                'HAPLO',
                'SNP_POS',
                'SNP_REF_COUNTS',
                'SNP_ALT_COUNTS',
                'BAF',
            ],
        )

        big_bb['CORRECTED_READS'] = np.NAN
        # For each sample, correct read counts to account for differences in coverage (as in HATCHet)
        # (i.e., multiply read counts by total-reads-normal/total-reads-sample)
        rc = pd.read_table(args['totalcounts'], header=None, names=['SAMPLE', '#READS'])
        normal_name = all_names[0]
        nreads_normal = rc[rc.SAMPLE == normal_name].iloc[0]['#READS']
        if nonormalFlag:
            for sample, df in big_bb.groupby('SAMPLE'):
                big_bb.loc[big_bb.SAMPLE == sample, 'RD'] = (
                    df['TOTAL_READS'] / (df['END'] - df['START']) / np.mean(df['TOTAL_READS'] / (df['END'] - df['START']))
                )
        else:
            for sample in rc.SAMPLE.unique():
                if sample == normal_name:
                    continue
                nreads_sample = rc[rc.SAMPLE == sample].iloc[0]['#READS']
                correction = nreads_normal / nreads_sample
                my_bb = big_bb[big_bb.SAMPLE == sample]

                # Correct the tumor reads propotionally to the total reads in corresponding samples
                big_bb.loc[big_bb.SAMPLE == sample, 'CORRECTED_READS'] = (my_bb.TOTAL_READS * correction).astype(np.int64)

                # Recompute RDR according to the corrected tumor reads
                big_bb.loc[big_bb.SAMPLE == sample, 'RD'] = (
                    big_bb.loc[big_bb.SAMPLE == sample, 'CORRECTED_READS']
                    / big_bb.loc[big_bb.SAMPLE == sample, 'NORMAL_READS']
                )

            if 'NORMAL_READS' not in big_bb:
                log('# NOTE: adding NORMAL_READS column to bb file', level='INFO')
                big_bb['NORMAL_READS'] = (big_bb.CORRECTED_READS / big_bb.RD).astype(np.uint32)
                
        big_bb.END = big_bb.END + 1
        big_bb.to_csv(outfile, index=False, sep='\t')

    exit(0)

    # Query command to extract CHROM, POS, and GT fields from the VCF
    query_command = [
        'bcftools', 'query',
        '-f', '%CHROM\t%POS\t[%GT]\n',
        phase
    ]

    process = subprocess.Popen(query_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    if process.returncode != 0:
        raise Exception(f"Error executing bcftools query: {stderr.decode('utf-8')}")

    # Parse stdout into a pandas dataframe
    data = [line.split() for line in stdout.decode('utf-8').strip().split('\n')]
    all_phase = pd.DataFrame(data, columns=['CHR', 'POS', 'GT'])

    # Split GT column into A and B columns
    all_phase[['A', 'B']] = all_phase['GT'].str.split('|', expand=True)
    all_phase.drop('GT', axis=1, inplace=True)

    log(msg=all_snps.head(), level='STEP')
    log(msg=all_phase.head(), level='STEP')

    def read_ranges_from_file(filename):
        df = pd.read_csv(filename, sep='\t', header=None, comment='#',
                        names=['CHR', 'source', 'feature', 'START', 'END', 'score', 'strand', 'frame', 'attribute'])
        ranges = df[['CHR', 'START', 'END']]
        return ranges
    # Read ranges from the haplotype file
    haplotype_blocks = read_ranges_from_file(haplotype_file)

    haplo_phases = {}
    for chrom in haplotype_blocks.CHR.unique():
        # Initialize the result dictionary
        range_dict = {}

        # Initialize pointers for all_phases and haplotype_blocks
        all_phases_index = 0
        ranges_index = 0

        # Get the number of rows in all_phases and haplotype_blocks
        num_all_phases = len(all_phase)
        num_ranges = len(haplotype_blocks)

        # Initialize a list to store indices for the current range
        current_range_indices = []

        while all_phases_index < num_all_phases and ranges_index < num_ranges:
            # Get current range
            start = int(haplotype_blocks.loc[ranges_index, 'START'])
            end = int(haplotype_blocks.loc[ranges_index, 'END'])
            
            # Get current POS
            pos = int(all_phase.loc[all_phases_index, 'POS'])
            
            if start <= pos <= end:
                # If POS is within the current range, add the index to the list
                current_range_indices.append(all_phases_index)
                all_phases_index += 1
            elif pos < start:
                # If POS is less than the start of the current range, move to the next POS
                all_phases_index += 1
            else:
                # If POS is greater than the end of the current range, store the current range indices and move to the next range
                range_dict[(start, end)] = all_phase.loc[current_range_indices]
                current_range_indices = []
                ranges_index += 1

        # Add the last range indices if it hasn't been added yet
        if current_range_indices:
            start = haplotype_blocks.loc[ranges_index, 'START']
            end = haplotype_blocks.loc[ranges_index, 'END']
            range_dict[(start, end)] = all_phase.loc[current_range_indices]

        # drop the keys from range_dict if the value has no rows
        range_dict = {k: v for k, v in range_dict.items() if len(v) > 0}
        haplo_phases[chrom] = range_dict

    exit(0)
    # Define the input files
    input_vcf = phase
    # haplotype_file = 'haplotype_file.txt'  # Filename for the table

    # Function to read the ranges from the haplotype file using pandas
    def read_ranges_from_file(filename):
        df = pd.read_csv(filename, sep='\t', header=None, comment='#', 
                        names=['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
        ranges = df[['chrom', 'start', 'end']].head(10).apply(tuple, axis=1).tolist()
        return ranges

    # Function to process the input VCF and extract the concatenated GT sequences for all ranges
    def process_vcf(input_vcf, ranges):
        # Sort ranges by chromosome and start position for efficient processing
        # ranges.sort()
        
        # Dictionary to store concatenated sequences for each range
        results = {range_: '' for range_ in ranges}
        
        # Query command to extract CHROM, POS, and GT fields from the VCF
        query_command = [
            'bcftools', 'query',
            '-f', '%CHROM\t%POS\t%GT\n',
            input_vcf
        ]
        
        process = subprocess.Popen(query_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            print(f'Error: {stderr.decode()}')
            return None
        
        # Process the output to extract the first elements of GT within the specified ranges
        current_range_index = 0
        range_count = len(ranges)
        
        for line in stdout.decode().strip().split('\n'):
            chrom, pos, gt = line.split('\t')
            pos = int(pos)
            
            # Only process GT values that are 0|1 or 1|0
            if gt not in {"0|1", "1|0"}:
                continue
            
            # Check if the current position falls within the current range
            while current_range_index < range_count and (chrom > ranges[current_range_index][0] or (chrom == ranges[current_range_index][0] and pos > ranges[current_range_index][2])):
                current_range_index += 1
            
            if current_range_index >= range_count:
                break
            
            # If the position is within the range, add the first element of GT to the results
            if ranges[current_range_index][0] == chrom and ranges[current_range_index][1] <= pos <= ranges[current_range_index][2]:
                results[ranges[current_range_index]] += gt.split('|')[0]
        
        return results

    # Read ranges from the haplotype file
    haplotype_blocks = read_ranges_from_file(haplotype_file)

    # Process the input VCF file and extract the concatenated GT sequences for all ranges
    results = process_vcf(input_vcf, haplotype_blocks)

    # Print the results
    for range_, sequence in results.items():
        chrom, start, end = range_
        log(msg=f'Range {chrom}:{start}-{end} -> {sequence}\n', level='STEP')
