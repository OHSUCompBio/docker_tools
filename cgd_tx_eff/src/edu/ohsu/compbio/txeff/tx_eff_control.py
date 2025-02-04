'''
Created on May 18, 2022

@author: pleyte
'''
import argparse
import datetime
import logging.config
import sys
import time

import biocommons.seqrepo
import hgvs

from edu.ohsu.compbio.txeff.tx_eff_annotate import TxEffAnnotate
from edu.ohsu.compbio.txeff.tx_eff_annovar import TxEffAnnovar
from edu.ohsu.compbio.txeff.tx_eff_ccds import TxEffCcds
from edu.ohsu.compbio.txeff.tx_eff_hgvs import TxEffHgvs
from edu.ohsu.compbio.txeff.tx_eff_vcf import TxEffVcf
from edu.ohsu.compbio.txeff.util.tfx_log_config import TfxLogConfig
from edu.ohsu.compbio.txeff.util.tx_eff_pysam import PysamTxEff

VERSION = '0.7.8'

def _parse_args():
    '''
    Validate and return command line arguments.
    '''
    parser = argparse.ArgumentParser(description='Read variants from a VCF, add Annovar variant effects, correct nomenclature using HGVS, and write out a new VCF.')

    parser.add_argument('-i', '--in_vcf', 
                help='Input VCF', 
                type=argparse.FileType('r'), 
                required=True)

    parser.add_argument('-c', '--ccds_map', 
                help='Input CSV or GFF with Annovar-to-CCDS mappings.', 
                type=argparse.FileType('r'), 
                required=True)
    
    parser.add_argument('--annovar_variant_function', 
                        help='Annovar variant_function file',
                        type=argparse.FileType('r'),
                        required=True)

    parser.add_argument('--annovar_exonic_variant_function', 
                        help='Annovar exonic_variant_function file',
                        type=argparse.FileType('r'),
                        required=True)

    parser.add_argument('--reference_fasta',
                        help='Reference genome, in FASTA format.  The associated index is expected to be in the same directory.',
                        required=True)

    parser.add_argument('-o', '--out_vcf',
                    help='Output VCF', 
                    type=argparse.FileType('w'), 
                    required=True)
    
    parser.add_argument('-b', '--benchmark',
                    help='Write benchmarking information to file', 
                    action='store_true')
    
    parser.add_argument('-s', '--sequence_source',
                    help='Method for looking up reference sequences (can be a url, path, or "ncbi")')
    
    parser.add_argument('-t', '--threads', default = 1, type=int,
                    help="Number of threads to use", )
    
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    
    args = parser.parse_args()
    
    if args.sequence_source \
       and not (args.sequence_source.lower().startswith('http') or args.sequence_source.startswith("/") or args.sequence_source == 'ncbi'):
        raise ValueError(f'--sequence_source is expected to be a url, directory, or "ncbi": "{args.sequence_source}"')

    if args.threads > 1 and args.benchmark:
        raise ValueError('Benchmarking is not possible when multiple threads are used')
    
    return args

def _get_time(real_seconds, user_seconds):
    """
    Take the real time and user time in seconds and return a string that breaks the time into hours, minutes, and seconds. 
    """
    
    dt_real = datetime.timedelta(seconds = real_seconds)
    dt_user = datetime.timedelta(seconds = user_seconds)
    
    return str(dt_real), str(dt_user)
    
def _main():
    '''
    main function
    '''
    print(f"tfx_cgd {VERSION} is starting...")
    print("Python version: " + sys.version.replace('\n', ''))
    print(f"biocommons.seqrepo version: {biocommons.seqrepo.__version__}")
    print(f"hgvs version: {hgvs.__version__}")
    
    logging.config.dictConfig(TfxLogConfig().log_config)
    logger = logging.getLogger("edu.ohsu.compbio.txeff.tx_eff_control")

    if logging.root.handlers[0].stream: 
        output = str(logging.root.handlers[0].stream.name)
    else:
        output = logging.root.handlers[0].baseFilename
        
    print(f"Log level={logging.root.getEffectiveLevel()}, output={output}")

    start_time_real = time.perf_counter()
    start_time_user = time.process_time()

    args = _parse_args()
    
    print(f"Thread count: {args.threads}")
    
    # Use tx_eff_annovar to read annovar records
    annovar_records = TxEffAnnovar().get_annovar_records(args.annovar_variant_function.name, args.annovar_exonic_variant_function.name)

    # Load the reference genome with pysam.
    pysam = PysamTxEff(args.reference_fasta)

    # Load the refseq to ccds mappings
    tx_eff_ccds = TxEffCcds(args.ccds_map.name)
    
    # Look for additional transcripts in the HGVA/UTA database and merge them with the annovar records.
    with TxEffHgvs(pysam = pysam, 
                   refseq_ccds_map = tx_eff_ccds.refseq_to_ccds_map,
                   sequence_source = args.sequence_source,
                   threads = args.threads,
                   benchmark = args.benchmark) as tx_eff_hgvs:
        merged_transcripts = tx_eff_hgvs.get_updated_hgvs_transcripts(annovar_records)
    
    # Close the reference FASTA
    pysam.my_fasta.close()

    # Add additional annotations to each variant
    tx_eff_annotate = TxEffAnnotate(args.sequence_source)
    tx_eff_annotate.annotate(merged_transcripts)
    
    # Use tx_eff_vcf to write the transcript effects to a VCF
    TxEffVcf(VERSION, args.in_vcf.name, args.out_vcf.name).create_vcf(merged_transcripts)
    
    print(f"Wrote {len(merged_transcripts)} transcripts to {args.out_vcf.name}")

    stop_time_real = time.perf_counter()
    stop_time_user = time.process_time()
    real_time, user_time = _get_time(stop_time_real - start_time_real, stop_time_user - start_time_user) 
    
    print(f"Time: real {real_time}, user {user_time}")
    logger.info(f"Time: real {real_time}, user {user_time}")

if __name__ == '__main__':
    _main()