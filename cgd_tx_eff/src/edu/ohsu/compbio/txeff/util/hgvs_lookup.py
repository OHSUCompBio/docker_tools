'''
Use the HGVS python library to lookup a single variant. This script allows us to test and debug 
HGVS processing of a variant without having to run a VCF through tx_eff_control.py.  

Created on Mar. 21, 2023

@author: pleyte
'''
from argparse import ArgumentParser
import argparse
import logging.config
import os

import hgvs.dataproviders.uta

from edu.ohsu.compbio.txeff.tx_eff_ccds import TxEffCcds
from edu.ohsu.compbio.txeff.tx_eff_hgvs import TxEffHgvs
from edu.ohsu.compbio.txeff.util.tfx_log_config import TfxLogConfig
from edu.ohsu.compbio.txeff.util.tx_eff_pysam import PysamTxEff
from edu.ohsu.compbio.txeff.variant_transcript import VariantTranscript

class HgvsLookup(object):
    def __init__(self):
        '''
        Constructor
        '''
        self.logger = logging.getLogger(__name__)

        if os.environ.get('HGVS_SEQREPO_DIR') is None:
            self.logger.warning("The HGVS_SEQREPO_DIR environment variable is not defined. The remote seqrepo database will be used.")
        else:
            logging.info(f"Using SeqRepo {os.environ.get('HGVS_SEQREPO_DIR')}")
                
        if os.environ.get('UTA_DB_URL') is None:
            logging.warning("The UTA_DB_URL environment variable is not defined. The remote UTA database will be used.")
        else:
            logging.info(f"Using UTA Database at {os.environ.get('UTA_DB_URL')}")

    def find_genotype(self, genotype: str, reference_fasta, refseq_ccds_file):
        '''
        Search UTA for transcripts associated with the genotype (e.g. 7-12345-C-G)
        '''
        chromosome, position, ref, alt = genotype.split('-')
        variant = VariantTranscript(chromosome, position, ref, alt)

        # Reference file   
        pysam = PysamTxEff(reference_fasta)
        
        # Load the refseq to ccds mappings
        tx_eff_ccds = TxEffCcds(refseq_ccds_file)

        with TxEffHgvs(pysam = pysam, refseq_ccds_map = tx_eff_ccds.refseq_to_ccds_map) as tx_eff_hgvs:
            transcripts = tx_eff_hgvs._lookup_hgvs_transcripts([variant])

        logging.info(f"Found {len(transcripts)} transcripts associated with {genotype}")
        
        self.print_result_transcripts(transcripts)

    def find_gene(self, gene:str):
        hdp = hgvs.dataproviders.uta.connect()
        gene = hdp.get_gene_info(gene)
        hdp.close()
        print(f'Gene: {gene}')

    def print_result_transcripts(self, transcripts: list):
        for transcript in transcripts:
            print(f'Variant: {transcript.chromosome}-{transcript.position}-{transcript.reference}-{transcript.alt}')
            print(f'  Gene: {transcript.hgnc_gene}')
            print(f'  g.: {transcript.sequence_variant}')
            print(f'  c.: {transcript.refseq_transcript}:{transcript.hgvs_c_dot}')
            print(f'  p.: {transcript.protein_transcript}:{transcript.hgvs_p_dot_one} / {transcript.hgvs_p_dot_three}')
            print(f'  Splicing: {transcript.splicing if transcript.splicing else "No"}')
            print('--------')
        
if __name__ == '__main__':
    logging.config.dictConfig(TfxLogConfig().utility_config)
    parser = ArgumentParser(description='Use the HGVS python library to lookup a single variant')
    parser.add_argument("-v", "--variant_genotype", help="Variant genotype (e.g. 7-123456-C-G)", type=str, required=False)
    parser.add_argument("-g", "--gene", help="Gene name", type=str, required=False)
    parser.add_argument('--reference_fasta',
                        help='Reference genome, in FASTA format.  The associated index is expected to be in the same directory.',
                        required=True)
    parser.add_argument('-c', '--ccds_map', 
                help='Input CSV or GFF with Annovar-to-CCDS mappings.', 
                type=argparse.FileType('r'), 
                required=True)    
    
    args = parser.parse_args()
    
    hgvs_lookup = HgvsLookup()
    
    if args.variant_genotype is not None:
        hgvs_lookup.find_genotype(args.variant_genotype, args.reference_fasta, args.ccds_map.name)
    elif args.gene is not None:
        hgvs_lookup.find_gene(args.gene)
    else:
        print("Please specify genotype or gene")
