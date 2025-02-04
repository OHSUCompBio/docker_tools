import os
import unittest

from edu.ohsu.compbio.txeff.tx_eff_ccds import TxEffCcds
from edu.ohsu.compbio.txeff.tx_eff_hgvs import TxEffHgvs
from edu.ohsu.compbio.txeff.util.tx_eff_pysam import PysamTxEff
from edu.ohsu.compbio.txeff.variant_transcript import VariantTranscript

FILE_HG37_REFERENCE_FASTA = '/opt/bioinformatics/Broad/Homo_sapiens_assembly19.fasta'
REFSEQ_CCDS_MAP = '../../../../../test-data/refseq_ccds_map.csv'

class TxEffHgvsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._pysam = PysamTxEff(FILE_HG37_REFERENCE_FASTA)        
        cls._tx_eff_ccds = TxEffCcds(REFSEQ_CCDS_MAP)
        
    def setUp(self):
        unittest.TestCase.setUp(self)
        os.environ.pop('HGVS_SEQREPO_DIR', None)
        os.environ.pop('HGVS_SEQREPO_URL', None)
        
    def test__merge_annovar_with_hgvs(self):
        '''
        Merge two lists of transcripts and get only those that are common to both lists
        '''
        annovar_a = self._create_variant('1', 1, 'A', 'T', 'aaa.1', gene='gene_a', variant_effect='stoploss', variant_type='exonic')
        annovar_b = self._create_variant('2', 2, 'A', 'T', 'bbb.1', gene='gene_b', variant_effect='nonsynonymous SNV', variant_type='intronic')

        hgvs_b = self._create_variant('2', 2, 'A', 'T', 'bbb.1', c_dot='c.1A>T', p_dot_one='p.L1P', p_dot_three='p.Leu1Pro', protein_transcript='NP_001')
        hgvs_c = self._create_variant('3', 3, 'A', 'T', 'ccc.1', c_dot='c.3A>T', p_dot_one='p.L3P', p_dot_three='p.Leu3Pro', protein_transcript='NP_003')

        merged_transcripts, unmerged_transcripts = TxEffHgvs(pysam = self._pysam, refseq_ccds_map = self._tx_eff_ccds.refseq_to_ccds_map)._merge_annovar_with_hgvs([annovar_a, annovar_b], [hgvs_b,hgvs_c])
        
        self.assertEqual(len(merged_transcripts), 1, 'Only one transcript is common to both lists')
        self.assertEqual(len(unmerged_transcripts), 2, 'Two transcripts are common to both lists')
    
    def test__get_unmatched_annovar_transcripts(self):
        '''
        Given two dictionary objects, create a list of values based on the the symmetric difference of the keys  
        '''        
        t1 = self._create_variant('1', 1, 'A', 'T', 'NM_1.1')
        t2 = self._create_variant('2', 2, 'A', 'T', 'NM_2.2')
        t3 = self._create_variant('3', 3, 'A', 'T', 'NM_3.3')
        t4 = self._create_variant('4', 4, 'A', 'T', 'NM_4.4')
        t5 = self._create_variant('5', 5, 'A', 'T', 'NM_5.5')
        t6 = self._create_variant('6', 6, 'A', 'T', 'NM_6.6')
        
        # Dict a contains 1,2,3,4    
        a = {t1.get_label():t1, t2.get_label():t2, t3.get_label():t3, t4.get_label():t4}
        
        # Dict b contains 3,4,5,6
        b = {t3.get_label():t3, t4.get_label():t4, t5.get_label():t5, t6.get_label():t6}
        
        # Unmatched will be 1,2,5,6
        c = TxEffHgvs(pysam = self._pysam, refseq_ccds_map = self._tx_eff_ccds.refseq_to_ccds_map)._get_unmatched_annovar_transcripts(a,b)
        
        self.assertEqual(len(c), 4, 'Symmetric difference') 
        
        # Use genotype labels rather than try to compare objects
        labels = [x.get_label() for x in c]
        self.assertTrue(t1.get_label() in labels, 't1 is unmatched')
        self.assertTrue(t2.get_label() in labels, 't2 is unmatched')

        self.assertTrue(t3.get_label() not in labels, 't3 not in list')
        self.assertTrue(t4.get_label() not in labels, 't4 not in list')
        
        self.assertTrue(t5.get_label() in labels, 't5 is unmatched')
        self.assertTrue(t6.get_label() in labels, 't6 is unmatched')
        
    def test__get_the_best_transcripts(self):
        '''
        Given three versions of a transcript return then one with the most values and the latest version
        '''
        merged_transcripts = [self._create_variant('1', 1, 'A', 'T', 'NM_000.1', gene='AAA', p_dot_one='p.L1P', p_dot_three='p.Leu1Pro', protein_transcript='NP_123.1', variant_effect='X'),
                              self._create_variant('1', 1, 'A', 'T', 'NM_000.3', gene='AAA'),
                              self._create_variant('1', 1, 'A', 'T', 'NM_000.2', gene='AAA', p_dot_one='p.L1P', p_dot_three='p.Leu1Pro', protein_transcript='NP_123.1', variant_effect='Y')]

        best = TxEffHgvs(pysam = self._pysam, refseq_ccds_map = self._tx_eff_ccds.refseq_to_ccds_map)._get_the_best_transcripts(merged_transcripts)
        
        self.assertEqual(len(best), 1, 'Only one transcript expected')
        self.assertEqual(best[0].get_label(), '1-1-A-T-NM_000.2', 'Most complete, lastest version')
    
    def test__get_the_best_transcripts_onlyOnePdot_zeroPoints(self):
        vt = VariantTranscript('1', 123, 'C', 'G')
        vt.hgvs_p_dot_one = 'p.E447del'
        vt.protein_transcript = 'NP_123.1'
        vt.hgvs_amino_acid_position = 447
        vt.hgvs_p_dot_three = ''        
        self.assertEqual(vt._get_self_score(), 0, "Score without p3")
        
    def test___get_self_score_allProteinFields_onePoint(self):
        vt = VariantTranscript('1', 123, 'C', 'G')
        vt.hgvs_p_dot_one = 'p.E447del'
        vt.protein_transcript = 'NP_123.1'
        vt.hgvs_p_dot_three = 'p.Glu447del'
        vt.hgvs_amino_acid_position = 447
        self.assertEqual(vt._get_self_score(), 1, "Score with all protein fields")
        
    def test___get_self_score_pDotNoProteinTranscript_zeroPoints(self):
        vt = VariantTranscript('1', 123, 'C', 'G')
        vt.hgvs_p_dot_one = 'p.E447del'
        vt.hgvs_p_dot_three = 'p.Glu447del'
        vt.hgvs_amino_acid_position = 447
        vt.protein_transcript = ''                
        self.assertEqual(vt._get_self_score(), 0, "Score without protein transcript")
        
    def test__correct_indel_coords(self):
        '''
        Given a genotype the _correct_indel_coords function returns the nomenclature describing the variant change (the g., minus the g. prefix) 
        ''' 
        # Substitution 
        genotype = '1-123-G-A'
        chromosome, position, ref, alt = genotype.split('-')
        pos_part = TxEffHgvs(pysam = self._pysam, refseq_ccds_map = self._tx_eff_ccds.refseq_to_ccds_map)._correct_indel_coords(chromosome, int(position), ref, alt)
        self.assertEqual(pos_part, '123G>A', "g. is incorrect")

        # Insertion 
        genotype = '1-123-G-AC'
        chromosome, position, ref, alt = genotype.split('-')
        pos_part = TxEffHgvs(pysam = self._pysam, refseq_ccds_map = self._tx_eff_ccds.refseq_to_ccds_map)._correct_indel_coords(chromosome, int(position), ref, alt)
        self.assertEqual(pos_part, '123_124insC', "g. is incorrect")

        # Multi-nucleotide insertion
        genotype = '7-33054117-T-TAAGA'
        chromosome, position, ref, alt = genotype.split('-')
        pos_part = TxEffHgvs(pysam = self._pysam, refseq_ccds_map = self._tx_eff_ccds.refseq_to_ccds_map)._correct_indel_coords(chromosome, int(position), ref, alt)
        self.assertEqual(pos_part, '33054117_33054118insAAGA', "g. is incorrect")

        # Deletion
        genotype = '1-123-AC-G'
        chromosome, position, ref, alt = genotype.split('-')
        pos_part = TxEffHgvs(pysam = self._pysam, refseq_ccds_map = self._tx_eff_ccds.refseq_to_ccds_map)._correct_indel_coords(chromosome, int(position), ref, alt)
        self.assertEqual(pos_part, '124del', "g. is incorrect")

        # Multi-nucleotide deletion
        genotype = '9-2039776-ACAG-A'
        chromosome, position, ref, alt = genotype.split('-')
        pos_part = TxEffHgvs(pysam = self._pysam, refseq_ccds_map = self._tx_eff_ccds.refseq_to_ccds_map)._correct_indel_coords(chromosome, int(position), ref, alt)
        self.assertEqual(pos_part, '2039777_2039779del', "g. is incorrect")

        # Duplication
        genotype = '3-142274739-A-AT'
        chromosome, position, ref, alt = genotype.split('-')
        pos_part = TxEffHgvs(pysam = self._pysam, refseq_ccds_map = self._tx_eff_ccds.refseq_to_ccds_map)._correct_indel_coords(chromosome, int(position), ref, alt)
        self.assertEqual(pos_part, '142274749dup', "g. is incorrect")

        # Multi-nucleotide duplication
        genotype = '5-79950724-G-GCCGCAGCGC'
        chromosome, position, ref, alt = genotype.split('-')
        pos_part = TxEffHgvs(pysam = self._pysam, refseq_ccds_map = self._tx_eff_ccds.refseq_to_ccds_map)._correct_indel_coords(chromosome, int(position), ref, alt)
        self.assertEqual(pos_part, '79950725_79950733dup', "g. is incorrect")

        # Indel
        genotype = '5-112175461-CAGTTCACTTGA-CGTC'
        chromosome, position, ref, alt = genotype.split('-')
        pos_part = TxEffHgvs(pysam = self._pysam, refseq_ccds_map = self._tx_eff_ccds.refseq_to_ccds_map)._correct_indel_coords(chromosome, int(position), ref, alt)
        self.assertEqual(pos_part, '112175462_112175472delinsGTC', "g. is incorrect")
        
        # Repeat
        genotype = '1-153924029-AG-A'
        chromosome, position, ref, alt = genotype.split('-')
        pos_part = TxEffHgvs(pysam = self._pysam, refseq_ccds_map = self._tx_eff_ccds.refseq_to_ccds_map)._correct_indel_coords(chromosome, int(position), ref, alt)
        self.assertEqual(pos_part, '153924033del', "g. is incorrect")
        
    def _create_variant(self, chromosome, pos, ref, alt, transcript, gene=None, c_dot=None, p_dot_one=None, p_dot_three=None, variant_effect=None, variant_type=None, protein_transcript=None):
        '''
        Helper method for creating variant-transcript objects
        ''' 
        variant_transcript = VariantTranscript(chromosome, pos, ref, alt)
        variant_transcript.refseq_transcript = transcript
        variant_transcript.hgnc_gene = gene
        variant_transcript.hgvs_c_dot = c_dot
        variant_transcript.hgvs_p_dot_one = p_dot_one
        variant_transcript.hgvs_p_dot_three = p_dot_three
        variant_transcript.variant_effect = variant_effect
        variant_transcript.variant_type = variant_type
        variant_transcript.protein_transcript = protein_transcript
        
        return variant_transcript

    def test___configure_sequence_source_noParameters_useNcbi(self):
        sequence_source = None
        tx_eff_hgvs = TxEffHgvs(pysam = self._pysam, refseq_ccds_map = self._tx_eff_ccds.refseq_to_ccds_map, sequence_source = sequence_source)        
        self.assertEqual(tx_eff_hgvs._configure_sequence_source(), 1, "Return value should indicate NCBI will be used")
    
    def test___configure_sequence_source_requestNcbi_useNcbi(self):
        sequence_source = 'ncbi'
        tx_eff_hgvs = TxEffHgvs(pysam = self._pysam, refseq_ccds_map = self._tx_eff_ccds.refseq_to_ccds_map, sequence_source = sequence_source)        
        self.assertEqual(tx_eff_hgvs._configure_sequence_source(), 1, "Return value should indicate NCBI will be used")
        
    def test___configure_sequence_source_invalid_raiseError(self):
        sequence_source = 'nonsense'
        tx_eff_hgvs = TxEffHgvs(pysam = self._pysam, refseq_ccds_map = self._tx_eff_ccds.refseq_to_ccds_map, sequence_source = sequence_source)
        self.assertRaises(ValueError,tx_eff_hgvs._configure_sequence_source)
        
    def test___configure_sequence_source_envUrl_urlConfigured(self):
        # URL overrides directory
        tx_eff_hgvs = TxEffHgvs(pysam = self._pysam, refseq_ccds_map = self._tx_eff_ccds.refseq_to_ccds_map)
        os.environ['HGVS_SEQREPO_URL'] = "http://localhost:5000/seqrepo"
        os.environ['HGVS_SEQREPO_DIR'] = "/path"
        self.assertEqual(tx_eff_hgvs._configure_sequence_source(), 2, "Return value should indicate URL will be used")
        self.assertIsNone(os.environ.get('HGVS_SEQREPO_DIR'), "HGVS_SEQREPO_DIR should be unset")
        
    def test___configure_sequence_source_envDirectory_directoryConfigured(self):
        tx_eff_hgvs = TxEffHgvs(pysam = self._pysam, refseq_ccds_map = self._tx_eff_ccds.refseq_to_ccds_map)
        os.environ['HGVS_SEQREPO_DIR'] = "/path"
        self.assertEqual(tx_eff_hgvs._configure_sequence_source(), 3, "Return value should indicate URL will be used")
    
    def test___configure_sequence_source_argUrlAndEnvUrl_argUrlConfigured(self):
        # The argument (tx_eff_hgvs._sequence_source) overrides any environment variable
        url = "http://host1/seqrepo"
        tx_eff_hgvs = TxEffHgvs(pysam = self._pysam, refseq_ccds_map = self._tx_eff_ccds.refseq_to_ccds_map, sequence_source = url)
        
        os.environ['HGVS_SEQREPO_URL'] = "http://host2/seqrepo"
        self.assertEqual(tx_eff_hgvs._configure_sequence_source(), 4, "Return value should indicate URL will be used")
        self.assertEqual(os.environ.get('HGVS_SEQREPO_URL'), url, "HGVS_SEQREPO_URL should be set to the argument path")
    
    def test___configure_sequence_source_argDirectory_argDirectoryConfigured(self):
        tx_eff_hgvs = TxEffHgvs(pysam = self._pysam, refseq_ccds_map = self._tx_eff_ccds.refseq_to_ccds_map, sequence_source = "/path0")
        
        os.environ['HGVS_SEQREPO_URL'] = "http://host/seqrepo"
        self.assertEqual(tx_eff_hgvs._configure_sequence_source(), 5, "Return value should indicate file repository will be used")
        self.assertIsNone(os.environ.get('HGVS_SEQREPO_URL'), "HGVS_SEQREPO_URL should be unset")
            
    def test___configure_sequence_source_argDirectoryAndEnvDirectory_argDirectoryConfigured(self):
        # The argument (tx_eff_hgvs._sequence_source) overrides any environment variable
        path = "/path0"
        tx_eff_hgvs = TxEffHgvs(pysam = self._pysam, refseq_ccds_map = self._tx_eff_ccds, sequence_source = path)
        os.environ['HGVS_SEQREPO_DIR'] = "/path1"
        os.environ['HGVS_SEQREPO_URL'] = "http://host/seqrepo"
        self.assertEqual(tx_eff_hgvs._configure_sequence_source(), 5, "Return value should indicate file repository will be used")
        self.assertEqual(os.environ.get('HGVS_SEQREPO_DIR'), path, "HGVS_SEQREPO_DIR should be set to the argument path")
        self.assertIsNone(os.environ.get('HGVS_SEQREPO_URL'), "HGVS_SEQREPO_URL should be unset")
