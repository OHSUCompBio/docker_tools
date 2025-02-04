'''
Created on Jan 28, 2025

@author: pleyte
'''
import unittest
from edu.ohsu.compbio.txeff.variant_transcript import VariantTranscript

def create_variant(chromosome, pos, ref, alt, transcript, gene=None, c_dot=None, p_dot_one=None, p_dot_three=None, variant_effect=None, variant_type=None, protein_transcript=None):
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
    
class VariantTranscriptTest(unittest.TestCase):

    def test___lt___sameGenotypes_notLessThanEach(self):
        vt0 = create_variant('1', 1, 'A', 'T', 'NM_123.1')
        vt1 = create_variant('1', 1, 'A', 'T', 'NM_123.1')
        
        self.assertFalse(vt0 < vt1, 'variant one is not less than variant two')
        
    def test___lt___differentGenotypes_compareAsString(self):
        vt0 = create_variant('1', 1, 'A', 'C', 'NM_111.1', gene = 'GENE', c_dot = 'c.38G>A')
        vt1 = create_variant('X', 100, 'T', 'G', 'NM_333.3')
        
        # vt0 has more fields filled in so it's score is greater than vt1
        self.assertGreater(vt0._get_self_score(), vt1._get_self_score(), "vt0's score > vt1")
        
        # But vt0 and vt1 have different genotypes so instead of comparing scores, they are compared as strings 
        self.assertLess(vt0, vt1, "vt0 as a string is < vt1")

    def test___lt___sameGenotypesSameScoreSameAccession_compareVersion(self):
        vt0 = create_variant('1', 1, 'A', 'C', 'NM_111.1')
        vt1 = create_variant('1', 1, 'A', 'C', 'NM_111.2')
        
        self.assertEquals(vt0._get_self_score(), vt1._get_self_score(), "equal scores")
        
        # When accessions are the same the higher version ranks higher  
        self.assertLess(vt0, vt1, "accession version .1 is less than accession version .2")

    def test___lt___sameGenotypesSameScoreDifferentAccession_compareAccession(self):
        vt0 = create_variant('1', 1, 'A', 'C', 'NM_111.1')
        vt1 = create_variant('1', 1, 'A', 'C', 'NM_222.2')        
        
        self.assertEquals(vt0._get_self_score(), vt1._get_self_score(), "equal scores")
        
        # When accessions are different the lower accession ranks higher because it is believed to have been curated for longer.  
        self.assertGreater(vt0, vt1, "accession 1 greater than than accession 2")

    def test___lt___sameVariantDifferentAccessionType_compareAsString(self):
        vt0 = create_variant('1', 1, 'A', 'C', 'NM_111.1')
        vt1 = create_variant('1', 1, 'A', 'C', 'CCDS222.2')
         
        # because "NM_" > "CCDS"
        self.assertGreater(vt0, vt1, "vt0 as string is > vt1")
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()