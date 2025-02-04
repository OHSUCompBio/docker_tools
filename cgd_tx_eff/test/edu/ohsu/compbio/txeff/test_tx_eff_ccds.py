import unittest
from edu.ohsu.compbio.txeff.tx_eff_ccds import TxEffCcds

REFSEQ_CCDS_MAP = '../../../../../test-data/refseq_ccds_map.csv'

class TxEffCcdsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._tx_eff_ccds = TxEffCcds(REFSEQ_CCDS_MAP)

    def test__read_mappings(self):
        self.assertIsNotNone(self._tx_eff_ccds.refseq_to_ccds_map, "dict created")
    
    def test_refseq_to_ccds_map(self):
        self.assertEqual(self._tx_eff_ccds.refseq_to_ccds_map.get('NM_014246.4'), 'CCDS14076.1', 'mapping found')