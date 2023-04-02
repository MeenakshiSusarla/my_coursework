import unittest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

class TestIdentifyClosestBreed(unittest.TestCase):

    def test_identify_closest_breed(self):
    
        # Define the input files for the test
        database_file = 'test_database.fa'
        query_file = 'test_query.fa'
    
        # Create test sequences
        database_seq1 = Seq('AGCTAGCTAGCT')
        database_seq2 = Seq('AGCTAGCTAGCT')
        database_seq3 = Seq('AGCTAGGTAGCT')
        query_seq = Seq('AGCTAGGTAGTT')
    
        # Create test database file
        SeqIO.write([database_seq1, database_seq2, database_seq3], database_file, 'fasta')
    
        # Create test query file
        SeqIO.write(query_seq, query_file, 'fasta')
    
        # Call the function to identify the closest breed
        closest_breed, closest_distance = identify_closest_breed(database_file, query_file)
    
        # Check the expected output
        self.assertEqual(closest_breed, 'SeqRecord1')
        self.assertAlmostEqual(closest_distance, 2.0)
if name == 'main':
unittest.main()


import unittest
from Bio import AlignIO
from Bio.Align import PairwiseAligner
import os

class TestClosestBreed(unittest.TestCase):

    def setUp(self):
        # Create test alignment file
        with open('test_alignment.fas', 'w') as f:
            f.write('>Seq1\nATCG\n>Seq2\nACGT\n>Seq3\nAGCT\n')
        
    def tearDown(self):
        # Remove test alignment file and output file
        os.remove('test_alignment.fas')
        os.remove('output.txt')
    
    def test_closest_breed(self):
    
        # Call the function to perform pairwise sequence alignment and count differences
        closest_breed('test_alignment.fas')
    
        # Check if the output file is created
        self.assertTrue(os.path.isfile('output.txt'))
    
        # Check the contents of the output file
        with open('output.txt', 'r') as f:
            output = f.read()
            self.assertIn('Alignment between Seq1 and Seq2', output)
            self.assertIn('>Seq1\nATCG\n', output)
            self.assertIn('>Seq2\nACGT\n', output)
            self.assertIn('Number of differences: 2', output)
            self.assertIn('Alignment between Seq1 and Seq3', output)
            self.assertIn('>Seq1\nATCG\n', output)
            self.assertIn('>Seq3\nAGCT\n', output)
            self.assertIn('Number of differences: 2', output)
            self.assertIn('Alignment between Seq2 and Seq3', output)
            self.assertIn('>Seq2\nACGT\n', output)
            self.assertIn('>Seq3\nAGCT\n', output)
            self.assertIn('Number of differences: 3', output)

if name == 'main':
unittest.main()


import unittest
import os
import matplotlib.pyplot as plt
from scipy.stats import shapiro

class TestCalculatePValue(unittest.TestCase):

    def setUp(self):
        # Create test alignment file
        with open('test_alignment.fas', 'w') as f:
            f.write('>Seq1\nATCG\n>Seq2\nACGT\n>Seq3\nAGCT\n')
        
    def tearDown(self):
        # Remove test alignment file and output file
        os.remove('test_alignment.fas')
        os.remove('output.txt')
    
    def test_calculate_pvalue(self):
    
        # Call the function to perform Shapiro-Wilk test and plot histogram
        calculate_pvalue('test_alignment.fas')
    
        # Check if the plot is displayed
        self.assertTrue(plt.fignum_exists(1))
    
        # Check if the returned p-value is correct
        self.assertAlmostEqual(shapiro([4, 3, 2, 5, 4, 3, 2, 1])[1], 0.0020490046565088034, places=6)

if name == 'main':
unittest.main()

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from tempfile import NamedTemporaryFile
from my_script import phylo_tree

# Create a temporary file containing some DNA sequences
with NamedTemporaryFile(mode='w', delete=False) as tmp_file:
    seqs = [SeqRecord(Seq('ATCGATCG', generic_dna), id='seq1'),
            SeqRecord(Seq('GCTAGCTA', generic_dna), id='seq2'),
            SeqRecord(Seq('CCCCCCCC', generic_dna), id='seq3'),
            SeqRecord(Seq('GGGGGGGG', generic_dna), id='seq4')]
    SeqIO.write(seqs, tmp_file.name, 'fasta')

    # Call the function and assert that it produces a tree
    phylo_tree(tmp_file.name)