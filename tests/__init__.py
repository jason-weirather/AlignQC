import unittest
import os, sys, shutil
from tempfile import mkdtemp, gettempdir
from alignqc import alignqc
from glob import glob

def _initialize_example_data(self):
        self.current = os.path.dirname(os.path.abspath(__file__))
        self.bam = os.path.join(self.current,'..','example_data','chr21chr22chrM.bam')
        self.fa = os.path.join(self.current,'..','example_data','chr21chr22chrM.fa.gz')
        self.gtf = os.path.join(self.current,'..','example_data','chr21chr22chrM.gencode25.gtf.gz')
        self.gpd = os.path.join(self.current,'..','example_data','chr21chr22chrM.gencode25.gpd.gz')
        return

class TestAnalysisMinimal(unittest.TestCase):
    """ Test a minimal analysis run """
    @classmethod
    def setUpClass(self):
        """Get the temporary directors and run the command"""
        _initialize_example_data(self)
        self.dirpath = mkdtemp(prefix="weirathe.",dir=gettempdir())
        self.outdir = os.path.join(self.dirpath,'output')
        cmd = ['alignqc','analyze',
               self.bam,
               '--no_genome','--no_transcriptome',
               '--output_folder',self.outdir]
        alignqc.external_cmd(cmd)
        return
    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.dirpath)
        return
    def test_PNG(self):
        """See if png figures got made"""
        png_files = glob(os.path.join(self.outdir,'plots')+'/*.png')
        self.assertEqual(len(png_files),5)
        """Test that the right amount of png images were made"""
        for file in png_files:
           self.assertTrue(os.stat(file).st_size > 8000)
           """Test image is not too small and that its there"""
    def test_PDF(self):
        """See if pdf figures got made"""
        pdf_files = glob(os.path.join(self.outdir,'plots')+'/*.pdf')
        self.assertEqual(len(pdf_files),5)
        """Test that the right amount of png images were made"""
        for file in pdf_files:
           self.assertTrue(os.stat(file).st_size > 2000)
           """Test image is not too small and that its there"""

class TestAnalysisComplete(unittest.TestCase):
    """ Test a minimal analysis run """
    @classmethod
    def setUpClass(self):
        """Get the temporary directors and run the command"""
        _initialize_example_data(self)
        self.dirpath = mkdtemp(prefix="weirathe.",dir=gettempdir())
        self.outdir = os.path.join(self.dirpath,'output')
        self.outxhtml = os.path.join(self.dirpath,'out.xhtml')
        self.outportable = os.path.join(self.dirpath,'portable.xhtml')
        cmd = ['alignqc','analyze',
               self.bam,
               '-g',self.fa,
               '-t',self.gtf,
               '--context_error_stopping_point','900',
               '--alignment_error_max_length','100000',
               '--output_folder',self.outdir,
               '-o',self.outxhtml,
               '--portable_output',self.outportable]
        alignqc.external_cmd(cmd)
        return
    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.dirpath)
    def test_PNG(self):
        """See if png figures got made"""
        png_files = glob(os.path.join(self.outdir,'plots')+'/*.png')
        self.assertEqual(len(png_files),13)
        """Test that the right amount of png images were made"""
        for file in png_files:
           self.assertTrue(os.stat(file).st_size > 8000)
           """Test image is not too small and that its there"""
    def test_PDF(self):
        """See if pdf figures got made"""
        pdf_files = glob(os.path.join(self.outdir,'plots')+'/*.pdf')
        self.assertEqual(len(pdf_files),13)
        """Test that the right amount of png images were made"""
        for file in pdf_files:
           self.assertTrue(os.stat(file).st_size > 3000)
           """Test image is not too small and that its there"""
    def test_dump(self):
        """Test the dump command"""
        ofile = os.path.join(self.outdir,'best.bed')
        cmd = ['alignqc','dump',self.outxhtml,'-e','best.sorted.bed','-o',ofile]
        alignqc.external_cmd(cmd)
        cnt = 0
        with open(ofile) as inf:
            for line in inf: cnt += 1
        self.assertTrue(cnt==1083)
    def test_compare(self):
        """Test the compare command"""
        odir = os.path.join(self.outdir,'compare')
        cmd = ['alignqc','compare',self.outxhtml,self.outxhtml,self.outxhtml,'-o',odir]
        ofile = os.path.join(odir,'stats_table.txt')
        alignqc.external_cmd(cmd)
        cnt = 0
        with open(ofile) as inf:
            for line in inf: cnt += 1
        self.assertTrue(cnt==4)
if __name__ == '__main__':
    unittest.main()