from setuptools import setup, find_packages
from codecs import open
from os import path

this_folder = path.abspath(path.dirname(__file__))
with open(path.join(this_folder,'README.md'),encoding='utf-8') as inf:
  long_description = inf.read()

setup(
  name='AlignQC',
  version='1.3.0',
  description='Python tools for working with biological sequence data',
  long_description=long_description,
  url='https://github.com/jason-weirather/AlignQC',
  author='Jason L Weirather',
  author_email='jason.weirather@gmail.com',
  license='Apache License, Version 2.0',
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'License :: OSI Approved :: Apache Software License'
  ],
  keywords='bioinformatics, sequence, alignment',
  packages=['alignqc'],
  entry_points = {
    'console_scripts':['alignqc=alignqc.alignqc:main']
  }
)
  
