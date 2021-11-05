from setuptools import setup, find_packages

setup(
    name="kmer-table",
    version='0.1.0',
    description='k-mer table',
    packages=find_packages(),
    entry_points="""
      [console_scripts]
      kmer-table = kmer_table.kmer_table:main
    """,
)