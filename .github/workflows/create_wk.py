TEMPLATE="""name: __NAME__ Testing

on:
  workflow_dispatch:
  push:
    branches:
      - main
      - dev
      - refactoring
    paths:
      - bioconvert/__NAME__.py
      - .github/workflows/__NAME__.yml
  pull_request:
    branches-ignore: []
    paths:
      - bioconvert/__NAME__.py
  schedule:
    - cron: '0 0 13 * *'

jobs:
  build-linux:
    runs-on: ubuntu-latest
    strategy:
      max-parallel: 5
      matrix:
        python: [3.8, 3.9, '3.10']
      fail-fast: false

    steps:

    - name: install graphviz and curl
      run: |
        sudo apt-get update
        sudo apt-get install -y graphviz-dev
        sudo apt-get install libcurl4-gnutls-dev
    - name: checkout git repo
      uses: actions/checkout@v2

    - name: conda/mamba
      uses: mamba-org/provision-with-micromamba@main
      with:
        cache-downloads: true
        environment-file: false
        environment-name: installation
        channels:
          conda-forge,bioconda,defaults,r
        extra-specs: |
          python=${{ matrix.python }}
          easydev
          biosniff
          colorlog
          deeptools
          gffread
          pandas
          biopython>=1.70
          mappy
          matplotlib-base
          networkx
          pyyaml
          pysam
          pyexcel
          pyexcel-ods3
          pyexcel-xls
          pyexcel-xlsx
          pyBigWig
          py2bit
          statsmodels
          tqdm
          bamtools
          bcftools
          bedtools
          bedops
          dsrc
          go==1.10.3
          goalign
          gotree
          mosdepth
          pbzip2
          pigz
          plink
          sambamba
          samtools>=1.9
          seqtk
          seqkit
          squizz
          sra-tools
          ucsc-wigtobigwig
          ucsc-twobittofa
          ucsc-fatotwobit
          ucsc-bedgraphtobigwig
          ucsc-bigwigtobedgraph
          wiggletools
          sed
          mawk
          xlrd

    #    mamba install -c conda-forge -c bioconda --quiet -y samtools bedtools bamtools mosdepth pbzip2 pigz dsrc sambamba squizz

    - name: Install with pip
      shell: bash -l {0}
      run: |
        pip install .[testing]
    - name: testing
      shell: bash -l {0}
      run: |
        pytest -n 1  --cov-report term --cov=bioconvert.__NAME__ test/test___NAME__.py
"""


for x in ['abi2fasta','abi2fastq','abi2qual','bam2bedgraph','bam2bigwig','bam2cov','bam2cram','bam2fasta','bam2fastq','bam2json','bam2sam','bam2tsv','bam2wiggle','bcf2vcf','bcf2wiggle','bed2wiggle','bedgraph2bigwig','bedgraph2cov','bedgraph2wiggle','bigbed2bed','bigbed2wiggle','bigwig2bedgraph','bigwig2wiggle','bplink2plink','bplink2vcf','bz22gz','clustal2fasta','clustal2nexus','clustal2phylip','clustal2stockholm','cram2bam','cram2fasta','cram2fastq','cram2sam','csv2tsv','csv2xls','dsrc2gz','embl2fasta','embl2genbank','fasta2clustal','fasta2faa','fasta2fasta_agp','fasta2fastq','fasta2genbank','fasta2nexus','fasta2phylip','fasta2twobit','fasta_qual2fastq','fastq2fasta_qual','fastq2fasta','fastq2qual','genbank2embl','genbank2fasta','genbank2gff3','gfa2fasta','gff22gff3','gff32gff2','gff32gtf', 'gz2bz2','gz2dsrc','json2yaml','maf2sam','newick2nexus','newick2phyloxml','nexus2clustal','nexus2fasta','nexus2newick','nexus2phylip','nexus2phyloxml','ods2csv','pdb2faa', 'phylip2clustal','phylip2fasta','phylip2nexus','phylip2stockholm','phylip2xmfa','phyloxml2newick','phyloxml2nexus','plink2bplink','plink2vcf','sam2bam','sam2cram','sam2paf','scf2fasta','scf2fastq','sra2fastq','stockholm2clustal','stockholm2phylip','tsv2csv','twobit2fasta','vcf2bcf','vcf2bed','vcf2bplink','vcf2plink','vcf2wiggle','wig2bed','xls2csv','xlsx2csv','xmfa2phylip','yaml2json']:
    with open(f"{x}.yml", "w") as fout:
        fout.write(TEMPLATE.replace("__NAME__", x))

