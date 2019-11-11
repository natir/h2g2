# H2G2: Heterozygotie to Genome Graph

To reproduce this analysis run:

```
conda env create -f conda_env.yml
conda activate h2g2
jupyter notebook h2g2.ipynb
```

When you are in your jupyter notebook check your are in the good env Python conda env: h2g2


## python module install and usage

Python 3 is required

```
git clone git@github.com:natir/h2g2.git
cd h2g2
pip install -r requirements.txt
python -m h2g2 split -v file1.vcf file2.vcf … -r reference.fasta
```
