<div align='center'>
<strong><font face='Times New Roman'><h1>Nomenclature of gene IDs in metaRLK</h1></font></strong>
</div>
<br>

<center><font size="4" face='Times New Roman'>nomenclature of gene IDs of chromosome-level plant based on the results of Comparison with model plants</font></center>
<br>
<br>


<strong><font size="5" face='Times New Roman'>Introduction</font></strong>

<font size="4" face='Times New Roman'>In the metaRLK database (https://metarlk.biocloud.top/), based on the classical blastp program, we propose a format that can unify gene IDs of various plant genomes. In order to facilitate users to apply this format to rename other plant gene IDs , the following script is provided here.

Author: Qian Liu

Contact email: 18229048546@163.com</font>

<strong><font size="5" face='Times New Roman'>Requirements</font></strong>

```bash
pip3 install pandas re os argparse biopython itertools
```

<strong><font size="5" face='Times New Roman'>Usage</font></strong>

<font size="4" face='Times New Roman'>*Step1 Run blastp*</font>

<font size="4" face='Times New Roman'>For the Consideration of improving the conventional blastp running speed, we recommend using diamond(https://github.com/bbuchfink/diamond) as an alternative,for example</font>

```bash
diamond blastp --db model_plant_protein.dmnd -q seqfile -o all_genome_protein_matches_fmt6
```

<font size="4" face='Times New Roman'>Here model_plant_protein.dmnd represents the diamond-formatted database of the full model plant's protein(we provide diamond-formatted databases of two model plants here, Arabidopsis thaliana and Oryza sativa. Users can build their own diamond-formatted database), seqfile represents the full protein annotation file of the target genome, and all_genome_protein_matches_fmt6 is the running result.</font>

<font size="4" face='Times New Roman'>In the metaRLK, we recommend to align monocot genome protein annotation file to the Oryza sativa database, and align the other genome protein annotation file to the Arabidopsis thaliana database.</font>

<font size="4" face='Times New Roman'>*Step2 Get the relationship between gene IDs and their chromosome position*</font>

```bash
python3 get_all_genes_genome.py gtf_file
```

<font size="4" face='Times New Roman'>Please note that in gtf_file, each chromosome must be represented by a number. For example, the first chromosome is represented as 1</font>

<font size="4" face='Times New Roman'>*Step3 Nomenclature of gene IDs*</font>

```bash
python3 get_all_genome_new_id_split_diamond_result.py seqfile all_genome_protein_matches_fmt6 genes_genome uniprot_abbreviation_upper mode
```
<font size="4" face='Times New Roman'>Among the parameters, seqfile is the full protein annotation file of the target genome in the above blastp process, and all_genome_protein_matches_fmt6 is the running result of the above blastp process, genes_genom is the running result of the get_all_genes_genom.py script above, and uniprot_abbreviation_upper represents the abbreviation of the target plant. Here we recommend using uniprot databases' species abbreviation model to determine this parameter. The mode parameter has two values: A and O. When the target plant is a monocot plant, please use O. In other cases, please use A.</font>

<strong><font size="5" face='Times New Roman'>Help</font></strong>
<br>
<br>
<font size="4" face='Times New Roman'>Users can get the running method of each script through the -h or --help parameter,for example:</font>

```bash
python3 get_all_genes_genome.py -h
usage: get_all_genes_genome.py [-h] gtf_file

get location of genes on chromosomes

positional arguments:
  gtf_file    Genome gtf annotaion file

optional arguments:
  -h, --help  show this help message and exit
####################################################################################
python3 get_all_genome_new_id_split_diamond_result.py -h
usage: get_all_genome_new_id_split_diamond_result.py [-h]
                                                     seqfile
                                                     all_genome_protein_matches_fmt6
                                                     genes_genome
                                                     uniprot_abbreviation_upper
                                                     mode

Rename the gene IDs of chromosome-level plant based on the results of
Comparison with model plants

positional arguments:
  seqfile               Genome protein annotaion file
  all_genome_protein_matches_fmt6
                        The fmt6 results of Comparison with model plants
  genes_genome          The chromosome location of genes
  uniprot_abbreviation_upper
                        Species id abbreviation
  mode                  Comparison model

optional arguments:
  -h, --help            show this help message and exit
```

