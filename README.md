# ProtistOAZ
The OAZ diversity research in protist genomes

To run this code you need to clone the repo, download the [archive](link) with data and extract it into the repo root directory or download and prepare iMicrobe project files (see imicrobe_data section)

## Requirements
- Python >= 3.6  
- HMMER3  
- FastME >= 2.0  
- MSAPROBS  
### Python packages:
&nbsp; - eagle >= 0.0.1

## imicrobe_data
If you want to start the analysis from iMicrobe project source files download ... files and place them to imicrobe_data/source/ directory.  
Run source data preparation:
```
python -m imicrobe_data --biosamples-ids <...> --fna-path <...> --18S-path <...>
```
You can do this step interactively via Jupyter notebook (imicrobe_data/source/preparation.ipynb)  

## homologs
To start OAZ homologs search with default OAZ representative sequences run the command below:
```
python -m homologs
```
If you have not downloaded the archive with data OAZ representative sequences from repo_data/oaz_repr.fasta will be used
If you want to use custom OAZ representative sequences use 'oaz-repr-path' option:
```
python -m homologs --oaz-repr-path <...>
```
You can do this step interactively via Jupyter notebook (homologs/...ipynb)  

## translation
Have you downloaded the archive with data run command below:
```
python -m translation
```
If you have not or you want to use your homologs run results type:
```
python -m translation --oaz-diversity-path <...>
```
You can do this step interactively via Jupyter notebook (translation/...ipynb)  

## rna
