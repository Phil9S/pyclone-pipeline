# pyclone-pipeline

pipeline to quickly deploy and compute clonal populations using pyclone

### Setup env
```
conda create -c conda-forge -n pyclone-vi --file requirements.txt --yes
```

### Install pyclone
```
pip install git+https://github.com/Roth-Lab/pyclone-vi.git
```

### Run pipeline

```
Rscript pyclone_processing_script.R input_list.txt
```
