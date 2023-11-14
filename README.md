# pyclone-pipeline

pipeline to quickly deploy and compute clonal populations using [pyclone](https://github.com/Roth-Lab/pyclone-vi).

### Setup env
```
conda create -c conda-forge -n pyclone-vi --file requirements.txt --yes
```

### Activate env
```
conda activate pyclone-vi
```

### Install pyclone
```
pip install git+https://github.com/Roth-Lab/pyclone-vi.git
```

### Run pipeline

```
Rscript pyclone_processing_script.R input_list.txt [output_directory]
```

## Parameters

 - For related samples, sample ids with matching patient ids are processed together.
 - If all samples are independent, sample_id should equal patient_id and be unique.
 - Required files should be provided as absolute file paths.
 - Purity should be specified as a float between 0.0 - 1.0.

Output directory by default is the current working directory. A different output directory can be specified as the second argument.
Additional parameters can be changed in the main script.
