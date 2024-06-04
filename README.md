magicBatch

## Overview

This is a modified version of the MAGIC algorithm that allows for batch-corrected imputation of multi-batch datasets.

## Installation

1) install the magicBatch python package using one of two methods. 

	From within R (preferred):

	```{r, eval=FALSE, message=FALSE, warning=FALSE, results = 'hide'}
	python_path <- system("which python3", intern = TRUE)
	system(paste(python_path, "-m pip install magicBatch"))
	```

	From the command line:

	```{}
	pip install magicBatch
	```
	And copy the output of

	```{}
	which python
	```
	as this will be needed to invoke the correct python runtime from within R.

2) install the magicBatch R package:

	```{r, eval=FALSE, message=FALSE, warning=FALSE, results = 'hide'}
	devtools::install_github("kbrulois/magicBatch")
	```

## Usage

```{r, eval=FALSE, message=FALSE, warning=FALSE, results = 'hide'}
python_path <- system("which python3", intern = TRUE) 
#python_path <- "path/to/venv/bin/python"
#python_path <- "/path/to/python3.x"

sce <- readRDS(url('https://stacks.stanford.edu/file/druid:cf352cg6610/PLN123_SCE.rds'))

MAGIC_w_correction <- magicBatch(data = as.matrix(t(logcounts(sce))), 
                                 mar_mat_input = reducedDim(sce, "MNN_correction"),
                                 t_param = 6, 
                                 python_command = python_path)
```

## Tutorial

Please see package website for full documentation:

https://kbrulois.github.io/magicBatch/articles/batch_correction.html

