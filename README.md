Batch-corrected Imputation================

## Overview

This is a modified version of the MAGIC algorithm that allows for batch-corrected imputation of multi-batch datasets.

## Installation

First install the magicBatch python package from within R by specifying the python_path using one of the following:

```{r, eval=FALSE, message=FALSE, warning=FALSE, results = 'hide'}

python_path <- system("which python3", intern = TRUE) #or "path/to/venv/bin/python" or "/path/to/python3.x"
system(paste(python_path, "-m pip install magicBatch"))
```

Next install the magicBatch R package:

```{r, eval=FALSE, message=FALSE, warning=FALSE, results = 'hide'}
devtools::install_github("kbrulois/magicBatch")

```

## Tutorial

Please see package website for full documentation:

https://kbrulois.github.io/magicBatch/articles/batch_correction.html

