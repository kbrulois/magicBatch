Markov Affinity-based Graph Imputation of Cells (MAGIC)
-------------------------------------------------------
MAGIC is an imputation algorthm developed by Van Dijk, David, et al. "Recovering Gene Interactions from Single-Cell Data Using Data Diffusion." Cell (2018).

https://www.ncbi.nlm.nih.gov/pubmed/29961576

This implementation includes an R wrapper and a flexible approach for imputation across multiple batches. 



##### Installation and dependencies for the R version
The R version can be installed using:
                        
        install.packages("devtools")
        devtools::install_github("kbrulois/magic", ref = "magicBatch")

