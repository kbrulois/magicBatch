


#' Impute data using the MAGIC algorithm
#' 
#' This is an implementation of the Marcov Affinity-based Graph Imputation of Cells (MAGIC) algorithm described in 
#' Van Dijk, David et al, modified for speed and flexibility. If the mar_mat_input argument is used, the diffusion operator is computed independently of the data to be imputed. 
#' This allows any low-dimensional representation of the data, including batch-corrected data, to be directly used to calculate the powered Marcov 
#' affinity matrix. 
#' 
#' @param data An expression matrix where cells correspond to rows and genes correspond to columns
#' @param mar_mat_input A matrix where cells correspond to rows and components or features correspond to columns. If left unspecified, the Marcov matrix calculation is initialized with PCA of data.
#' @param select_features A vector of features to use for imputation.
#' @param import_mar_mat Whether to return the Markov matrix
#' @param pca An integer specifying the number of PCA components that should be used
#' @param t_param An integer or a vector of integers to be used to power the marcov affinity matrix
#' @param k The number of nearest neighbors used to construct the knn graph
#' @param ka This controls the standard deviation used in the Gaussian kernel width for a given cell, which is set to the distance to the ka-th nearest neighbor.
#' @param epsilon Epsilon parameter used in MAGIC
#' @param rescale_percent Percentile to rescale data to after imputation
#' @param rescale_method A string passed to the rescale_method argument of the rescale_data function. Two methods are available: "adaptive" or "classic" See rescale_data function for details.
#' @param n_diffusion_components Number of diffusion map components to compute. If set to 0, this diffusion map will not be computed.
#' @param python_command A character string passed to the "command" argument of the system2 function in order to invoke python. E.g. "/usr/local/bin/python3" on a Mac.
#' @return A \code{list} that includes the following elements:
#' \describe{
#'   \item{imputed_data}{A cell by gene \code{matrix} of the imputed gene expression values.}
#'   \item{diffusion_map}{A cell by diffusion map component \code{matrix}.}
#'   \item{marcov_matrix}{A cell by cell \code{matrix} of the unpowered markov affinity matrix.}
#' }
#'
#' @author Kevin Brulois
#' @export

magicBatch <- function(data,
                       select_features = NULL,
                       mar_mat_input = NULL,
                       import_mar_mat = FALSE,
                       use_numpy = TRUE,
                       pca = 20,
                       t_param = c(2,4,6),
                       n_diffusion_components = 0,
                       k = 9, 
                       ka = 3, 
                       epsilon = 1, 
                       rescale_percent = 90,
                       rescale_method = "adaptive",
                       python_command = system("which python3", intern = TRUE)) {
  
  assertthat::assert_that(!(!is.null(select_features) & !is.null(mar_mat_input)))
  
  on.exit({print(paste("removing temporary files"))
    try({file.remove(to.remove)})
    print(paste("done") )})
  
  data <- as.matrix(data)
  
  num_rows <- nrow(data)
  name_rows <- rownames(data)
  
  if(!is.null(select_features)) {
    assertthat::assert_that(class(select_features) %in% c("character", "integer", "numeric", "logical"))
    if(is.character(select_features)) {
      select_features_sub <- select_features[select_features %in% colnames(data)]
      features_removed <- setdiff(select_features, select_features_sub)
      if(length(features_removed > 0)) {
        message("Warning the following select_features do not match any colnames of data: \n",
                paste(features_removed, collapse = " "))
      }
      inds <- colnames(data) %in% select_features_sub
    } else {
      assertthat::are_equal(length(select_features), ncol(data))
      if(is.logical(select_features)) {
        inds <- select_features
      }
      else {
        inds <- 1:ncol(data) %in% select_features
      }
    }
  } 
  
  message("exporting data to python")
  
  path2MagScript <- system.file("exec", "MAGIC.py", package="magicBatch")
  
  path2Mar_matOut <- tempfile("mar_mat", fileext = ".csv")
  to.remove <-path2Mar_matOut
  
  path2Dif_mapOut <- tempfile("dif_map", fileext = ".csv")
  if(n_diffusion_components != 0) {
    to.remove <- c(to.remove, path2Dif_mapOut)
  }
  
  path2Data_Out <- tempfile("dat_mat", fileext = ".csv")
  if(use_numpy) {
    to.remove <- c(to.remove, path2Data_Out)
  }
  
  magParams <- paste("-o1", path2Mar_matOut,
                     "-o2", path2Dif_mapOut,
                     "-o3", path2Data_Out,
                     "-c", n_diffusion_components,
                     "-p", pca,
                     "-k", k,
                     "-ka", ka,
                     "-e", epsilon,
                     "-r", rescale_percent)
  
  if(!is.null(mar_mat_input)) {
    assertthat::are_equal(nrow(data), nrow(mar_mat_input))
    path2Aff_matIn <- tempfile("mar_mat_in", fileext = ".csv")
    data.table::fwrite(as.data.frame(mar_mat_input), file = path2Aff_matIn, row.names = FALSE)
    magParams <- paste(magParams, "-m", path2Aff_matIn)
    to.remove <- c(to.remove, path2Aff_matIn)
  }
  if((is.null(mar_mat_input) | use_numpy) & is.null(select_features)) {
    path2DataIn <- tempfile("dat_mat_in", fileext = ".csv")
    data.table::fwrite(as.data.frame(data), file = path2DataIn, row.names = FALSE)
    magParams <- paste(magParams, "-d", path2DataIn)
    to.remove <- c(to.remove, path2DataIn)
  }
  if(!is.null(select_features)) {
    path2DataIn <- tempfile("dat_mat_in", fileext = ".csv")
    data.table::fwrite(as.data.frame(rbind(fzJKy2SHmdwIPQ = inds, data)), file = path2DataIn, row.names = TRUE)
    magParams <- paste(magParams, "-d", path2DataIn)
    to.remove <- c(to.remove, path2DataIn)
  }
  if(use_numpy) {
    magParams <- paste(magParams, "-t", paste(t_param, collapse = "_"))
  }
  
  message("computing Marcov matrix")
  output <- system2(python_command, args= paste(path2MagScript, magParams), stdout=TRUE)
  print(paste(output))
  message("importing data from python to R")
  if(import_mar_mat) {
    L <- data.table::fread(path2Mar_matOut)
    L <- Matrix::sparseMatrix(i = L$row + 1, j = L$col + 1, x = L$data, dims = c(num_rows, num_rows), dimnames = list(name_rows, name_rows))
  } else {
    L <- NULL
  }
  if(n_diffusion_components != 0) {
    diffusion_map <- as.matrix(data.table::fread(path2Dif_mapOut, header = TRUE)[,-1])
    rownames(diffusion_map) <- rownames(data)
  } else {
    diffusion_map <- NULL
  }
  
  if(use_numpy) {
    message("importing imputed data from python to R")
    imputed_data <- as.matrix(data.table::fread(path2Data_Out, header = TRUE))
    cols <- ncol(imputed_data) / length(t_param)
    imputed_data <- Map(function(x) {
      temp <- imputed_data[,1:cols + cols * (x - 1)]
      dimnames(temp) <- dimnames(data)
      temp
    }, seq_along(t_param))
    names(imputed_data) <- paste0("t", t_param)
  } else {
    message("imputing data")
    data <- as(as.matrix(data), "dgCMatrix")
    imputed_data <- multi_t_fast_impute(data, L, t_param)
  }
  if(rescale_percent != 0 & rescale_method %in% c("classic", "adaptive")) {
    imputed_data <- lapply(imputed_data, function(x) rescale_data(data = data, 
                                                                  imputed_data = x, 
                                                                  rescale_percent = rescale_percent, 
                                                                  rescale_method = rescale_method))
  }
  return(list(imputed_data = imputed_data,
              marcov_mat = L,
              diffusion_map = diffusion_map))
}

#' Impute data using the MAGIC algorithm
#' 
#' This is a helper function that imputes data using a pre-computed Marcov matrix. It minimizing the number of matrix multiplications needed to 
#' compute a given set of matrix powering operations. The exponents (t_param) must be integers from 1 to 36.
#' 
#' @param data An expression matrix where cells correspond to rows and genes correspond to columns
#' @param marcov_mat A Markov matrix in dgCMatrix format where cells correspond to rows and to columns
#' @param t_param A vector of integers ranging from 1 to 32
#' @return A \code{list} that includes the following elements:
#' \describe{
#'   \item{imputed_data}{A cell by gene \code{matrix} of the imputed gene expression values.}
#'   \item{diffusion_map}{A cell by diffusion map component \code{matrix}.}
#'   \item{affinity_matrix}{A cell by cell \code{matrix} of the unpowered markov affinity matrix.}
#' }
#'
#' @author Kevin Brulois
#' @export
#' 
multi_t_fast_impute <- function(data, L, t_param) {
  
  t_max <- max(t_param)
  
  if(t_max >= 2) L_2 <- L %*% L
  if(t_max >= 4) L_4 <- L_2 %*% L_2
  if(t_max >= 8) L_8 <- L_4 %*% L_4
  if(t_max >= 16) L_16 <- L_8 %*% L_8
  if(t_max >= 32) L_32 <- L_16 %*% L_16
  if(sum(t_param %in% c(3,7,11,15,19,23,27,31,35)) > 0) L_3 <- L_2 %*% L
  
  dif_operator <- expression(L,
                             L_2,
                             L_2 %*% L,
                             L_4,
                             L_4 %*% L,
                             L_4 %*% L_2,
                             L_4 %*% L_3,
                             L_8,
                             L_8 %*% L,
                             L_8 %*% L_2,
                             L_8 %*% L_3,
                             L_8 %*% L_4,
                             L_8 %*% L_4 %*% L,
                             L_8 %*% L_4 %*% L_2,
                             L_8 %*% L_4 %*% L_3,
                             L_16,
                             L_16 %*% L,
                             L_16 %*% L_2,
                             L_16 %*% L_3,
                             L_16 %*% L_4,
                             L_16 %*% L_4 %*% L,
                             L_16 %*% L_4 %*% L_2,
                             L_16 %*% L_4 %*% L_3,
                             L_16 %*% L_8,
                             L_16 %*% L_8 %*% L,
                             L_16 %*% L_8 %*% L_2,
                             L_16 %*% L_8 %*% L_3,
                             L_16 %*% L_8 %*% L_4,
                             L_16 %*% L_8 %*% L_4 %*% L,
                             L_16 %*% L_8 %*% L_4 %*% L_2,
                             L_16 %*% L_8 %*% L_4 %*% L_3,
                             L_32,
                             L_32 %*% L,
                             L_32 %*% L_2,
                             L_32 %*% L_3,
                             L_32 %*% L_4
  )
  
  imputed_data <- lapply(t_param, function(x) as.matrix(eval(dif_operator[[x]]) %*% data))
  names(imputed_data) <- paste0("t", t_param)
  return(imputed_data)
  
}


#' Rescale Data
#' 
#' This function rescales imputed data using one of two methods. The "classic" method rescales the imputed data as in the original implementation.
#' The "adaptive" method adjusts rescaling factors for each gene using only cells with a non-zero expression value. 
#' 
#' @param data The original unimputed expression matrix where cells correspond to rows and genes correspond to columns
#' @param imputed_data Unscaled imputed data where cells corrsponnd to rows and genes correspond to columns
#' @param rescale_percent Percentile used to rescale imputed data.
#' @param rescale_method A string specifying one of three rescale methods: "adaptive" or "classic" The "classic" method performs rescaling as in the original version of MAGIC: rescaling factors are computed using the given percentile of the original and imputed data or the max if the data is zero at the given percentile. The "adaptive" method adjusts the rescale percentile for each gene between recale_percent and 100 based on the number of non-zero data points.
#' @return A cell by gene \code{matrix} with the rescaled imputed gene expression values.
#' 
#' @author Kevin Brulois
#' @export
#' 
rescale_data <- function(data, imputed_data, rescale_percent, rescale_method) {
  
  rescale_percent <- rescale_percent / 100
  num_col <- ncol(data)
  num_row <- nrow(data)
  
  if(rescale_method == "adaptive") {
    z <- apply(data, 2, function(x) sum(x == 0)/num_row)
    nz <- 1 - z
    rescale_percent_adj <- rescale_percent * nz + z
    M99 <- list()
    for(i in 1:num_col) M99[[i]] <- quantile(data[,i], rescale_percent_adj[i])
    M99 <- do.call(c, M99)
    maxes <- apply(data, 2, max)
    inds <- M99 == 0
    M99[inds] <- maxes[inds]
    M99_new <- list()
    for(i in 1:num_col) M99_new[[i]] <- quantile(imputed_data[,i], rescale_percent_adj[i])
    M99_new <- do.call(c, M99_new)
    maxes <- apply(imputed_data, 2, max)
    inds <- M99_new == 0
    M99_new[inds] <- maxes[inds]
  } 
  
  if(rescale_method == "classic") {
    M99 <- apply(data, 2, function(x) quantile(x, rescale_percent))
    maxes <- apply(data, 2, max)  
    inds <- M99 == 0
    M99[inds] <- maxes[inds]
    M99_new <- apply(imputed_data, 2, function(x) quantile(x, rescale_percent))
    maxes <- apply(data, 2, max)  
    inds <- M99_new == 0
    M99_new[inds] <- maxes[inds]
  }
  ratio = M99 / M99_new
  t(t(imputed_data) * ratio)
}


