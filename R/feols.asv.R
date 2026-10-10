#' Average Shared Variance (ASV) Decomposition for Two-Way Fixed-Effects Models
#' Estimated with `fixest::feols`
#'
#' @param mod A fitted `fixest` model estimated using `fixest::feols`. The model
#'   must include at least one fixed effect and should contain at least two
#'   explanatory variables unless predictor groups are supplied through `iv`.
#' @param iv Optional named list defining groups of explanatory variables. Each
#'   element should contain the names of one or more predictors in the model.
#'   ASV contributions are calculated for the specified predictor groups. If
#'   `NULL` (default), predictors are evaluated individually.
#' @param commonality Logical. If `TRUE`, return the commonality decomposition
#'   for all non-empty combinations of predictors or predictor groups. If
#'   `FALSE` (default), return the ASV decomposition.
#'
#' @details This function implements the average shared variance (ASV) approach
#'   for linear models fitted with `fixest::feols`, including models with
#'   two-way fixed effects. It fits all non-empty combinations of explanatory
#'   variables (or groups) while retaining the fixed effects specified in the
#'   original model, and uses each model's within \(R^2\) to quantify
#'   explained variation.
#'
#'   ASV allocates explained variation among predictors by assigning shared
#'   fractions equally among the predictors involved in each fraction. The
#'   resulting individual contributions are additive and order-independent;
#'   their sum equals the full model's within \(R^2\), subject to numerical
#'   precision and model-estimation constraints. This provides an interpretable
#'   allocation of explained variation when predictors share explanatory
#'   information.
#'
#'   The optional `commonality = TRUE` output provides the underlying
#'   commonality decomposition, partitioning explained variation into unique
#'   and shared fractions. Negative commonality fractions can occur.
#'
#'   Interaction terms are not currently supported. The model formula must
#'   include at least one fixed effect, and the data used to fit the model must
#'   be available when the function is called.
#'
#' @return A list of class `feolsasv` containing:
#'   \item{Within.R2}{Within \(R^2\) of the full model.}
#'   \item{hierarchical.partitioning}{ASV contributions for each predictor or
#'   predictor group, including unique contribution, allocated shared
#'   contribution, total individual contribution, and percentage of total
#'   explained variation. This legacy component name is retained for
#'   compatibility.}
#'   \item{commonality.analysis}{When `commonality = TRUE`, a matrix reporting
#'   the commonality fraction and its percentage of the total for each
#'   non-empty predictor combination, plus a total row.}
#'   \item{variables}{Names of the predictors or predictor groups analysed.}
#'   \item{type}{Indicates whether ASV or commonality analysis was returned.}
#'   \item{decomposition.error}{Difference between the sum of commonality
#'   fractions and the full model's within \(R^2\).}
#'
#' @author Jiangshan Lai \email{lai@njfu.edu.cn}
#'
#' @references
#' \itemize{
#' \item Lai, J., Zou, Y., Zhang, J. & Peres-Neto, P. (2022). Generalizing
#'   hierarchical and variation partitioning in multiple regression and
#'   canonical analyses using the rdacca.hp R package. *Methods in Ecology
#'   and Evolution*, 13(4), 782--788. doi:10.1111/2041-210X.13800
#' \item Lai, J., Zou, Y., Zhang, S., Zhang, X. & Mao, L. (2022). glmm.hp: an R
#'   package for computing individual effect of predictors in generalized
#'   linear mixed models. *Journal of Plant Ecology*, 15(6), 1302--1307.
#'   doi:10.1093/jpe/rtac096
#' \item Lai, J., Zhu, W., Cui, D. & Mao, L. (2023). Extension of the glmm.hp
#'   package to zero-inflated generalized linear mixed models and multiple
#'   regression. *Journal of Plant Ecology*, 16(6), rtad038.
#'   doi:10.1093/jpe/rtad038
#' \item Chevan, A. & Sutherland, M. (1991). Hierarchical partitioning.
#'   *The American Statistician*, 45, 90--96. doi:10.1080/00031305.1991.10475776
#' \item Nimon, K., Oswald, F. L. & Roberts, J. K. (2013). Yhat: Interpreting
#'   regression effects. R package version 2.0.0.
#' \item Nimon, K., Ho, L. S. T. & Ane, C. (2014). A linear-time algorithm for
#'   Gaussian and non-Gaussian trait evolution models. *Systematic Biology*,
#'   63(3), 397--408.
#' }
#'
#' @export
#' @examples
#' library(fixest)
#'
#' # Use base_did, a balanced panel dataset supplied with fixest.
#' # It contains 104 individuals observed over 10 periods.
#' data("base_did", package = "fixest")
#'
#' # Create the treatment-by-post interaction used in a difference-in-differences
#' # model. The main effects of treat and post are absorbed by id and period
#' # fixed effects, respectively.
#' base_did$did <- base_did$treat * base_did$post
#'
#' # Estimate a two-way fixed-effects model.
#' mod <- feols(y ~ x1 + did | id + period,
#'              data = base_did, cluster = ~id)
#'
#' # Average shared variance (ASV) decomposition of within R-squared.
#' asv <- feols.asv(mod)
#' asv
#'
#' # Optional commonality decomposition.
#' ca <- feols.asv(mod, commonality = TRUE)
#' ca
#'
#' # Evaluate grouped predictor contributions.
#' iv <- list(Control = "x1", Treatment = "did")
#' feols.asv(mod, iv = iv)

feols.asv <- function(mod, iv = NULL, commonality = FALSE) {
  
  # =========================================================
  # 0. Model checking
  # =========================================================
  
  if (!inherits(mod, "fixest")) {
    stop("feols.asv only supports fixest::feols models.")
  }
  
  
  # =========================================================
  # 1. Extract data and formula
  # =========================================================
  
  dat <- eval(
    mod$call$data,
    envir = parent.frame()
  )
  
  fml <- mod$call$fml
  
  fml_chr <- paste(
    deparse(fml),
    collapse = ""
  )
  
  parts <- strsplit(
    fml_chr,
    "\\|"
  )[[1]]
  
  if (length(parts) != 2) {
    stop(
      "feols.asv currently requires at least one fixed effect."
    )
  }
  
  
  # Linear part of the formula
  linear_fml <- as.formula(
    trimws(parts[1])
  )
  
  
  # ---------------------------------------------------------
  # IMPORTANT:
  # Keep the complete response expression.
  # For example, log(Euros) must remain log(Euros),
  # rather than being reduced to Euros.
  # ---------------------------------------------------------
  
  response <- paste(
    deparse(linear_fml[[2]]),
    collapse = ""
  )
  
  
  # Extract fixed effects
  fixef_vars <- all.vars(
    as.formula(
      paste(
        "~",
        trimws(parts[2])
      )
    )
  )
  
  
  # Extract predictor terms
  ivname <- attr(
    terms(linear_fml),
    "term.labels"
  )
  
  
  if (length(ivname) < 2 && is.null(iv)) {
    stop(
      "Analysis not conducted. ",
      "At least two predictors are required."
    )
  }
  
  
  # Interaction terms are not supported at present
  if (any(grepl(":", ivname))) {
    stop(
      "Interaction terms are not currently supported. ",
      "Please create interaction variables explicitly."
    )
  }
  
  
  # =========================================================
  # 2. Predictor grouping
  # =========================================================
  
  if (!is.null(iv)) {
    
    nvar <- length(iv)
    
    if (nvar < 2) {
      stop(
        "At least two predictor groups are required."
      )
    }
    
    ilist <- names(iv)
    
    if (is.null(ilist)) {
      
      names(iv) <- paste0(
        "X",
        seq_len(nvar)
      )
      
    } else {
      
      whichnoname <- which(
        is.na(ilist) |
          ilist == ""
      )
      
      if (length(whichnoname) > 0) {
        
        names(iv)[whichnoname] <-
          paste0(
            "X",
            whichnoname
          )
      }
    }
    
    iv.name <- names(iv)
    
  } else {
    
    nvar <- length(ivname)
    
    iv.name <- ivname
  }
  
  
  # =========================================================
  # 3. Generate all predictor combinations
  # =========================================================
  
  totalN <- 2^nvar - 1
  
  binarymx <- matrix(
    0,
    nrow = nvar,
    ncol = totalN
  )
  
  for (i in seq_len(totalN)) {
    
    binarymx[, i] <- as.integer(
      intToBits(i)[seq_len(nvar)]
    )
  }
  
  
  # =========================================================
  # 4. Full-model Within R2
  # =========================================================
  
  full_within_R2 <- as.numeric(
    fixest::fitstat(
      mod,
      "wr2"
    )
  )
  
  
  # =========================================================
  # 5. Fit all predictor combinations
  # =========================================================
  
  commonM <- matrix(
    NA_real_,
    nrow = totalN,
    ncol = 3
  )
  
  colnames(commonM) <- c(
    "Model",
    "R2",
    "Commonality"
  )
  
  
  for (i in seq_len(totalN)) {
    
    selected <- which(
      binarymx[, i] == 1
    )
    
    
    # -------------------------------------------------------
    # Ungrouped predictors
    # -------------------------------------------------------
    
    if (is.null(iv)) {
      
      tmp.name <- iv.name[selected]
      
    }
    
    
    # -------------------------------------------------------
    # Grouped predictors
    # -------------------------------------------------------
    
    else {
      
      selected_groups <- iv.name[selected]
      
      tmp.name <- unique(
        unlist(
          iv[selected_groups]
        )
      )
    }
    
    
    rhs <- paste(
      tmp.name,
      collapse = " + "
    )
    
    
    # -------------------------------------------------------
    # Construct reduced model
    # -------------------------------------------------------
    
    fml_new <- as.formula(
      paste(
        response,
        "~",
        rhs,
        "|",
        paste(
          fixef_vars,
          collapse = " + "
        )
      )
    )
    
    
    modnew <- fixest::feols(
      fml_new,
      data = dat
    )
    
    
    commonM[i, "Model"] <- i
    
    commonM[i, "R2"] <- as.numeric(
      fixest::fitstat(
        modnew,
        "wr2"
      )
    )
  }
  
  
  # =========================================================
  # 6. Calculate commonality coefficients
  #
  # C(S) =
  # sum[(-1)^(|S|-|T|) * R2(T)]
  #
  # =========================================================
  
  for (i in seq_len(totalN)) {
    
    selected <- which(
      binarymx[, i] == 1
    )
    
    n_selected <- length(selected)
    
    subsetN <- 2^n_selected - 1
    
    ccsum <- 0
    
    
    for (j in seq_len(subsetN)) {
      
      local_bits <- as.integer(
        intToBits(j)[seq_len(n_selected)]
      )
      
      selected_local <- which(
        local_bits == 1
      )
      
      global_bits <- rep(
        0,
        nvar
      )
      
      global_bits[
        selected[
          selected_local
        ]
      ] <- 1
      
      
      # Find the corresponding model
      model_index <- which(
        apply(
          binarymx,
          2,
          function(x)
            all(x == global_bits)
        )
      )
      
      
      if (length(model_index) != 1) {
        stop(
          "Internal error while calculating commonality."
        )
      }
      
      
      sign <- (-1)^(
        n_selected -
          sum(local_bits)
      )
      
      
      ccsum <- ccsum +
        sign *
        commonM[
          model_index,
          "R2"
        ]
    }
    
    
    commonM[
      i,
      "Commonality"
    ] <- ccsum
  }
  
  
  # =========================================================
  # 7. Order models by number of predictors
  # =========================================================
  
  orderList <- order(
    colSums(binarymx),
    seq_len(totalN)
  )
  
  
  commonM <- commonM[
    orderList,
    ,
    drop = FALSE
  ]
  
  
  binarymx <- binarymx[
    ,
    orderList,
    drop = FALSE
  ]
  
  
  # =========================================================
  # 8. Generate glmm.asv-style commonality names
  # =========================================================
  
  rowNames <- character(
    totalN + 1
  )
  
  
  for (i in seq_len(totalN)) {
    
    nbits <- sum(
      binarymx[, i]
    )
    
    selected <- which(
      binarymx[, i] == 1
    )
    
    
    if (nbits == 1) {
      
      rowNames[i] <- paste(
        "Unique to",
        iv.name[selected]
      )
      
    } else {
      
      vars <- iv.name[selected]
      
      if (nbits == 2) {
        
        rowNames[i] <- paste(
          "Common to",
          vars[1],
          "and",
          vars[2]
        )
        
      } else {
        
        rowNames[i] <- paste(
          "Common to",
          paste(
            vars[-length(vars)],
            collapse = ", "
          ),
          "and",
          vars[length(vars)]
        )
      }
    }
  }
  
  
  rowNames[
    totalN + 1
  ] <- "Total"
  
  
  # =========================================================
  # 9. Commonality analysis output
  # =========================================================
  
  totalRSquare <- sum(
    commonM[
      ,
      "Commonality"
    ]
  )
  
  
  outputcommonM <- matrix(
    NA_real_,
    nrow = totalN + 1,
    ncol = 2
  )
  
  
  colnames(outputcommonM) <- c(
    "Fractions",
    "% Total"
  )
  
  
  rownames(outputcommonM) <- rowNames
  
  
  outputcommonM[
    seq_len(totalN),
    "Fractions"
  ] <- round(
    commonM[
      ,
      "Commonality"
    ],
    4
  )
  
  
  outputcommonM[
    seq_len(totalN),
    "% Total"
  ] <- round(
    100 *
      commonM[
        ,
        "Commonality"
      ] /
      totalRSquare,
    2
  )
  
  
  outputcommonM[
    totalN + 1,
    "Fractions"
  ] <- round(
    totalRSquare,
    4
  )
  
  
  outputcommonM[
    totalN + 1,
    "% Total"
  ] <- 100
  
  
  # =========================================================
  # 10. Average shared variance (ASV) allocation
  # =========================================================
  
  VariableImportance <- matrix(
    NA_real_,
    nrow = nvar,
    ncol = 4
  )
  
  
  for (i in seq_len(nvar)) {
    
    VariableImportance[
      i,
      3
    ] <-
      sum(
        binarymx[i, ] *
          (
            commonM[
              ,
              "Commonality"
            ] /
              colSums(binarymx)
          )
      )
  }
  
  
  # Unique contribution under ASV
  VariableImportance[
    ,
    1
  ] <-
    commonM[
      seq_len(nvar),
      "Commonality"
    ]
  
  
  # Shared contribution allocated equally under ASV
  VariableImportance[
    ,
    2
  ] <-
    VariableImportance[
      ,
      3
    ] -
    VariableImportance[
      ,
      1
    ]
  
  
  # Percentage of total ASV contribution
  total <- sum(
    VariableImportance[
      ,
      3
    ]
  )
  
  
  VariableImportance[
    ,
    4
  ] <-
    100 *
    VariableImportance[
      ,
      3
    ] /
    total
  
  
  dimnames(
    VariableImportance
  ) <- list(
    iv.name,
    c(
      "Unique",
      "Shared",
      "Individual",
      "I.perc(%)"
    )
  )
  
  
  # =========================================================
  # 11. Decomposition check
  # =========================================================
  
  decomposition_error <-
    totalRSquare -
    full_within_R2
  
  
  # =========================================================
  # 12. Return results
  # =========================================================
  
  outputList <- list()
  
  
  outputList$Within.R2 <-
    full_within_R2
  
  
  if (commonality) {
    
    outputList$commonality.analysis <-
      outputcommonM
    
  } else {
    
    outputList$hierarchical.partitioning <-
      VariableImportance
  }
  
  
  outputList$variables <-
    iv.name
  
  
  outputList$type <-
    ifelse(
      commonality,
      "commonality.analysis",
      "hierarchical.partitioning"
    )
  
  
  outputList$decomposition.error <-
    decomposition_error
  
  
  class(outputList) <-
    "feolsasv"
  
  
  return(outputList)
}

