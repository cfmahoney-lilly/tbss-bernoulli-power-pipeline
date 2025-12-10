#' Propensity Score Matching with Fixed Matching Ratio
#'
#' Performs propensity score matching using the MatchIt package to create
#' matched cohorts for treatment effect estimation. The function matches
#' treated and comparator patients based on propensity scores derived from
#' covariates, then calculates the exposure proportion in the matched sample.
#'
#' @param input_data A data frame or tibble containing patient data as returned
#'   by \code{plasmode()}. Must contain columns: \code{id}, \code{treatment_group},
#'   \code{outcome}, and covariate columns. The data should have one row per
#'   patient-outcome combination.
#' @param match_args A named list of arguments to pass to \code{MatchIt::matchit()}.
#'   Must include at minimum:
#'   \describe{
#'     \item{formula}{Formula specifying treatment ~ covariates (e.g., 
#'       \code{treatment_group ~ age + sex})}
#'     \item{method}{Matching method (e.g., "nearest", "optimal", "full")}
#'   }
#'   Additional arguments may include: \code{distance}, \code{ratio}, 
#'   \code{caliper}, \code{replace}, etc. See \code{?MatchIt::matchit} for details.
#'
#' @return A data frame containing only matched patients with the following columns:
#'   \describe{
#'     \item{id}{Patient identifiers}
#'     \item{treatment_group}{Treatment assignment (0 or 1)}
#'     \item{outcome}{Outcome codes}
#'     \item{...}{All covariate columns from input}
#'     \item{exp_p}{Numeric proportion of treated patients in matched sample}
#'   }
#'   The returned data maintains the one-row-per-patient-outcome structure.
#'
#' @details
#' The matching process:
#' \enumerate{
#'   \item Extracts unique patient records (removing duplicate outcomes per patient)
#'   \item Estimates propensity scores based on specified covariates
#'   \item Performs matching according to specified method and parameters
#'   \item Filters original data to matched patients only
#'   \item Calculates exposure proportion: treated / (treated + comparator)
#'   \item Adds exposure proportion as \code{exp_p} column to all rows
#' }
#'
#' **Requirements:**
#' \itemize{
#'   \item The \code{MatchIt} package must be installed
#'   \item \code{input_data} must contain both treatment groups (0 and 1)
#'   \item Covariate columns specified in formula must exist in data
#'   \item At least one patient from each group must remain after matching
#' }
#'
#' @export
ps_match <- function(input_data, match_args) {
  
  # Input validation
  
  # Check input_data is provided
  if (missing(input_data)) {
    stop("Argument 'input_data' is missing with no default.", call. = FALSE)
  }
  
  # Check match_args is provided
  if (missing(match_args)) {
    stop("Argument 'match_args' is missing with no default.", call. = FALSE)
  }
  
  # Check input_data is a data frame
  if (!is.data.frame(input_data)) {
    stop("'input_data' must be a data frame or tibble.", call. = FALSE)
  }
  
  # Check input_data is not empty
  if (nrow(input_data) == 0) {
    stop("'input_data' is empty (contains no rows).", call. = FALSE)
  }
  
  # Check required columns exist
  required_cols <- c("id", "treatment_group", "outcome")
  missing_cols <- setdiff(required_cols, names(input_data))
  
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "'input_data' is missing required column(s): %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  # Check match_args is a list
  if (!is.list(match_args)) {
    stop("'match_args' must be a list.", call. = FALSE)
  }
  
  # Check match_args is not empty
  if (length(match_args) == 0) {
    stop("'match_args' cannot be empty.", call. = FALSE)
  }
  
  # Check match_args has names
  if (is.null(names(match_args)) || any(names(match_args) == "")) {
    stop("All elements in 'match_args' must be named.", call. = FALSE)
  }
  
  # Check formula exists in match_args
  if (!"formula" %in% names(match_args)) {
    stop(
      "'match_args' must include 'formula' element (e.g., treatment_group ~ age + sex).",
      call. = FALSE
    )
  }
  
  # Check formula is a formula object
  if (!inherits(match_args$formula, "formula")) {
    stop("'match_args$formula' must be a formula object.", call. = FALSE)
  }
  
  # Check method exists in match_args
  if (!"method" %in% names(match_args)) {
    stop(
      "'match_args' must include 'method' element (e.g., 'nearest', 'optimal', 'full').",
      call. = FALSE
    )
  }
  
  # Check MatchIt package is available
  if (!requireNamespace("MatchIt", quietly = TRUE)) {
    stop(
      "Package 'MatchIt' is required but not installed. Install it with: install.packages('MatchIt')",
      call. = FALSE
    )
  }
  
  # Check both treatment groups exist
  unique_treatments <- unique(input_data$treatment_group)
  if (!all(c(0, 1) %in% unique_treatments)) {
    stop(
      sprintf(
        "'input_data' must contain both treatment groups (0 and 1).\nFound: %s",
        paste(unique_treatments, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  # Extract covariate names from formula
  formula_terms <- all.vars(match_args$formula)
  covariate_cols <- setdiff(formula_terms, "treatment_group")
  
  # Check that covariates exist in data
  missing_covars <- setdiff(covariate_cols, names(input_data))
  if (length(missing_covars) > 0) {
    stop(
      sprintf(
        "Covariates specified in formula not found in 'input_data': %s\nAvailable columns: %s",
        paste(missing_covars, collapse = ", "),
        paste(names(input_data), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  # Prepare data for matching
  
  # Extract unique patient records (one row per patient)
  # Keep all columns except outcome for matching
  cov <- input_data %>%
    dplyr::select(-outcome) %>%
    dplyr::distinct()
  
  # Check for duplicate patient IDs
  n_patients <- dplyr::n_distinct(cov$id)
  if (nrow(cov) != n_patients) {
    stop(
      sprintf(
        "Data contains patients with inconsistent covariate values.\nExpected %d unique patients, found %d rows after removing outcomes.",
        n_patients,
        nrow(cov)
      ),
      call. = FALSE
    )
  }
  
  # Convert to data frame (MatchIt requirement)
  cov_df <- as.data.frame(cov)
  
  # Check for missing values in covariates
  missing_in_covars <- sapply(cov_df[covariate_cols], function(x) sum(is.na(x)))
  if (any(missing_in_covars > 0)) {
    covars_with_na <- names(missing_in_covars[missing_in_covars > 0])
    warning(
      sprintf(
        "Covariates contain missing values:\n%s",
        paste(sprintf("  %s: %d missing", covars_with_na, missing_in_covars[covars_with_na]), 
              collapse = "\n")
      ),
      call. = FALSE
    )
  }
  
  # Perform propensity score matching 
  
  # Combine match_args with data
  args_list <- c(match_args, list(data = cov_df))
  
  # Call matchit with error handling
  m_out <- tryCatch(
    {
      do.call(MatchIt::matchit, args_list)
    },
    error = function(e) {
      stop(
        sprintf(
          "MatchIt::matchit() failed with error:\n%s\n\nCheck your match_args specification.",
          e$message
        ),
        call. = FALSE
      )
    }
  )
  
  # Extract matched data
  m_data <- tryCatch(
    {
      MatchIt::match.data(m_out)
    },
    error = function(e) {
      stop(
        sprintf(
          "MatchIt::match.data() failed with error:\n%s",
          e$message
        ),
        call. = FALSE
      )
    }
  )
  
  # Check if any matches were found
  if (nrow(m_data) == 0) {
    stop(
      "Matching resulted in zero matched patients. Consider:\n  - Relaxing caliper restrictions\n  - Changing matching method\n  - Checking data quality",
      call. = FALSE
    )
  }
  
  # Extract matched patient IDs
  matched_ids <- m_data$id
  
  # Filter outcomes to matched patients only
  matched_outcomes <- input_data %>%
    dplyr::filter(id %in% matched_ids)
  
  # Check both treatment groups remain after matching
  matched_groups <- unique(matched_outcomes$treatment_group)
  if (!all(c(0, 1) %in% matched_groups)) {
    warning(
      sprintf(
        "Matching eliminated one treatment group. Remaining group(s): %s",
        paste(matched_groups, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  # Calculate exposure proportion
  
  # Count unique patients by treatment group
  group_counts <- matched_outcomes %>%
    dplyr::distinct(id, treatment_group) %>%
    dplyr::count(treatment_group, name = "n_patients")
  
  # Extract counts
  matched_treat_n <- group_counts %>%
    dplyr::filter(treatment_group == 1) %>%
    dplyr::pull(n_patients)
  
  matched_comp_n <- group_counts %>%
    dplyr::filter(treatment_group == 0) %>%
    dplyr::pull(n_patients)
  
  # Handle missing groups
  if (length(matched_treat_n) == 0) matched_treat_n <- 0
  if (length(matched_comp_n) == 0) matched_comp_n <- 0
  
  # Calculate exposure proportion
  matched_p <- matched_treat_n / (matched_treat_n + matched_comp_n)
  
  # Add exposure proportion to data
  final_matched <- matched_outcomes %>%
    dplyr::mutate(exp_p = matched_p)
  
  # Summary message
  message(sprintf(
    "Matching complete: %d treated, %d comparator (%.1f%% treated)",
    matched_treat_n,
    matched_comp_n,
    matched_p * 100
  ))
  
  return(final_matched)
}

