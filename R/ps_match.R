#' @title Propensity score matching with fixed matching ratio.
#' @description Propensity score matching with fixed matching ratio for use in 
#'   TBSS with Bernoulli probability model. Arguments are specified in the file
#'   match_params.R
#' @return Matched data in identical format to the initial data with id, 
#'   covariate, treatment, and outcome columns.
#' @param input_data Data with patient id, treatment status, covariates, and 
#'   outcomes, as returned by plasmode().
#' @param match_args List of arguments for matching algorithm as specified in 
#'   match_params.R.
ps_match <- function(input_data, match_args) {
  cov <- input_data %>%
    dplyr::select(-outcome) %>%
    distinct()

  cov_df <- as.data.frame(cov)

  args_df <- c(match_args, list(data = cov_df))

  m_out <- do.call(matchit, args_df)

  m_data <- match_data(m_out)
  matched_ids <- m_data$id

  # separate into treatment and comparator incident outcomes
  matched_outcomes <- input_data %>%
    filter(id %in% (matched_ids))

  matched_treat_n <- matched_outcomes %>%
    filter(treatment_group == 1) %>%
    summarise(unique_ids = n_distinct(id)) %>%
    pull(unique_ids)

  matched_comp_n <- matched_outcomes %>%
    filter(treatment_group == 0) %>%
    summarise(unique_ids = n_distinct(id)) %>%
    pull(unique_ids)

  matched_p <- matched_treat_n / (matched_treat_n + matched_comp_n)

  final_matched <- matched_outcomes %>%
    mutate(exp_p = matched_p)

  return(final_matched)
}
