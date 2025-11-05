

# utility function to output the string 'x (p%)' where p=100*x/total
# (for use in summary tables)
format_value_and_pct <- function(x, total) {
  paste0(as.character(x),
         ' (',
         round(100*x/total, 1) %>% format(nsmall=1) %>% as.character(),
         '%)') }

# utility function to format gene + allele columns
format_allele <- function(x) {
  y <- gsub('Homsap ', '', x)
  y <- gsub(', or', ',', y)
  y <- gsub(' or', ',', y)
  y <- gsub(' \\(see comment\\)', '', y)
  y <- gsub(' F', '', y)
  y <- gsub(' \\(F\\)', '', y)
  return(y)
}

# utility function to format gene + allele columns
format_gene <- function(x, locus='IGH') {
  z <- stringr::str_extract_all(x, paste0(locus, 'V[0-9]*-[0-9]*')) %>%
    lapply(function(l) {unique(l) %>% paste(collapse=', ')}) %>%
    unlist()
  return(z)
}

# FIXME: need to deal with case where no PBs exist
get_cell_type_table <- function(df, caption) {
  df_sum <- df %>% group_by(PTID, Visit, Cell_Type) %>%
    summarise(n_cells = n()) %>%
    ungroup() %>%
    tidyr::pivot_wider(names_from='Cell_Type', values_from='n_cells', values_fill=0)
  if ('PB' %in% names(df_sum)) {
    df_sum <- df_sum %>%
      mutate(n_cells = MBC + PB) %>%
      mutate(MBC = format_value_and_pct(MBC, n_cells),
             PB = format_value_and_pct(PB, n_cells),
             PTID = ifelse(PTID=='NotSorted', 'Not Sorted', PTID),
             Visit = ifelse(Visit==0, '', Visit)) %>%
      select(PTID, Visit, MBC, PB)
    sumrow <- data.frame(PTID = c('Total'),
                         Visit = c(''),
                         MBC = sum(df$Cell_Type == 'MBC'),
                         PB = sum(df$Cell_Type == 'PB'),
                         n_cells = nrow(df)) %>%
      mutate(MBC = format_value_and_pct(MBC, n_cells),
             PB = format_value_and_pct(PB, n_cells)) %>%
      select(PTID, Visit, MBC, PB)
    cell_type_table <- df_sum %>% rbind(sumrow) %>%
      setNames(c('PTID', 'Visit', 'No. Memory B Cells', 'No. Plasmablasts')) %>%
      kable(format='latex', booktabs=T, align='ccrr', caption=caption) %>%
      kable_styling(latex_options = "HOLD_position") %>%
      row_spec(nrow(df_sum)+1, bold=T) %>%
      collapse_rows(columns = 1:2, latex_hline = "major")
  } else {
    df_sum <- df_sum %>%
      mutate(n_cells = MBC) %>%
      mutate(MBC = format_value_and_pct(MBC, n_cells),
             # PB = format_value_and_pct(PB, n_cells),
             PTID = ifelse(PTID=='NotSorted', 'Not Sorted', PTID),
             Visit = ifelse(Visit==0, '', Visit)) %>%
      select(PTID, Visit, MBC)
    sumrow <- data.frame(PTID = c('Total'),
                         Visit = c(''),
                         MBC = sum(df$Cell_Type == 'MBC'),
                         # PB = sum(df$Cell_Type == 'PB'),
                         n_cells = nrow(df)) %>%
      mutate(MBC = format_value_and_pct(MBC, n_cells),
             # PB = format_value_and_pct(PB, n_cells)
      ) %>%
      select(PTID, Visit, MBC)
    cell_type_table <- df_sum %>% rbind(sumrow) %>%
      setNames(c('PTID', 'Visit', 'No. Memory B Cells')) %>%
      kable(format='latex', booktabs=T, align='ccrr', caption=caption) %>%
      kable_styling(latex_options = "HOLD_position") %>%
      row_spec(nrow(df_sum)+1, bold=T) %>%
      collapse_rows(columns = 1:2, latex_hline = "major")
  }
  return(cell_type_table)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

my_column_summary <- function(x) {
  y <- table(x)
  paste0(names(y), ' (', y, ')', collapse=', ')
}

my_max <- function(v) {
  if (all(is.na(v))) {
    return(NA)
  } else {
    return(max(v, na.rm=T) %>% round(1))
  }
}

my_paste_unique <- function(v) {
  if (all(is.na(v))) {
    return(NA)
  } else {
    v_filtered <- v[!is.na(v)]
    return(sort(unique(v_filtered) %>% paste(collapse=', ')))
  }
}

visc_multinom_ci <- function(n, conf.level=0.95, sides='two.sided', method='wilson') {
  DescTools::MultinomCI(n,
                        conf.level = conf.level,
                        sides = sides,
                        method = method) %>%
    as.data.frame() %>%
    mutate(var = NA, method = 'multinomial')
  select(mean = est, var, lower = lwr.ci, upper = upr.ci, method)
}

optimize_single_gene_pseudo_likelihood <- function(P, # vector of empirical binomial proportions (one per subject, all for the same gene)
                                                   N, # vector of observed counts (one per subject, same ordering as P, all for the same gene)
                                                   optim_method="L-BFGS-B") {

  # drop cases with fewer than two observations for fitting the model
  # (avoids infinite values in compute_pseudo_likelihood)
  P_obs <- P[N>0]
  N_obs <- N[N>0]

  # define the pseudo likelihood
  compute_pseudo_likelihood <- function(param) {
    mu <- param[1] # mean parameter
    sigma2 <- param[2] # variance parameter
    V <- sigma2 + mu*(1-mu)/N_obs # divided by N_obs from original
    pseudo_likelihood <- sum( (P_obs-mu)^2/V + log(V) )
    return(pseudo_likelihood)
  }

  # optimize the pseudo likelihood
  init_mu <- mean(P_obs)
  init_sigma2 <- max(var(P_obs), 0.0000001, na.rm = T) # force to be non-zero
  est <- optim(par = c(init_mu, init_sigma2), # initial estimates
               compute_pseudo_likelihood,
               method = optim_method,
               lower = c(0, 0.00000001),
               upper = c(1, 1))
  mu_hat   <- est$par[1]
  sigma2_hat <- est$par[2]
  return(c(mu_hat, sigma2_hat))
}

get_pl_est_and_bootstrap_ci <- function(P,
                                        N,
                                        n_bootstrap = 1000,
                                        conf.level = 0.95,
                                        optim_method = "L-BFGS-B") {

  # drop participants with no observations for fitting the model
  P_obs <- P[N>0]
  N_obs <- N[N>0]

  if (sum(P_obs) == 0) {

    print('No observations of this gene in this dataset')
    prop_res <- data.frame(empirical_mean = NA,
                           mean = NA,
                           lower = NA,
                           upper = NA,
                           method = NA)

  } else if (length(N_obs) < 3) {

    print('Not enough samples to use bootstrap, reporting empirical frequencies only')
    prop_res <- data.frame(empirical_mean = mean(P_obs),
                           mean = mean(P_obs),
                           lower = NA,
                           upper = NA,
                           method = 'Empirical Mean',
                           n_seq = sum(P_obs * N_obs),
                           n_samples = length(N_obs),
                           n_bootstrap = NA)

  } else {

    # get the main estimate
    est <- optimize_single_gene_pseudo_likelihood(P_obs, N_obs, optim_method=optim_method)
    est_prop <- est[1]
    est_sigma2 <- est[2]

    # get two-sided bootstrap CIs
    alpha <- 1 - conf.level
    n <- length(P_obs)
    bootstrap_result <- matrix(NA, n_bootstrap, 2)
    for (i in 1:n_bootstrap){

      # sample individuals
      sample_indices <- sample(x = 1:n, replace = TRUE, prob = NULL)
      N_bootstrap <- N_obs[sample_indices]

      # sample within individuals
      P_bootstrap <- c()
      for (j in sample_indices) {
        P_j_bootstrap <- rbinom(n = 1, size = N_obs[j], prob = P_obs[j]) / N_obs[j]
        P_bootstrap <- c(P_bootstrap, P_j_bootstrap)
      }
      P_bootstrap <- as.vector(P_bootstrap)

      bootstrap_estimate <- optimize_single_gene_pseudo_likelihood(P_bootstrap,
                                                                   N_bootstrap,
                                                                   optim_method)
      bootstrap_result[i, 1:2] <- c(bootstrap_estimate[1], bootstrap_estimate[2])

    }
    prop_lower <- quantile(bootstrap_result[,1], probs = alpha/2)
    prop_upper <- quantile(bootstrap_result[,1], probs = 1-alpha/2)
    sigma2_lower <- quantile(bootstrap_result[,2], probs = alpha/2)
    sigma2_upper <- quantile(bootstrap_result[,2], probs = 1-alpha/2)

    prop_res <- data.frame(empirical_mean = mean(P_obs),
                           mean = est_prop,
                           var = est_sigma2,
                           lower = prop_lower,
                           upper = prop_upper,
                           method = 'Fitted Mean (PL Estimate)',
                           n_seq = sum(P_obs * N_obs),
                           n_samples = length(N_obs),
                           n_bootstrap = n_bootstrap)
    row.names(prop_res) <- NULL
  }

  return(prop_res)
}

compute_gene_prevalence <- function(gene_vec,
                                    uid_vec,
                                    method = 'pl',
                                    drop_distal = TRUE) {

  # break ties in gene column apart
  df <- data.frame(gene = gene_vec, uid = uid_vec) %>%
    separate_longer_delim(gene, delim = ',') %>%
    filter(!is.na(gene))

  if (drop_distal) {
    df <- df %>% filter(!str_detect(gene, 'D'))
  }

  if (method == 'multinom') {

    df <- df %>%
      group_by(gene) %>%
      summarize(n = n(), .groups = 'drop') %>%
      distinct() %>%
      mutate(visc_multinom_ci(n)) %>%
      ungroup() %>%
      mutate(gene = factor(gene))

  } else if (method == 'pl') {

    df <- df %>%
      group_by(uid, gene) %>%
      summarize(n_gene = n(), .groups = 'drop')

    # add in "zero observed" counts for combinations of participants and genes that aren't currently in data
    df <- df %>%
      pivot_wider(names_from = gene,
                  names_prefix = 'GENE_',
                  values_from = n_gene,
                  values_fill = 0) %>%
      pivot_longer(starts_with('GENE_'),
                   names_prefix = 'GENE_',
                   names_to = 'gene',
                   values_to = 'n_gene')

    df <- df %>%
      group_by(uid) %>%
      mutate(P = n_gene / sum(n_gene),
             N = sum(n_gene)) %>%
      ungroup() %>%
      distinct() %>%
      group_by(gene) %>%
      mutate(get_pl_est_and_bootstrap_ci(P,N)) %>%
      ungroup() %>%
      mutate(gene = factor(gene))
  }

  return(df)

}
