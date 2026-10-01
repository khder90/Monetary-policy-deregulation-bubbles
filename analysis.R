# Monetary announcement shocks and bubble confidence indicators
# Run from the directory containing data.xlsx with: Rscript analysis.R

options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20260921)

packages <- c("readxl", "openxlsx", "dplyr", "tidyr", "purrr", "ggplot2",
              "moments", "tseries", "sandwich", "lmtest", "car", "mgcv",
              "strucchange")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install these CRAN packages first: ",
                                   paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages(invisible(lapply(packages, library, character.only = TRUE)))

input <- "data.xlsx"
if (!file.exists(input)) stop("Place data.xlsx beside analysis.R and run the script from that directory.")
dir.create("results", showWarnings = FALSE)
dir.create(file.path("results", "figures"), showWarnings = FALSE)

outcomes <- c("SPI", "MPI", "LPI", "SNI", "MNI", "LNI")
positive <- c("SPI", "MPI", "LPI")
shocks <- c("MP_MED_Z", "CBI_MED_Z")
controls <- c("CFNAI_L1_Z", "INF_L1_Z", "DREG_L1_Z")
state <- "DREG_L1_Z"
H <- 12L
bootstrap_reps <- 499L
block_length <- 12L
outputs <- list()
save_result <- function(name, x) { outputs[[name]] <<- as.data.frame(x); invisible(x) }
standardize <- function(x) {
  x <- as.numeric(x); s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
to_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))
  as.Date(x)
}
lag_by <- function(x, k, group = NULL, ahead = FALSE) {
  fun <- if (ahead) function(z) dplyr::lead(z, k) else function(z) dplyr::lag(z, k)
  if (is.null(group)) return(fun(x))
  ave(x, group, FUN = fun)
}
prepare_us <- function(d, components = c("PC1_HF", "MP_MED", "CBI_MED", "MP_PM", "CBI_PM")) {
  d <- d[order(to_date(d$DATE)), , drop = FALSE]
  d$DATE <- to_date(d$DATE)
  for (pair in list(c("CFNAI", "CFNAI_L1_Z"), c("INF_YOY", "INF_L1_Z"),
                    c("RISK", "RISK_L1_Z")))
    d[[pair[2]]] <- standardize(dplyr::lag(d[[pair[1]]]))
  # Preserve the observed first lag rather than recomputing it after filtering.
  d$DREG_L1_Z <- standardize(d$DREG_L1)
  for (v in c("CFNAI", "INF_YOY", "RISK", "DREG", components))
    if (v %in% names(d)) d[[paste0(v, "_Z")]] <- standardize(d[[v]])
  d
}
us <- prepare_us(as.data.frame(readxl::read_excel(input, sheet = "US_MP_CBI")))
us_u <- prepare_us(as.data.frame(readxl::read_excel(input, sheet = "US_U1U4")),
                   components = paste0("U", 1:4))
eu <- as.data.frame(readxl::read_excel(input, sheet = "EU_Panel"))
eu$DATE <- to_date(eu$DATE)
eu <- eu[order(eu$ISO, eu$DATE), , drop = FALSE]
eu <- eu |>
  dplyr::group_by(ISO) |>
  dplyr::mutate(dplyr::across(dplyr::all_of(c("DREG", "DREG_L1", "PC1_HF",
      "MP_MED", "CBI_MED", "MP_PM", "CBI_PM")),
      standardize, .names = "{.col}_Z")) |>
  dplyr::ungroup()
save_result("QC_US", data.frame(variable = names(us), N = vapply(us, function(z) sum(!is.na(z)), 0L)))

dynamic <- function(d, y, s, p, group = NULL) {
  g <- if (is.null(group)) NULL else d[[group]]
  terms <- character()
  for (j in seq_len(p)) for (v in c(y, s)) {
    nm <- paste0("L", j, "_", v)
    d[[nm]] <- lag_by(d[[v]], j, g)
    terms <- c(terms, nm)
  }
  list(data = d, terms = terms, group = g)
}
fit_lp <- function(d, y, s = shocks, x = controls, p = 1L, h = 0L,
                   group = NULL, cumulative = FALSE, transition = NULL,
                   interaction = FALSE) {
  z <- dynamic(d, y, s, p, group)
  dd <- z$data
  if (cumulative) {
    future <- sapply(0:h, function(k) lag_by(dd[[y]], k, z$group, TRUE))
    dd$Y <- ifelse(rowSums(is.na(future)) == 0, rowSums(future), NA_real_)
  } else dd$Y <- lag_by(dd[[y]], h, z$group, TRUE)
  main <- s
  if (interaction) for (v in s) {
    nm <- paste0("INT_", v)
    dd[[nm]] <- dd[[v]] * dd[[state]]
    main <- c(main, nm)
  }
  if (!is.null(transition)) {
    dd$FSTATE <- plogis(transition$gamma * (dd[[state]] - transition$c))
    main <- "FSTATE"
    for (v in s) {
      lo <- paste0(v, "_LOW"); hi <- paste0(v, "_HIGH")
      dd[[lo]] <- dd[[v]] * (1 - dd$FSTATE)
      dd[[hi]] <- dd[[v]] * dd$FSTATE
      main <- c(main, lo, hi)
    }
  }
  rhs <- unique(c(main, x, z$terms))
  keep <- complete.cases(dd[, c("Y", rhs), drop = FALSE])
  dd <- dd[keep, , drop = FALSE]
  if (nrow(dd) < length(rhs) + 10L) stop("Insufficient observations: ", y, " h=", h)
  fit <- stats::lm(stats::reformulate(rhs, "Y"), data = dd)
  V <- sandwich::NeweyWest(fit, lag = max(h + 1L, p), prewhite = FALSE, adjust = TRUE)
  list(fit = fit, V = V, data = dd, rhs = rhs)
}
coefs <- function(obj, terms) {
  b <- coef(obj$fit)[terms]; se <- sqrt(diag(obj$V)[terms])
  data.frame(term = terms, estimate = as.numeric(b), se = as.numeric(se),
    t = as.numeric(b / se), p = as.numeric(2 * pt(abs(b / se), obj$fit$df.residual,
                                                  lower.tail = FALSE)),
    lo90 = as.numeric(b - 1.645 * se), hi90 = as.numeric(b + 1.645 * se),
    lo95 = as.numeric(b - 1.96 * se), hi95 = as.numeric(b + 1.96 * se))
}
linear_path <- function(y, d, s = shocks, x = controls, p = 1L, max_h = H,
                        group = NULL) purrr::map_dfr(0:max_h, function(h) {
  obj <- fit_lp(d, y, s, x, p, h, group)
  cbind(coefs(obj, s), outcome = y, h = h, N = nobs(obj$fit))
})
common_lag <- function(d, ys = outcomes, s = shocks, x = controls, pmax = 12L) {
  purrr::map_dfr(seq_len(pmax), function(p) {
    b <- sum(vapply(ys, function(y) BIC(fit_lp(d, y, s, x, p, 1L)$fit), 0.0))
    data.frame(p = p, Total_BIC = b, Models = length(ys))
  })
}
lag_selection <- common_lag(us)
p_base <- lag_selection$p[which.min(lag_selection$Total_BIC)]
save_result("Lag_Selection_BIC", lag_selection)
save_result("Chosen_Lag", data.frame(Selected_p = p_base))

desc_vars <- c(outcomes, "CFNAI", "INF_YOY", "RISK", "DREG", "DREG_L1",
               "PC1_HF", "MP_MED", "CBI_MED", "MP_PM", "CBI_PM")
save_result("US_Descriptive", purrr::map_dfr(desc_vars, function(v) {
  x <- us[[v]]
  data.frame(Variable = v, N = sum(!is.na(x)), Mean = mean(x, na.rm = TRUE),
    SD = sd(x, na.rm = TRUE), P90 = as.numeric(quantile(x, .9, na.rm = TRUE)),
    Max = max(x, na.rm = TRUE), Skewness = moments::skewness(x, na.rm = TRUE),
    Zero_Share = mean(abs(x) < 1e-12, na.rm = TRUE))
}))
save_result("Shock_Descriptive", purrr::map_dfr(c("PC1_HF", "MP_MED", "CBI_MED",
                                                 "MP_PM", "CBI_PM"), function(v) {
  x <- us[[v]]
  data.frame(Shock = v, N = sum(!is.na(x)), Mean = mean(x, na.rm = TRUE),
    SD = sd(x, na.rm = TRUE), Min = min(x, na.rm = TRUE), Max = max(x, na.rm = TRUE),
    Nonzero_Share = mean(abs(x) > 1e-12, na.rm = TRUE))
}))
cor_vars <- c(shocks, "PC1_HF_Z", controls[1:2], "RISK_L1_Z", state)
cor_matrix <- cor(us[, cor_vars], use = "pairwise.complete.obs")
save_result("US_Correlations", data.frame(Variable = rownames(cor_matrix), cor_matrix,
                                           row.names = NULL, check.names = FALSE))
unit_vars <- c(outcomes, "CFNAI", "INF_YOY", "RISK", "DREG", "MP_MED",
               "CBI_MED", "MP_PM", "CBI_PM", "PC1_HF")
unit_test <- function(x, v) {
  x <- na.omit(as.numeric(x))
  get_test <- function(fun) tryCatch(fun(x), error = function(e) NULL)
  adf <- get_test(function(z) tseries::adf.test(z, k = min(12L, floor((length(z) - 1)^(1/3)))))
  pp <- get_test(tseries::pp.test)
  kpss <- get_test(function(z) tseries::kpss.test(z, null = "Level"))
  get <- function(test, field) if (is.null(test)) NA_real_ else as.numeric(test[[field]])
  data.frame(Variable = v, ADF_Statistic = get(adf, "statistic"), ADF_p = get(adf, "p.value"),
    PP_Statistic = get(pp, "statistic"), PP_p = get(pp, "p.value"),
    KPSS_Statistic = get(kpss, "statistic"), KPSS_p = get(kpss, "p.value"))
}
save_result("Unit_Root_Tests", purrr::map_dfr(unit_vars, function(v) unit_test(us[[v]], v)))
vif_fit <- lm(LPI ~ MP_MED_Z + CBI_MED_Z + CFNAI_L1_Z + INF_L1_Z +
              DREG_L1_Z + MP_MED_Z:DREG_L1_Z + CBI_MED_Z:DREG_L1_Z, us)
vif_values <- car::vif(vif_fit)
save_result("VIF_Baseline", if (is.matrix(vif_values))
  data.frame(Term = rownames(vif_values), vif_values, row.names = NULL) else
  data.frame(Term = names(vif_values), VIF = as.numeric(vif_values)))
for (kind in c("ACF", "PACF")) save_result(paste0(kind, "_Bubbles"),
  purrr::map_dfr(outcomes, function(v) {
    x <- na.omit(us[[v]])
    obj <- if (kind == "ACF") acf(x, lag.max = 24L, plot = FALSE) else
      pacf(x, lag.max = 24L, plot = FALSE)
    data.frame(Outcome = v, Lag = as.numeric(obj$lag),
               Value = as.numeric(obj$acf), Bound = 1.96 / sqrt(length(x)))
  }))
save_result("US_Linear_LP_H12", purrr::map_dfr(outcomes, linear_path, d = us, p = p_base))
save_result("US_Linear_LP_H24", purrr::map_dfr(outcomes, linear_path, d = us,
                                               p = p_base, max_h = 24L))
save_result("Cumulative_Horizon_Tests", purrr::map_dfr(c(3L, 6L, 12L), function(h)
  purrr::map_dfr(outcomes, function(y) {
    obj <- fit_lp(us, y, p = p_base, h = h, cumulative = TRUE)
    cbind(coefs(obj, shocks), outcome = y, cumulative_horizon = h, N = nobs(obj$fit))
  })))

interaction_path <- function(y, d, s = shocks, x = controls, p = p_base,
                             group = NULL, max_h = H) {
  q <- quantile(d[[state]], c(.1, .25, .5, .75, .9), na.rm = TRUE)
  purrr::map_dfr(0:max_h, function(h) {
    obj <- fit_lp(d, y, s, x, p, h, group, interaction = TRUE)
    b <- coef(obj$fit); V <- obj$V
    purrr::map_dfr(s, function(v) purrr::map_dfr(seq_along(q), function(i) {
      nm <- paste0("INT_", v); a <- as.numeric(q[i])
      est <- b[v] + a * b[nm]
      se <- sqrt(max(V[v, v] + a^2 * V[nm, nm] + 2 * a * V[v, nm], 0))
      data.frame(outcome = y, h = h, shock = v, State_Quantile = names(q)[i],
        State_Value = a, estimate = as.numeric(est), se = se,
        lo95 = as.numeric(est - 1.96 * se), hi95 = as.numeric(est + 1.96 * se),
        p = 2 * pnorm(abs(est / se), lower.tail = FALSE))
    }))
  })
}
save_result("Interaction_Marginal", purrr::map_dfr(outcomes, interaction_path, d = us))

smooth_vc <- function(d, y, h, p = p_base, s = shocks, n_grid = 60L) {
  z <- dynamic(d, y, s, p)
  dd <- z$data
  dd$Y <- dplyr::lead(dd[[y]], h)
  rhs <- unique(c(s, state, "CFNAI_L1_Z", "INF_L1_Z", z$terms))
  dd <- dd[complete.cases(dd[, c("Y", rhs), drop = FALSE]), , drop = FALSE]
  f <- as.formula(paste("Y ~", paste(c(s, paste0("s(", state, ",bs='cr',k=6)"),
     paste0("s(", state, ",by=", s, ",bs='cr',k=6)"),
     "CFNAI_L1_Z", "INF_L1_Z", z$terms), collapse = "+")))
  fit <- mgcv::gam(f, data = dd, method = "REML", select = TRUE)
  grid <- seq(quantile(dd[[state]], .05), quantile(dd[[state]], .95), length.out = n_grid)
  purrr::map_dfr(s, function(v) purrr::map_dfr(grid, function(a) {
    base <- as.data.frame(as.list(stats::setNames(rep(0, length(rhs)), rhs)))
    base[[state]] <- a; raised <- base; raised[[v]] <- 1
    L <- predict(fit, raised, type = "lpmatrix") - predict(fit, base, type = "lpmatrix")
    est <- as.numeric(L %*% coef(fit))
    se <- sqrt(max(as.numeric(L %*% fit$Vp %*% t(L)), 0))
    data.frame(outcome = y, h = h, shock = v, State_Value = a,
               estimate = est, se = se, lo95 = est - 1.96 * se, hi95 = est + 1.96 * se)
  }))
}
save_result("Smooth_Varying_Coefficient", purrr::map_dfr(outcomes, function(y)
  purrr::map_dfr(c(0L, 3L, 6L, 12L), function(h) smooth_vc(us, y, h))))

threshold_profile <- function(d, p = p_base, h = 1L, group = NULL,
                              grid_n = 61L, ys = outcomes) {
  grid <- seq(quantile(d[[state]], .15, na.rm = TRUE),
              quantile(d[[state]], .85, na.rm = TRUE), length.out = grid_n)
  purrr::map_dfr(grid, function(cut) {
    rss <- 0; nn <- 0L
    for (y in ys) {
      z <- dynamic(d, y, shocks, p, group); dd <- z$data
      dd$Y <- lag_by(dd[[y]], h, z$group, TRUE)
      dd$HIGH <- as.numeric(dd[[state]] > cut)
      ints <- paste0(shocks, "_HIGH")
      for (i in seq_along(shocks)) dd[[ints[i]]] <- dd[[shocks[i]]] * dd$HIGH
      rhs <- unique(c(shocks, "HIGH", ints, controls, z$terms))
      dd <- dd[complete.cases(dd[, c("Y", rhs), drop = FALSE]), , drop = FALSE]
      if (nrow(dd) < 50L) next
      fit <- lm(reformulate(rhs, "Y"), data = dd)
      rss <- rss + sum(residuals(fit)^2); nn <- nn + nobs(fit)
    }
    data.frame(Threshold = cut, SSR = rss, Total_N = nn)
  })
}
profile <- threshold_profile(us)
cutoff <- profile$Threshold[which.min(profile$SSR)]
save_result("Threshold_Profile", profile)
save_result("Threshold_Selected", data.frame(Threshold_Z = cutoff))

moving_blocks <- function(d, len = block_length) {
  n <- nrow(d); starts <- sample(seq_len(n - len + 1L), ceiling(n / len) + 2L, TRUE)
  pieces <- lapply(seq_along(starts), function(i) {
    z <- d[starts[i] + seq_len(len) - 1L, , drop = FALSE]
    z$.block_id <- i; z
  })
  dplyr::bind_rows(pieces)[seq_len(n), , drop = FALSE]
}
set.seed(20260921)
boot_cut <- replicate(bootstrap_reps, {
  pr <- try(threshold_profile(moving_blocks(us), group = ".block_id", grid_n = 41L),
            silent = TRUE)
  if (inherits(pr, "try-error")) NA_real_ else pr$Threshold[which.min(pr$SSR)]
})
save_result("Threshold_Bootstrap", data.frame(Bootstrap_Threshold = boot_cut))
ci_cut <- quantile(boot_cut, c(.025, .5, .975), na.rm = TRUE)
save_result("Threshold_CI", data.frame(Lower_95 = unname(ci_cut[1]),
                                        Median = unname(ci_cut[2]), Upper_95 = unname(ci_cut[3])))

transition <- function(d, cut, gamma) {
  d$FSTATE <- plogis(gamma * (d[[state]] - cut))
  terms <- "FSTATE"
  for (v in shocks) {
    lo <- paste0(v, "_LOW"); hi <- paste0(v, "_HIGH")
    d[[lo]] <- d[[v]] * (1 - d$FSTATE)
    d[[hi]] <- d[[v]] * d$FSTATE
    terms <- c(terms, lo, hi)
  }
  list(data = d, terms = terms)
}
gamma_profile <- purrr::map_dfr(c(.25, .5, 1, 2, 5, 10), function(gam) {
  rss <- sum(vapply(outcomes, function(y) {
    z <- dynamic(us, y, shocks, p_base); d <- z$data
    d$Y <- dplyr::lead(d[[y]], 1)
    t <- transition(d, cutoff, gam); d <- t$data
    rhs <- unique(c(t$terms, controls, z$terms))
    d <- d[complete.cases(d[, c("Y", rhs), drop = FALSE]), , drop = FALSE]
    sum(residuals(lm(reformulate(rhs, "Y"), d))^2)
  }, 0.0))
  data.frame(Gamma = gam, SSR = rss)
})
gamma_star <- gamma_profile$Gamma[which.min(gamma_profile$SSR)]
save_result("Gamma_Profile", gamma_profile)
save_result("Gamma_Selected", data.frame(Gamma = gamma_star))
transition_settings <- list(c = cutoff, gamma = gamma_star)

state_contrasts <- function(obj, s = shocks) purrr::map_dfr(s, function(v) {
  b <- coef(obj$fit); V <- obj$V
  lo <- paste0(v, "_LOW"); hi <- paste0(v, "_HIGH")
  a <- c(LOW = unname(b[lo]), HIGH = unname(b[hi]),
         HIGH_MINUS_LOW = unname(b[hi] - b[lo]))
  se <- c(LOW = sqrt(V[lo, lo]), HIGH = sqrt(V[hi, hi]),
          HIGH_MINUS_LOW = sqrt(max(V[hi, hi] + V[lo, lo] - 2 * V[hi, lo], 0)))
  data.frame(shock = v, regime = names(a), estimate = as.numeric(a), se = as.numeric(se),
    lo95 = as.numeric(a - 1.96 * se), hi95 = as.numeric(a + 1.96 * se),
    p = 2 * pnorm(abs(a / se), lower.tail = FALSE))
})
save_result("Smooth_Transition_LP", purrr::map_dfr(outcomes, function(y)
  purrr::map_dfr(0:H, function(h) {
    obj <- fit_lp(us, y, p = p_base, h = h, transition = transition_settings)
    cbind(state_contrasts(obj), outcome = y, h = h, N = nobs(obj$fit))
  })))

expected_sign <- function(y, s) {
  if (s %in% c("MP_MED_Z", "MP_PM_Z")) return(if (y %in% positive) -1 else 1)
  if (y %in% positive) 1 else -1
}
test_metrics <- function(est, se, df, direction) {
  t <- est / se
  c(t = t, p2 = 2 * pt(abs(t), df, lower.tail = FALSE),
    p1 = if (direction > 0) pt(t, df, lower.tail = FALSE) else pt(t, df))
}
study_test <- function(d, y, window, type = c("direct", "state"), p = p_base) {
  type <- match.arg(type)
  obj <- fit_lp(d, y, p = p, h = window, cumulative = TRUE,
                transition = if (type == "state") transition_settings else NULL)
  b <- coef(obj$fit); V <- obj$V
  purrr::map_dfr(shocks, function(v) {
    if (type == "direct") {
      est <- b[v]; se <- sqrt(V[v, v]); low <- high <- NA_real_
      direction <- expected_sign(y, v)
    } else {
      lo <- paste0(v, "_LOW"); hi <- paste0(v, "_HIGH")
      low <- b[lo]; high <- b[hi]; est <- high - low
      se <- sqrt(max(V[hi, hi] + V[lo, lo] - 2 * V[hi, lo], 0))
      direction <- if (y %in% positive) 1 else -1
    }
    z <- test_metrics(est, se, obj$fit$df.residual, direction)
    data.frame(Outcome = y, Window = paste0("0-", window), Horizon_End = window,
      Test_Type = paste(if (v == shocks[1]) "MP" else "CBI",
                        if (type == "direct") "Direct" else "High-Low"),
      Shock = v, Low_Effect = as.numeric(low), High_Effect = as.numeric(high),
      Estimate = as.numeric(est), Average_Effect = as.numeric(est / (window + 1)),
      SE = as.numeric(se), Lower_95 = as.numeric(est - 1.96 * se),
      Upper_95 = as.numeric(est + 1.96 * se), Expected_Sign = direction,
      t_Statistic = unname(z["t"]), Two_Sided_p = unname(z["p2"]),
      One_Sided_p = unname(z["p1"]), Direction_Match = sign(est) == direction,
      N = nobs(obj$fit))
  })
}
hyp <- purrr::map_dfr(c("direct", "state"), function(kind)
  purrr::map_dfr(outcomes, function(y)
    purrr::map_dfr(c(3L, 6L, 12L), function(h) study_test(us, y, h, kind))))
save_result("Hypothesis_Tests", hyp)
save_result("Hypothesis_Summary", hyp |>
  dplyr::group_by(Test_Type, Window) |>
  dplyr::summarise(Tests = dplyr::n(), Expected_Direction = sum(Direction_Match),
    Directional_Share = mean(Direction_Match),
    Supported_5pct = sum(Direction_Match & One_Sided_p < .05),
    Supported_10pct = sum(Direction_Match & One_Sided_p < .1), .groups = "drop"))

diagnostics <- purrr::map_dfr(c(1L, 3L, 6L), function(p)
  purrr::map_dfr(outcomes, function(y)
    purrr::map_dfr(c(0L, 3L, 6L, 12L), function(h) {
      fit <- fit_lp(us, y, p = p, h = h)$fit
      bg <- lmtest::bgtest(fit, order = 12)
      lb <- Box.test(residuals(fit), lag = 12, type = "Ljung-Box", fitdf = 0)
      bp <- lmtest::bptest(fit)
      cusum <- try({ strucchange::sctest(strucchange::efp(formula(fit),
                   data = model.frame(fit), type = "Rec-CUSUM"))$p.value }, silent = TRUE)
      data.frame(Outcome = y, Lag_Order = p, Horizon = h, N = nobs(fit),
        BG_p = bg$p.value, LjungBox_p = lb$p.value, BP_p = bp$p.value,
        CUSUM_p = if (inherits(cusum, "try-error")) NA_real_ else cusum)
    })))
save_result("Diagnostics_Full", diagnostics)
diag0 <- diagnostics |>
  dplyr::filter(Horizon == 0) |>
  tidyr::pivot_longer(c(BG_p, LjungBox_p), names_to = "Diagnostic_Test",
                       values_to = "p_value") |>
  dplyr::group_by(Lag_Order, Diagnostic_Test) |>
  dplyr::summarise(Rejections = sum(p_value < .05), Total = dplyr::n(),
    Rejection_Rate = mean(p_value < .05), Median_p = median(p_value), .groups = "drop")
save_result("Diagnostics_h0", diag0)

wild_band <- function(d, y, v, B = bootstrap_reps, len = block_length) {
  d$.row_id <- seq_len(nrow(d))
  models <- lapply(0:H, function(h) {
    obj <- fit_lp(d, y, p = p_base, h = h)
    leverage <- pmin(hatvalues(obj$fit), .999999)
    obj$resid_hc2 <- residuals(obj$fit) / sqrt(1 - leverage)
    obj
  })
  orig <- vapply(models, function(m) coef(m$fit)[v], 0.0)
  se0 <- vapply(models, function(m) sqrt(m$V[v, v]), 0.0)
  tstar <- matrix(NA_real_, B, H + 1L)
  for (b in seq_len(B)) {
    signs <- sample(c(-1, 1), ceiling(nrow(d) / len), replace = TRUE)
    for (j in seq_along(models)) {
      m <- models[[j]]; dd <- m$data
      dd$Ystar <- fitted(m$fit) + m$resid_hc2 * signs[ceiling(dd$.row_id / len)]
      star <- lm(reformulate(m$rhs, "Ystar"), data = dd)
      V <- sandwich::NeweyWest(star, lag = max(j, p_base), prewhite = FALSE, adjust = TRUE)
      tstar[b, j] <- (coef(star)[v] - orig[j]) / sqrt(V[v, v])
    }
  }
  max_t <- apply(abs(tstar), 1, max)
  critical <- unname(quantile(max_t, .95))
  list(bands = data.frame(outcome = y, shock = v, h = 0:H,
       estimate = orig, se = se0, lo95_pointwise = orig - 1.96 * se0,
       hi95_pointwise = orig + 1.96 * se0, Sim_Critical_95 = critical,
       Sim_Lower_95 = orig - critical * se0, Sim_Upper_95 = orig + critical * se0,
       Simultaneous_p = vapply(abs(orig / se0), function(t) mean(max_t >= t), 0.0)),
       critical = data.frame(Outcome = y, Shock = v, Bootstrap_Draws = B,
                             Block_Length = len, Critical_95 = critical))
}
set.seed(20260921)
wild <- lapply(outcomes, function(y) lapply(shocks, function(v) wild_band(us, y, v)))
wild <- purrr::flatten(wild)
save_result("Wild_Sim_Bands", dplyr::bind_rows(lapply(wild, `[[`, "bands")))
save_result("Wild_Critical", dplyr::bind_rows(lapply(wild, `[[`, "critical")))

# Four pure linear specifications and formal within-equation CBI–MP contrasts.
linear_validation <- purrr::map_dfr(c(6L, 1L), function(p)
  purrr::map_dfr(c(TRUE, FALSE), function(with_state) {
    x <- if (with_state) controls else setdiff(controls, state)
    spec <- if (with_state) "With lagged deregulation control" else "Without deregulation control"
    purrr::map_dfr(outcomes, function(y) {
      result <- linear_path(y, us, x = x, p = p)
      result$Lag_Order <- p; result$Control_Spec <- spec
      result
    })
  }))
linear_validation$Expected_Sign <- mapply(expected_sign, linear_validation$outcome,
                                          linear_validation$term)
linear_validation$Direction_Match <- sign(linear_validation$estimate) ==
                                      linear_validation$Expected_Sign
save_result("Linear_Validation_IRF", linear_validation)
validation_test <- purrr::map_dfr(c(6L, 1L), function(p)
  purrr::map_dfr(c(TRUE, FALSE), function(with_state) {
    x <- if (with_state) controls else setdiff(controls, state)
    spec <- if (with_state) "With lagged deregulation control" else "Without deregulation control"
    purrr::map_dfr(outcomes, function(y) purrr::map_dfr(c(3L, 6L, 12L), function(h) {
      m <- fit_lp(us, y, x = x, p = p, h = h, cumulative = TRUE)
      b <- coef(m$fit); V <- m$V
      sm <- expected_sign(y, shocks[1]); sc <- expected_sign(y, shocks[2])
      aligned_mp <- sm * b[shocks[1]]; aligned_cbi <- sc * b[shocks[2]]
      diff <- aligned_cbi - aligned_mp
      a <- c(-sm, sc); vv <- V[shocks, shocks, drop = FALSE]
      se <- sqrt(as.numeric(t(a) %*% vv %*% a))
      p1 <- pt(diff / se, m$fit$df.residual, lower.tail = FALSE)
      data.frame(Outcome = y, Lag_Order = p, Control_Spec = spec,
        Window = paste0("0-", h), Horizon_End = h,
        MP_Aligned_Average = as.numeric(aligned_mp / (h + 1)),
        CBI_Aligned_Average = as.numeric(aligned_cbi / (h + 1)),
        CBI_minus_MP_Average = as.numeric(diff / (h + 1)),
        SE_Difference = se, p_OneSided_CBI_Greater_MP = p1,
        CBI_Greater_MP = diff > 0,
        CBI_Greater_MP_5pct = diff > 0 && p1 < .05,
        CBI_Greater_MP_10pct = diff > 0 && p1 < .1,
        Absolute_Ratio_CBI_to_MP = if (abs(b[shocks[1]]) > 1e-10)
          abs(b[shocks[2]]) / abs(b[shocks[1]]) else NA_real_, N = nobs(m$fit))
    }))
  }))
save_result("Linear_Validation_Tests", validation_test)
save_result("Linear_Validation_Summary", validation_test |>
  dplyr::group_by(Lag_Order, Control_Spec, Window) |>
  dplyr::summarise(Number_of_Bubble_Indicators = dplyr::n(),
    CBI_Aligned_Greater_MP = sum(CBI_Greater_MP),
    Share_CBI_Greater_MP = mean(CBI_Greater_MP),
    CBI_Greater_MP_at_5pct = sum(CBI_Greater_MP_5pct),
    Median_CBI_to_MP_Absolute_Ratio = median(Absolute_Ratio_CBI_to_MP, na.rm = TRUE),
    .groups = "drop"))
save_result("Linear_Directional_Overall", linear_validation |>
  dplyr::group_by(Lag_Order, Control_Spec, term) |>
  dplyr::summarise(Total_Horizon_Outcome_Cells = dplyr::n(),
    Correct_Direction_Cells = sum(Direction_Match),
    Overall_Directional_Share = mean(Direction_Match), .groups = "drop"))

# Other reported sensitivity analyses. Calendar months remain contiguous.
save_result("Risk_Aversion_Robustness", purrr::map_dfr(outcomes, linear_path,
  d = us, x = c(controls, "RISK_L1_Z"), p = p_base))
save_result("Contemporaneous_Controls", purrr::map_dfr(outcomes, linear_path,
  d = us, x = c("CFNAI_Z", "INF_YOY_Z", "DREG_Z"), p = p_base))
save_result("Alternative_Lag_Orders", purrr::map_dfr(sort(unique(c(3L, 6L, 12L, p_base))),
  function(p) { z <- purrr::map_dfr(outcomes, linear_path, d = us, p = p)
                z$Lag_Order <- p; z }))
winsor <- function(x) { q <- quantile(x, c(.01, .99), na.rm = TRUE); pmax(q[1], pmin(q[2], x)) }
us_win <- us
for (v in shocks) us_win[[v]] <- winsor(us_win[[v]])
save_result("Winsorized_Shocks", purrr::map_dfr(outcomes, linear_path, d = us_win, p = p_base))
save_result("Alternative_PM_LP", purrr::map_dfr(outcomes, linear_path,
  d = us, s = c("MP_PM_Z", "CBI_PM_Z"), p = p_base))
save_result("Undecomposed_PC1_LP", purrr::map_dfr(outcomes, linear_path,
  d = us, s = "PC1_HF_Z", p = p_base))
save_result("U1_U4_Linear_LP", purrr::map_dfr(outcomes, linear_path,
  d = us_u, s = paste0("U", 1:4, "_Z"), p = p_base))
save_result("U1_U4_Marginal", purrr::map_dfr(outcomes, interaction_path,
  d = us_u, s = paste0("U", 1:4, "_Z"), p = p_base))
us_binary <- us
for (v in outcomes) us_binary[[paste0(v, "_BIN")]] <- as.numeric(us_binary[[v]] > 0)
save_result("Zero_Extensive_Margin", purrr::map_dfr(paste0(outcomes, "_BIN"),
  linear_path, d = us_binary, p = p_base))
save_result("Zero_Intensive_Margin", purrr::map_dfr(outcomes, function(y)
  purrr::map_dfr(0:H, function(h) {
    z <- dynamic(us, y, shocks, p_base); d <- z$data
    d$Y <- dplyr::lead(d[[y]], h)
    rhs <- unique(c(shocks, controls, z$terms))
    d <- d[complete.cases(d[, c("Y", rhs), drop = FALSE]) & d$Y > 0, , drop = FALSE]
    if (nrow(d) < 35L) return(NULL)
    fit <- lm(reformulate(rhs, "Y"), d)
    obj <- list(fit = fit, V = sandwich::NeweyWest(fit,
      lag = max(h + 1L, p_base), prewhite = FALSE, adjust = TRUE))
    cbind(coefs(obj, shocks), outcome = y, h = h, N = nobs(fit))
  })))

eu_lags <- purrr::map_dfr(c("FR", "DE", "IT"), function(iso) {
  d <- eu[eu$ISO == iso, , drop = FALSE]
  z <- common_lag(d, x = character(), pmax = 6L)
  z$ISO <- iso; z$Selected <- z$p == z$p[which.min(z$Total_BIC)]
  z
})
save_result("EU_Lag_Selection", eu_lags)
eu_linear <- eu_int <- eu_pm <- eu_pc1 <- list()
for (iso in c("FR", "DE", "IT")) {
  d <- eu[eu$ISO == iso, , drop = FALSE]
  p <- eu_lags$p[eu_lags$ISO == iso & eu_lags$Selected][1]
  label <- unique(d$COUNTRY)[1]
  add_id <- function(z) { z$ISO <- iso; z$COUNTRY <- label; z$Lag_Order <- p; z }
  eu_linear[[iso]] <- add_id(purrr::map_dfr(outcomes, linear_path, d = d,
                                            x = character(), p = p))
  eu_int[[iso]] <- add_id(purrr::map_dfr(outcomes, interaction_path, d = d,
                                         x = character(), p = p))
  eu_pm[[iso]] <- add_id(purrr::map_dfr(outcomes, linear_path, d = d,
    s = c("MP_PM_Z", "CBI_PM_Z"), x = character(), p = p))
  eu_pc1[[iso]] <- add_id(purrr::map_dfr(outcomes, linear_path, d = d,
    s = "PC1_HF_Z", x = character(), p = p))
}
save_result("EU_Linear_LP", dplyr::bind_rows(eu_linear))
save_result("EU_Interaction", dplyr::bind_rows(eu_int))
save_result("EU_PM_Robustness", dplyr::bind_rows(eu_pm))
save_result("EU_PC1_Robustness", dplyr::bind_rows(eu_pc1))

rolling <- purrr::map_dfr(c("MPI", "LPI"), function(y)
  purrr::map_dfr(shocks, function(v) {
    starts <- seq(1L, nrow(us) - 180L + 1L, by = 12L)
    purrr::map_dfr(starts, function(i) {
      d <- us[i:(i + 179L), , drop = FALSE]
      obj <- fit_lp(d, y, p = min(p_base, 6L), h = 6L)
      tab <- coefs(obj, v)
      data.frame(outcome = y, shock = v, Start_Date = min(d$DATE),
        End_Date = max(d$DATE), estimate = tab$estimate, se = tab$se,
        lo95 = tab$lo95, hi95 = tab$hi95)
    })
  }))
save_result("Rolling_Stability", rolling)

# Publication plots are computed from the saved coefficients, never entered by hand.
plot_irf <- function(z, name, shock_col = "term", facet = "outcome") {
  p <- ggplot2::ggplot(z, ggplot2::aes(x = h, y = estimate,
       colour = .data[[shock_col]], fill = .data[[shock_col]])) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo95, ymax = hi95),
                         alpha = .12, colour = NA) +
    ggplot2::geom_line(linewidth = .8) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::facet_wrap(stats::as.formula(paste("~", facet)), scales = "free_y") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(x = "Months after announcement", y = "Confidence-score response")
  ggplot2::ggsave(file.path("results", "figures", paste0(name, ".png")),
                  p, width = 12, height = 9, dpi = 600)
}
plot_irf(outputs$US_Linear_LP_H12, "linear_p1")
plot_irf(outputs$Smooth_Transition_LP[
  outputs$Smooth_Transition_LP$regime != "HIGH_MINUS_LOW", ], "smooth_transition",
  shock_col = "regime")
for (s in shocks) {
  z <- outputs$Interaction_Marginal |>
    dplyr::filter(shock == s, State_Quantile %in% c("10%", "50%", "90%"))
  plot_irf(z, paste0("interaction_", s), shock_col = "State_Quantile")
}
ggplot2::ggsave(file.path("results", "figures", "threshold_profile.png"),
  ggplot2::ggplot(profile, ggplot2::aes(Threshold, SSR)) +
    ggplot2::geom_line() + ggplot2::geom_vline(xintercept = cutoff, linetype = 2) +
    ggplot2::theme_minimal(base_size = 13), width = 10, height = 6, dpi = 600)

settings <- data.frame(setting = c("lag", "threshold", "gamma", "bootstrap_draws",
  "block_length", "US_first_date", "US_last_date", "input_file"),
  value = as.character(c(p_base, cutoff, gamma_star, bootstrap_reps, block_length,
                         min(us$DATE), max(us$DATE), input)))
save_result("Settings", settings)
openxlsx::write.xlsx(outputs, file.path("results", "all_results.xlsx"),
                     overwrite = TRUE)
saveRDS(outputs, file.path("results", "all_results.rds"))
writeLines(capture.output(sessionInfo()), file.path("results", "session_info.txt"))
cat("Complete. Results are in results/.\n")
