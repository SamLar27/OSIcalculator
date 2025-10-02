#' OSI calculator (hockey-stick or linear) from a single FEV1/Tif pair
#'
#' Provide a single pair of values (FEV1_PCT, Tif_PCT).
#'   Hockey-stick:
#'     - Pre-BD threshold = 70
#'     - Post-BD threshold = 85
#'   Linear:
#'     - Uses pooled pre-BD or post-BD coefficients (per your MI + Rubin pooling)
#'
#' @param FEV1_PCT numeric vector, FEV1 % predicted (NA allowed)
#' @param Tif_PCT  numeric vector, FEV1/FVC % predicted (NA allowed)
#' @param type     "Pre_BD", "Post_BD", or "Pooled" (default: "Pre_BD")
#' @param model    "hockey_stick" (default) or "linear"
#'
#' @return numeric vector of OSI values (for "Pooled": mean of pre & post).
#'         Returns NA for any input row containing NA values.
#' @export
OSI_calculator <- function(FEV1_PCT, Tif_PCT,
                           type  = c("Pre_BD", "Post_BD", "Pooled"),
                           model = c("hockey_stick", "linear")) {
  type  <- match.arg(type)
  model <- match.arg(model)

  ## ---- validation ----
  if (!is.numeric(FEV1_PCT) || !is.numeric(Tif_PCT)) {
    stop("FEV1_PCT and Tif_PCT must be numeric.", call. = FALSE)
  }
  if (length(FEV1_PCT) != length(Tif_PCT)) {
    stop("FEV1_PCT and Tif_PCT must have the same length.", call. = FALSE)
  }

  # Track NAs to propagate
  na_mask <- is.na(FEV1_PCT) | is.na(Tif_PCT)

  # Shared offset (rates per 365 days)
  log_offset <- log(365)

  ## ---- helpers ----
  compute_hockey <- function(FEV1, Tif, coef, covars, threshold, hs_name) {
    # constant part: intercept + offset + covariates
    lp_const <- coef[["(Intercept)"]] + log_offset
    for (nm in names(covars)) {
      if (nm %in% names(coef)) lp_const <- lp_const + coef[[nm]] * covars[[nm]]
    }

    # Tif coefficient (name differs between pre/post)
    tif_name <- names(coef)[grep("^Tif_.*_PCT_0W$", names(coef))]
    if (length(tif_name) != 1) {
      stop("Could not uniquely identify Tif coefficient name.", call. = FALSE)
    }

    # "Normal" predictor at 100/100
    hs_norm <- max(0, threshold - 100)
    lp_norm <- lp_const + coef[[hs_name]] * hs_norm + coef[[tif_name]] * 100
    risk_norm <- exp(lp_norm)

    # Real predictor with provided FEV1/Tif
    hs_real <- pmax(0, threshold - FEV1, na.rm = FALSE)
    lp_real <- lp_const + coef[[hs_name]] * hs_real + coef[[tif_name]] * Tif
    risk_real <- exp(lp_real)

    out <- as.numeric(risk_real / risk_norm)
    out[na_mask] <- NA_real_
    out
  }

  compute_linear <- function(FEV1, Tif, coef, covars, fev1_pat, tif_pat) {
    # constant part: intercept + offset + covariates
    lp_const <- coef[["(Intercept)"]] + log_offset
    for (nm in names(covars)) {
      if (nm %in% names(coef)) lp_const <- lp_const + coef[[nm]] * covars[[nm]]
    }

    # Locate the FEV1 and Tif coefficient names robustly
    fev1_name <- names(coef)[grep(fev1_pat, names(coef))]
    tif_name  <- names(coef)[grep(tif_pat,  names(coef))]
    if (length(fev1_name) != 1 || length(tif_name) != 1) {
      stop("Could not uniquely identify FEV1/Tif coefficient names for linear model.", call. = FALSE)
    }

    # normal (100/100)
    lp_norm <- lp_const + coef[[fev1_name]] * 100 + coef[[tif_name]] * 100
    risk_norm <- exp(lp_norm)

    # real values
    lp_real <- lp_const + coef[[fev1_name]] * FEV1 + coef[[tif_name]] * Tif
    risk_real <- exp(lp_real)

    out <- as.numeric(risk_real / risk_norm)
    out[na_mask] <- NA_real_
    out
  }

  ## ---- hockey-stick parameters (unchanged) ----
  if (model == "hockey_stick") {
    # Pre-BD (threshold 70)
    coef_pre <- c(
      "(Intercept)"               = -11.209318534,
      "hs_FEV1_70"                =   0.008718600,
      "Tif_preBD_PCT_0W"          =  -0.006450547,
      "SexMale"                   =  -0.195311773,
      "BMI"                       =   0.011065735,
      "Smoking_history_yes_noYes" =   0.305937906,
      "CRSCRSsNP"                 =   0.186434571,
      "CRSCRSwNP"                 =   0.198466095,
      "Allergic_RhinitisYes"      =   0.083795510,
      "GINA_step_numeric"         =   0.197548736,
      "ACQ_score_0W"              =   0.116025550,
      "Attack_history"            =   0.227946799,
      "BEC_log10"                 =   0.343092398,
      "FeNO_log10"                =   0.368922673
    )
    covars_pre <- list(
      SexMale = 0,
      BMI = 28.1,
      Smoking_history_yes_noYes = 0,
      CRSCRSsNP = 0,
      CRSCRSwNP = 0,
      Allergic_RhinitisYes = 0,
      GINA_step_numeric = 4,
      ACQ_score_0W = 2.4,
      Attack_history = 2,
      BEC_log10 = 8.397940,
      FeNO_log10 = 1.332438
    )

    # Post-BD (threshold 85)
    coef_post <- c(
      "(Intercept)"               = -11.156455138,
      "hs_FEV1_postBD_85"         =   0.006919494,
      "Tif_postBD_PCT_0W"         =  -0.005362333,
      "SexMale"                   =  -0.194883278,
      "BMI"                       =   0.008282615,
      "Smoking_history_yes_noYes" =   0.232018557,
      "CRSCRSsNP"                 =   0.217138590,
      "CRSCRSwNP"                 =   0.190478519,
      "Allergic_RhinitisYes"      =   0.075043557,
      "GINA_step_numeric"         =   0.327958915,
      "ACQ_score_0W"              =   0.112912777,
      "Attack_history"            =   0.215231020,
      "BEC_log10"                 =   0.296815830,
      "FeNO_log10"                =   0.224753254
    )
    covars_post <- list(
      SexMale = 0,
      BMI = 28.00000000,
      Smoking_history_yes_noYes = 0,
      CRSCRSsNP = 0,
      CRSCRSwNP = 0,
      Allergic_RhinitisYes = 1,
      GINA_step_numeric = 4.00000000,
      ACQ_score_0W = 2.40000000,
      Attack_history = 2.00000000,
      BEC_log10 = 8.38021124,
      FeNO_log10 = 1.37106786
    )

    if (type == "Pre_BD") {
      return(compute_hockey(FEV1_PCT, Tif_PCT, coef_pre,  covars_pre,  70, "hs_FEV1_70"))
    } else if (type == "Post_BD") {
      return(compute_hockey(FEV1_PCT, Tif_PCT, coef_post, covars_post, 85, "hs_FEV1_postBD_85"))
    } else { # Pooled
      pre  <- compute_hockey(FEV1_PCT, Tif_PCT, coef_pre,  covars_pre,  70, "hs_FEV1_70")
      post <- compute_hockey(FEV1_PCT, Tif_PCT, coef_post, covars_post, 85, "hs_FEV1_postBD_85")
      out <- (pre + post) / 2
      out[na_mask] <- NA_real_
      return(out)
    }
  }

  ## ---- linear parameters (separate pre- and post-BD) ----
  # PRE-BD linear (your previously pooled pre-BD coefficients)
  coef_pre_linear <- c(
    "(Intercept)" = -10.767819412223,
    "FEV1_preBD_PCT_0W" = -0.004937991998,
    "Tif_preBD_PCT_0W"  = -0.007575109211,
    "SexMale" = -0.194812009093,
    "BMI" = 0.011337184198,
    "Smoking_history_yes_noYes" = 0.307259669380,
    "CRSCRSsNP" = 0.185496792529,
    "CRSCRSwNP" = 0.199355069785,
    "Allergic_RhinitisYes" = 0.079450744564,
    "GINA_step_numeric" = 0.199465333811,
    "ACQ_score_0W" = 0.119925331430,
    "Attack_history" = 0.228782649735,
    "BEC_log10" = 0.346510714967,
    "FeNO_log10" = 0.362382297626
  )

  # Reference covariates for PRE-BD linear (as used before; Allergic_RhinitisYes = 0)
  covars_pre_linear <- list(
    SexMale = 0,
    BMI = 28.1,
    Smoking_history_yes_noYes = 0,
    CRSCRSsNP = 0,
    CRSCRSwNP = 0,
    Allergic_RhinitisYes = 0,
    GINA_step_numeric = 4,
    ACQ_score_0W = 2.4,
    Attack_history = 2,
    BEC_log10 = 8.39794001,
    FeNO_log10 = 1.34242268
  )

  # POST-BD linear (your newly provided pooled post-BD coefficients)
  coef_post_linear <- c(
    "(Intercept)" = -10.637319415759,
    "FEV1_postBD_PCT_0W" = -0.005175221748,
    "Tif_postBD_PCT_0W"  = -0.005795153302,
    "SexMale" = -0.195146389353,
    "BMI" = 0.008243115294,
    "Smoking_history_yes_noYes" = 0.233664965443,
    "CRSCRSsNP" = 0.220200395210,
    "CRSCRSwNP" = 0.191860234696,
    "Allergic_RhinitisYes" = 0.073584384851,
    "GINA_step_numeric" = 0.326309331302,
    "ACQ_score_0W" = 0.113578910571,
    "Attack_history" = 0.215323871453,
    "BEC_log10" = 0.297454219224,
    "FeNO_log10" = 0.221367585392
  )

  # Reference covariates for POST-BD linear (mapped from your labels):
  covars_post_linear <- list(
    SexMale = 0,
    BMI = 28.00000000,
    Smoking_history_yes_noYes = 0,
    CRSCRSsNP = 0,
    CRSCRSwNP = 0,
    Allergic_RhinitisYes = 1,
    GINA_step_numeric = 4.00000000,
    ACQ_score_0W = 2.40000000,
    Attack_history = 2.00000000,
    BEC_log10 = 8.38021124,
    FeNO_log10 = 1.37106786
  )

  # Choose linear branch by type
  if (type == "Pre_BD") {
    lin <- compute_linear(
      FEV1_PCT, Tif_PCT,
      coef = coef_pre_linear,
      covars = covars_pre_linear,
      fev1_pat = "FEV1_preBD_PCT_0W",
      tif_pat  = "^Tif_preBD_PCT_0W$"
    )
    lin[na_mask] <- NA_real_
    return(lin)
  } else if (type == "Post_BD") {
    lin <- compute_linear(
      FEV1_PCT, Tif_PCT,
      coef = coef_post_linear,
      covars = covars_post_linear,
      fev1_pat = "FEV1_postBD_PCT_0W",
      tif_pat  = "^Tif_postBD_PCT_0W$"
    )
    lin[na_mask] <- NA_real_
    return(lin)
  } else { # Pooled: average of pre- and post-BD linear predictions
    lin_pre <- compute_linear(
      FEV1_PCT, Tif_PCT,
      coef = coef_pre_linear,
      covars = covars_pre_linear,
      fev1_pat = "FEV1_preBD_PCT_0W",
      tif_pat  = "^Tif_preBD_PCT_0W$"
    )
    lin_post <- compute_linear(
      FEV1_PCT, Tif_PCT,
      coef = coef_post_linear,
      covars = covars_post_linear,
      fev1_pat = "FEV1_postBD_PCT_0W",
      tif_pat  = "^Tif_postBD_PCT_0W$"
    )
    out <- (lin_pre + lin_post) / 2
    out[na_mask] <- NA_real_
    return(out)
  }
}
