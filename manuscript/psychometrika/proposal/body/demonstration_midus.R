
# Setup ------------------------------------------------------------------------

library(tidyverse)
library(lavaan)
library(R2spa)
library(modelsummary)

midus1_name <- load(here::here("data/02760-0001-Data.rda"))
midus2_name <- load(here::here("data/04652-0001-Data.rda"))
midus3_name <- load(here::here("data/36346-0001-Data.rda"))
# Select MIDUS 1 variables (ID, age, personal growth and purpose in life items)
midus1_df <- get(midus1_name)[, c(
  "M2ID", "A1PRAGE_2019",
  "A1SF1K", "A1SF1L", "A1SF1N", # personal grwoth
  "A1SF1C", "A1SF1G", "A1SF1J" # purpose in life
)]
# Select MIDUS 2 variables (ID, personal growth and purpose in life items)
midus2_df <- get(midus2_name)[, c(
  "M2ID",
  "B1SE1AA", "B1SE1I", "B1SE1GG", 
  "B1SE1OO", "B1SE1E", "B1SE1QQ"
)]
# Select MIDUS 3 variables (ID, personal growth and purpose in life items)
midus3_df <- get(midus3_name)[, c(
  "M2ID",
  "C1SE1AA", "C1SE1I", "C1SE1GG", 
  "C1SE1OO", "C1SE1E", "C1SE1QQ"
)]
# Merge data
midus_merged <- midus1_df |> 
  left_join(midus2_df, by = "M2ID") |>
  left_join(midus3_df, by = "M2ID") |>
  rename(m2id = M2ID, age = A1PRAGE_2019,
         # personal growth
         learn1 = A1SF1K, newexp1 = A1SF1L, improve1 = A1SF1N,
         learn2 = B1SE1AA, newexp2 = B1SE1I, improve2 = B1SE1GG,
         learn3 = C1SE1AA, newexp3 = C1SE1I, improve3 = C1SE1GG, 
         # purpose of life
         purpose1 = A1SF1C, future1 = A1SF1G, done1 = A1SF1J, 
         purpose2 = B1SE1OO, future2 = B1SE1E, done2 = B1SE1QQ, 
         purpose3 = C1SE1OO, future3 = C1SE1E, done3 = C1SE1QQ) |>
  # Subset participants aged 40 or below at Wave 1
  filter(age <= 40) |>
  # Drop missing data
  drop_na() |>
  # Convert items from factor to numeric 
  mutate_at(vars(learn1:done3), as.numeric) |>
  mutate_at(vars(matches("learn|newexp|purpose")), ~ 7 - .)

# Partial strict model ---------------------------------------------------------

# Please fine more details in "exploration/midus/midus_personal_growth.qmd"

# Strict invariance model
pstrict_mod <- "
    eta1 =~ NA * learn1 + (lam1) * learn1 + (lam2) * newexp1 + (lam31) * improve1
		eta2 =~ NA * learn2 + (lam1) * learn2 + (lam2) * newexp2 + (lam3) * improve2
		eta3 =~ NA * learn3 + (lam1) * learn3 + (lam2) * newexp3 + (lam3) * improve3
		# Measurement intercepts
    learn1 ~ (nu1) * 1
    newexp1 ~ (nu2) * 1
    improve1 ~ (nu31) * 1
    learn2 ~ (nu1) * 1
    newexp2 ~ (nu2) * 1
    improve2 ~ (nu3) * 1
    learn3 ~ (nu1) * 1
    newexp3 ~ (nu2) * 1
    improve3 ~ (nu3) * 1
    # Unique factor variances
    learn1 ~~ (theta11) * learn1
    learn2 ~~ (theta1) * learn2
    learn3 ~~ (theta1) * learn3
    newexp1 ~~ (theta21) * newexp1
    newexp2 ~~ (theta2) * newexp2
    newexp3 ~~ (theta2) * newexp3
    improve1 ~~ (theta3) * improve1
    improve2 ~~ (theta3) * improve2
    improve3 ~~ (theta3) * improve3
    # Unique factor covariances
    learn1 ~~ learn2 + learn3
    learn2 ~~ learn3
    newexp1 ~~ newexp2 + newexp3
    newexp2 ~~ newexp3
    improve1 ~~ improve2 + improve3
    improve2 ~~ improve3
    # First factor variance and mean to 1 and 0
    eta1 ~~ 1 * eta1
    eta1 ~ 0
    # Other factor means to free
    eta2 ~ NA * 1
    eta3 ~ NA * 1
"
pstrict_fit <- cfa(pstrict_mod, midus_merged)
filter(parameterestimates(pstrict_fit), op == "=~")

# 2S-PA ------------------------------------------------------------------------

get_tspa_mod <- function(data, vc, fsA) {
  var <- colnames(vc)
  len <- nrow(vc)
  
  col <- colnames(data)
  fs <- paste0("fs_", var)
  colnames(vc) <- rownames(vc) <- fs
  
  # latent variables
  loadings <- paste0(fsA, " * ", fs)
  loadings_list <- split(loadings, factor(rep(var, each = len),
                                          levels = var))
  loadings_c <- lapply(loadings_list, function(x) {
    paste0(x, collapse = " + ")
  })
  latent_var_str <- paste(var, "=~", loadings_c)
  # error variances
  vc_in <- !upper.tri(vc)
  ev_rhs <- colnames(vc)[col(vc_in)[vc_in]]
  ev_lhs <- rownames(vc)[row(vc_in)[vc_in]]
  error_constraint_str <- paste0(ev_lhs, " ~~ ", vc[vc_in], " * ", ev_rhs)
  
  paste0(c(
    "# latent variables (indicated by factor scores)",
    latent_var_str,
    "# constrain the errors",
    error_constraint_str
  ),
  collpase = "\n")
}
fs_pgrowth <- get_fs(midus_merged, pstrict_mod, method = "Bartlett")
tspa_pgrowth <- get_tspa_mod(midus_merged, 
                             vc = attr(fs_pgrowth, "av_efs"), 
                             fsA = attr(fs_pgrowth, "fsA"))

tspa_mod <- "
# Latent variables (indicated by factor scores)
eta1 =~ 1 * fs_eta1
eta2 =~ 1 * fs_eta2
eta3 =~ 1 * fs_eta3

# Constrain the errors
fs_eta1 ~~ 0.270712395783511 * fs_eta1
fs_eta2 ~~ 0.0439357846685017 * fs_eta1
fs_eta3 ~~ 0.0179222578946954 * fs_eta1
fs_eta2 ~~ 0.536000083179472 * fs_eta2
fs_eta3 ~~ 0.129328073618597 * fs_eta2
fs_eta3 ~~ 0.521117152047182 * fs_eta3

# Constrain intercepts
fs_eta1 ~ 0 * 1
fs_eta2 ~ 0 * 1
fs_eta3 ~ 0 * 1

# Linear Growth Model
i =~ 1 * eta1 + 1 * eta2 + 1 * eta3
s =~ 0 * eta1 + 1 * eta2 + 2 * eta3

# Variance-covariances of intercepts and slopes
i ~~ phi1 * i
s ~~ s
i ~~ s

# Means of level and slope
i ~ 1
s ~ 1

# Fixed disturbances of latent outcomes to zero
eta1 ~ 0 * 1
eta2 ~ 0 * 1
eta3 ~ 0 * 1

# Constrain total variance at Wave 1
eta1 ~~ psi1 * eta1
psi1 + phi1 == 1
"
tspa_fit <- sem(tspa_mod, fs_pgrowth)
saveRDS(tspa_fit, "manuscript/psychometrika/proposal/body/tspa_fit.rds")

# JSEM -------------------------------------------------------------------------

jsem_mod <- "
# Factor loadings
eta1 =~ NA * learn1 + l1 * learn1 + l2 * newexp1 + l31 * improve1
eta2 =~ NA * learn2 + l1 * learn2 + l2 * newexp2 + l3 * improve2
eta3 =~ NA * learn3 + l1 * learn3 + l2 * newexp3 + l3 * improve3

# Measurement intercepts
learn1 ~ 5.341 * 1
newexp1 ~ nu2 * 1
improve1 ~ nu31 * 1
learn2 ~ 5.341 * 1
newexp2 ~ nu2 * 1
improve2 ~ nu3 * 1
learn3 ~ 5.341 * 1
newexp3 ~ nu2 * 1
improve3 ~ nu3 * 1

# Unique factor variances
learn1 ~~ (theta11) * learn1
learn2 ~~ (theta1) * learn2
learn3 ~~ (theta1) * learn3
newexp1 ~~ (theta21) * newexp1
newexp2 ~~ (theta2) * newexp2
newexp3 ~~ (theta2) * newexp3
improve1 ~~ (theta3) * improve1
improve2 ~~ (theta3) * improve2
improve3 ~~ (theta3) * improve3

# Unique factor covariances
learn1 ~~ learn2 + learn3
learn2 ~~ learn3
newexp1 ~~ newexp2 + newexp3
newexp2 ~~ newexp3
improve1 ~~ improve2 + improve3
improve2 ~~ improve3

# Linear growth model
i =~ 1 * eta1 + 1 * eta2 + 1 * eta3
s =~ 0 * eta1 + 1 * eta2 + 2 * eta3

# Variance-covariances of intercepts and slopes
i ~~ phi1 * i
s ~~ s
i ~~ s

# Means of level (fixed to zero for identification) and slope
i ~ 1 
s ~ 1

# Disturbances of latent outcomes
eta1 ~ 0 * 1
eta2 ~ 0 * 1
eta3 ~ 0 * 1

# Constrain total variance at Wave 1
eta1 ~~ psi1 * eta1
psi1 + phi1 == 1
"
jsem_fit <- sem(jsem_mod, midus_merged)
saveRDS(jsem_fit, "manuscript/psychometrika/proposal/body/jsem_fit.rds")

# Comparison -------------------------------------------------------------------

msummary(
  list(
    "JSEM Growth" = jsem_fit, 
    "2S-PA Growth" = tspa_fit
  ),
  output = "markdown",
  statistic = "conf.int",
  gof_omit = "IC", 
  coef_map = c(
    "i ~1 " = "Mean(Level)",
    "s ~1 " = "Mean(Slope)",
    "i ~~ i" = "Var(Level)",
    "s ~~ s" = "Var(Slope)",
    "i ~~ s" = "Cov(Level, Slope)"
  )
)

# Archive ----------------------------------------------------------------------
# (1) If the measurement model (e.g., constraining noninvariant items), is 
#     mispecified, 2S-PA and JSEM also yield biased results---totally reasonable
# (2) Loading estimates change in the JSEM even when we use the first loading
#     estimate from the partial strict model

# 2S-PA (misspecified) ---------------------------------------------------------

mis_pstrict_mod <- "
    eta1 =~ NA * learn1 + (lam1) * learn1 + (lam2) * newexp1 + (lam31) * improve1
		eta2 =~ NA * learn2 + (lam1) * learn2 + (lam2) * newexp2 + (lam3) * improve2
		eta3 =~ NA * learn3 + (lam1) * learn3 + (lam2) * newexp3 + (lam3) * improve3
		# Measurement intercepts
    learn1 ~ (nu1) * 1
    newexp1 ~ (nu2) * 1
    improve1 ~ (nu31) * 1
    learn2 ~ (nu1) * 1
    newexp2 ~ (nu2) * 1
    improve2 ~ (nu3) * 1
    learn3 ~ (nu1) * 1
    newexp3 ~ (nu2) * 1
    improve3 ~ (nu3) * 1
    # Unique factor variances
    learn1 ~~ (theta1) * learn1
    learn2 ~~ (theta1) * learn2
    learn3 ~~ (theta1) * learn3
    newexp1 ~~ (theta21) * newexp1
    newexp2 ~~ (theta2) * newexp2
    newexp3 ~~ (theta2) * newexp3
    improve1 ~~ (theta3) * improve1
    improve2 ~~ (theta3) * improve2
    improve3 ~~ (theta3) * improve3
    # Unique factor covariances
    learn1 ~~ learn2 + learn3
    learn2 ~~ learn3
    newexp1 ~~ newexp2 + newexp3
    newexp2 ~~ newexp3
    improve1 ~~ improve2 + improve3
    improve2 ~~ improve3
    # First factor variance and mean to 1 and 0
    eta1 ~~ 1 * eta1
    eta1 ~ 0
    # Other factor means to free
    eta2 ~ NA * 1
    eta3 ~ NA * 1
"
mis_fs_pgrowth <- get_fs(midus_merged, mis_pstrict_mod, method = "Bartlett")
mis_tspa_pgrowth <- get_tspa_mod(midus_merged, 
                                 vc = attr(mis_fs_pgrowth, "av_efs"), 
                                 fsA = attr(mis_fs_pgrowth, "fsA"))

mis_tspa_mod <- "
# Latent variables (indicated by factor scores)
eta1 =~ 1 * fs_eta1
eta2 =~ 1 * fs_eta2
eta3 =~ 1 * fs_eta3

# Constrain the errors
fs_eta1 ~~ 0.380595773255747 * fs_eta1
fs_eta2 ~~ 0.0345276093109072 * fs_eta1
fs_eta3 ~~ 0.00308371840280018 * fs_eta1
fs_eta2 ~~ 0.487058475519938 * fs_eta2
fs_eta3 ~~ 0.0536546180676673 * fs_eta2
fs_eta3 ~~ 0.47012592978719 * fs_eta3

# Constrain intercepts
fs_eta1 ~ 0 * 1
fs_eta2 ~ 0 * 1
fs_eta3 ~ 0 * 1

# Linear Growth Model
i =~ 1 * eta1 + 1 * eta2 + 1 * eta3
s =~ 0 * eta1 + 1 * eta2 + 2 * eta3

# Variance-covariances of intercepts and slopes
i ~~ phi1 * i
s ~~ s
i ~~ s

# Means of level and slope
i ~ 1
s ~ 1

# Fixed disturbances of latent outcomes to zero
eta1 ~ 0 * 1
eta2 ~ 0 * 1
eta3 ~ 0 * 1

# Constrain total variance at Wave 1
eta1 ~~ psi1 * eta1
psi1 + phi1 == 1
"
mis_tspa_fit <- sem(mis_tspa_mod, mis_fs_pgrowth)


# JSEM (misspecified) ----------------------------------------------------------

mis_jsem_mod <- "
# Factor loadings
eta1 =~ NA * learn1 + l1 * learn1 + l2 * newexp1 + l31 * improve1
eta2 =~ NA * learn2 + l1 * learn2 + l2 * newexp2 + l3 * improve2
eta3 =~ NA * learn3 + l1 * learn3 + l2 * newexp3 + l3 * improve3

# Measurement intercepts
learn1 ~ nu1 * 1
newexp1 ~ nu2 * 1
improve1 ~ nu31 * 1
learn2 ~ nu1 * 1
newexp2 ~ nu2 * 1
improve2 ~ nu3 * 1
learn3 ~ nu1 * 1
newexp3 ~ nu2 * 1
improve3 ~ nu3 * 1

# Unique factor variances
learn1 ~~ (theta1) * learn1
learn2 ~~ (theta1) * learn2
learn3 ~~ (theta1) * learn3
newexp1 ~~ (theta21) * newexp1
newexp2 ~~ (theta2) * newexp2
newexp3 ~~ (theta2) * newexp3
improve1 ~~ (theta3) * improve1
improve2 ~~ (theta3) * improve2
improve3 ~~ (theta3) * improve3

# Unique factor covariances
learn1 ~~ learn2 + learn3
learn2 ~~ learn3
newexp1 ~~ newexp2 + newexp3
newexp2 ~~ newexp3
improve1 ~~ improve2 + improve3
improve2 ~~ improve3

# Linear growth model
i =~ 1 * eta1 + 1 * eta2 + 1 * eta3
s =~ 0 * eta1 + 1 * eta2 + 2 * eta3

# Variance-covariances of intercepts and slopes
i ~~ phi1 * i
s ~~ s
i ~~ s

# Means of level (fixed to zero for identification) and slope
i ~ 0 * 1 
s ~ 1

# Disturbances of latent outcomes
eta1 ~ 0 * 1
eta2 ~ 0 * 1
eta3 ~ 0 * 1

# Constrain total variance at Wave 1
eta1 ~~ psi1 * eta1
psi1 + phi1 == 1
"
mis_jsem_fit <- sem(mis_jsem_mod, midus_merged)

# Comparison -------------------------------------------------------------------

msummary(
  list(
    "JSEM Growth" = jsem_fit, 
    "2S-PA Growth" = tspa_fit,
    "JSEM Growth (misspecified)" = mis_jsem_fit, 
    "2S-PA Grwoth (misspecified)" = mis_tspa_fit
  ),
  output = "markdown",
  statistic = "conf.int",
  coef_map = c(
    # "i ~1 " = "Mean(Level)",
    "s ~1 " = "Mean(Slope)",
    "i ~~ i" = "Var(Level)",
    "s ~~ s" = "Var(Slope)",
    "i ~~ s" = "Cov(Level, Slope)"
  )
)


# Additional second-order growth model -----------------------------------------

# Fix the first loading to obtain an interpretable metric

add_jsem_mod <- "
# Factor loadings
eta1 =~ 0.8477485 * learn1 + l2 * newexp1 + l31 * improve1
eta2 =~ 0.8477485 * learn2 + l2 * newexp2 + l3 * improve2
eta3 =~ 0.8477485 * learn3 + l2 * newexp3 + l3 * improve3

# Measurement intercepts
learn1 ~ nu1 * 1
newexp1 ~ nu2 * 1
improve1 ~ nu31 * 1
learn2 ~ nu1 * 1
newexp2 ~ nu2 * 1
improve2 ~ nu3 * 1
learn3 ~ nu1 * 1
newexp3 ~ nu2 * 1
improve3 ~ nu3 * 1

# Unique factor variances
learn1 ~~ (theta11) * learn1
learn2 ~~ (theta1) * learn2
learn3 ~~ (theta1) * learn3
newexp1 ~~ (theta21) * newexp1
newexp2 ~~ (theta2) * newexp2
newexp3 ~~ (theta2) * newexp3
improve1 ~~ (theta3) * improve1
improve2 ~~ (theta3) * improve2
improve3 ~~ (theta3) * improve3

# Unique factor covariances
learn1 ~~ learn2 + learn3
learn2 ~~ learn3
newexp1 ~~ newexp2 + newexp3
newexp2 ~~ newexp3
improve1 ~~ improve2 + improve3
improve2 ~~ improve3

# Linear growth model
i =~ 1 * eta1 + 1 * eta2 + 1 * eta3
s =~ 0 * eta1 + 1 * eta2 + 2 * eta3

# Variance-covariances of intercepts and slopes
i ~~ phi1 * i
s ~~ s
i ~~ s

# Means of level (fixed to zero for identification) and slope
i ~ 0 * 1 
s ~ 1

# Disturbances of latent outcomes
eta1 ~ 0 * 1
eta2 ~ 0 * 1
eta3 ~ 0 * 1

# Constrain total variance at Wave 1
eta1 ~~ psi1 * eta1
psi1 + phi1 == 1
"
add_jsem_fit <- sem(add_jsem_mod, midus_merged)
