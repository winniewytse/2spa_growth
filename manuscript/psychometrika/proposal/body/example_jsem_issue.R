# Proposal example (issue with JSEM)

# Setup ------------------------------------------------------------------------

library(tidyverse)
library(lavaan)
library(R2spa)

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


# Personal growth model (Wave 1) -----------------------------------------------

pg_mod <- "
pg1 =~ learn1 + newexp1 + improve1
"
pg_fit <- cfa(pg_mod, midus_merged, std.lv = TRUE)
parameterestimates(pg_fit)

# Purpose of life model (Wave 1) -----------------------------------------------

pl_mod <- "
pl2 =~ purpose2 + future2 + done2
"
pl_fit <- cfa(pl_mod, midus_merged, std.lv = TRUE)
parameterestimates(pl_fit)

# Joint model ------------------------------------------------------------------

joint_mod <- "
pg1 =~ learn1 + newexp1 + improve1
pl2 =~ purpose2 + future2 + done2

pl2 ~ pg1
"
joint_fit <- sem(joint_mod, midus_merged, std.lv = TRUE)
parameterestimates(joint_fit)
