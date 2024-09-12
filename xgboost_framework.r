install.packages("tidymodels")
install.packages("tidyverse")
install.packages("recipes")
install.packages("parsnip")
install.packages("rsample")
install.packages("yardstick")
install.packages("workflows")
install.packages("xgboost")
install.packages("finetune")
install.packages("lme4")
install.packages("purrr")
install.packages("stringr")

library(tidymodels)
library(tidyverse)
library(recipes)
library(parsnip)
library(rsample)
library(yardstick)
library(workflows)
library(xgboost)
library(finetune)
library(lme4)
library(purrr)
library(stringr)

# read data
data <- read.csv("CovidHospDataBrasil.csv", header = TRUE)

data <- data %>% 
  mutate(daysSinceStart = as.numeric(dateHosp - as.Date("2021-01-01")))

# clean data and get a vector of all Symptoms
data$sex <- as.factor(data$sex)

for (i in 1:length(data)) {
  if (grepl("date", colnames(data)[i])) {
    data[,i] = as.Date(data[,i])
  } else {
    if (i != 1 && i !=3 && i != 4) {
      data[,i] = as.logical(data[,i])
    }
  }
}
# get just the known information
known_data <- data[,-(21:29)]

# have known_data
# make TRUE and FALSE cols into 0's and 1's for xgb
for (i in 1:length(known_data)) {
  if (is.logical(known_data[,i])) {
    known_data[,i] = ifelse(known_data[,i] == TRUE, 1, 0)
  }
}

known_data$sex <- ifelse(known_data$sex == "M", 0, 1)
known_data$covidDeath <- ifelse(known_data$covidDeath == 1, "DIED", "SURVIVED")
known_data$covidDeath <- as.factor(known_data$covidDeath)

# remove patient_id and dateBirth
known_data <- known_data[, -(1:2)]


# partition into training, validation, and test
test_index <- sample(1:nrow(known_data), 0.2 * nrow(known_data))
test <- known_data[test_index,]
train <- known_data[-test_index,]

# create 10 folds for 10-fold cross validation
folds <- vfold_cv(known_data, strata = covidDeath)

folds

# make model recipe
rec <- recipe(
  covidDeath ~ age + sex + vaccine + fever + cough + daysSinceStart +
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic +
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity,
  data = known_data
) %>%
  step_dummy(all_nominal_predictors(), one_hot = TRUE)

prep(rec)
xgb_spec <-
  boost_tree(
    trees = tune(),
    min_n = tune(),
    mtry = tune(),
    tree_depth = tune(),
    learn_rate = tune(),
    loss_reduction = tune(),
    sample_size = tune(),
    stop_iter = tune()
  ) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

xgb_wf <- workflow(rec, xgb_spec)

split_category <- function(x) {
  x %>%
    stringr::str_split(", ") %>%
    purrr::map(stringr::str_remove_all, "[:punct:]") %>%
    purrr::map(stringr::str_squish) %>%
    purrr::map(stringr::str_to_lower) %>%
    purrr::map(stringr::str_replace_all, " ", "_")
}

doParallel::registerDoParallel()

xgb_rs <- tune_race_anova(
  xgb_wf,
  resamples = folds,
  grid = 200,
  metrics = metric_set(f_meas),
  control = control_race(verbose_elim = TRUE)
)

print(select_best(xgb_rs))

best_mtry <- 4
best_trees <- 343
best_treedepth <- 6
best_min_n <- 9    
best_learn_rate <- 1.15
best_loss_reduction <- 3108
best_sample_size <- 0.727
best_stop_iter <- 49


# make TRUE and FALSE cols into 0's and 1's for xgb
for (i in 1:length(data)) {
  if (is.logical(data[,i])) {
    data[,i] = ifelse(data[,i] == TRUE, 1, 0)
  }
}

data$sex <- ifelse(data$sex == "M", 0, 1)
data$covidDeath <- ifelse(data$covidDeath == 1, "DIED", "SURVIVED")
data$covidDeath <- as.factor(data$covidDeath)

xgb_spec <- boost_tree(
  trees = best_trees,
  min_n = best_min_n,
  mtry = best_mtry,
  tree_depth = best_treedepth,
  learn_rate = best_learn_rate,
  loss_reduction = best_loss_reduction,
  sample_size = best_sample_size,
  stop_iter = best_stop_iter
) %>%
  set_engine("xgboost") %>%
  set_mode("classification")
rec <- recipe(
  covidDeath ~ age + sex + vaccine + fever + cough + daysSinceStart +
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic +
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity,
  data = data
) %>%
  step_dummy(all_nominal_predictors(), one_hot = TRUE)

xgb_wf <- workflow() %>%
  add_recipe(rec) %>%
  add_model(xgb_spec)
final_fit <- fit(xgb_wf, data = data)

competition_data <- read.csv("CovidHospDataBrasil_test.csv", header = TRUE)

competition_data$sex <- as.factor(competition_data$sex)

for (i in 1:length(competition_data)) {
  if (grepl("date", colnames(competition_data)[i])) {
    competition_data[,i] = as.Date(competition_data[,i],
                                   format = "%d/%m/%Y")
  } else {
    if (i != 1 && i !=3 && i != 4) {
      competition_data[,i] = as.logical(competition_data[,i])
    }
  }
}

competition_data <- competition_data %>% 
  mutate(daysSinceStart = as.numeric(dateHosp - as.Date("2021-01-01")))

# make TRUE and FALSE cols into 0's and 1's for xgb
for (i in 1:length(competition_data)) {
  if (is.logical(competition_data[,i])) {
    competition_data[,i] = ifelse(competition_data[,i] == TRUE, 1, 0)
  }
}

competition_data$sex <- ifelse(competition_data$sex == "M", 0, 1)

predictions <- predict(final_fit, new_data = competition_data)
predictions <- ifelse(predictions == "DIED", TRUE, FALSE)

submission <- data.frame(
  Patient_id = competition_data$Patient_id,
  covidDeath = predictions
)
write.csv(submission, file = "submission.csv", row.names = FALSE)
