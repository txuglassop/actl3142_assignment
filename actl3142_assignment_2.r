library(tidyverse)
library(mgcv)
library(boot)
library(caret)
library(gam)
library(car)
library(broom)
library(glmnet)
library(psych)
library(ggplot2)
library(corrplot)
library(ranger)
library(rpart)
library(rpart.plot)
library(tree)
library(pROC)
library(gridExtra)
library(glmnetUtils)


# read data
data <- read.csv("CovidHospDataBrasil.csv", header = TRUE)

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

data <- transform(
  data, 
  "daysInHosp" = as.numeric(dateEndObs - dateHosp)
)

non_survivors <- data %>%
  filter(covidDeath) 

survivors <- data %>%
  filter(!covidDeath)

survivors$group <- "survivors"
non_survivors$group <- "nonsurvivors"

combined_data <- rbind(survivors, non_survivors)

ggplot(combined_data, aes(x = group, y = daysInHosp, fill = group)) +
  geom_boxplot(alpha = 0.5) +
  labs(
    title = "Days Spent in Hospital of Survivors and Non-Survivors",
    fill = "Group",
    x = "Group",
    y = "Days Spent in Hospital"
  ) +
  theme_minimal() +
  scale_fill_manual(
    values = c("red", "blue"),
    labels = c("Non-Survivors", "Survivors")
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

filtered_data <- data %>%
  filter(daysInHosp < 40)


non_survivors <- filtered_data %>%
  filter(covidDeath) 

survivors <- filtered_data %>%
  filter(!covidDeath)

survivors$group <- "survivors"
non_survivors$group <- "nonsurvivors"

combined_data <- rbind(survivors, non_survivors)

ggplot(combined_data, aes(x = group, y = daysInHosp, fill = group)) +
  geom_boxplot(alpha = 0.5) +
  labs(
    title = "Days Spent in Hospital of Survivors and Non-Survivors",
    fill = "Group",
    x = "Group",
    y = "Days Spent in Hospital"
  ) +
  theme_minimal() +
  scale_fill_manual(
    values = c("red", "blue"),
    labels = c("Non-Survivors", "Survivors")
  )

data <- filtered_data


################################## Part 1 ######################################

# in this section, we will compare different models using all available data
# with the goal of inference in mind. Models will be compared with direct
# methods and machine learning methods will be applied to models


### Partition data into train and test
test_index <- sample(1:nrow(data), 0.2 * nrow(data))
test <- data[test_index,]
train <- data[-test_index,]


####### The search for interaction terms

# visualise correlations between predictors
medical_predictors <- data %>% 
  dplyr::select(-covidDeath, -Patient_id, -dateBirth, -age, -sex, -dateHosp,
                -dateAdmIcu, -dateDisIcu, -dateEndObs, -daysInHosp)

# Calculate the Phi coefficient matrix
cor_matrix <- psych::polychoric(medical_predictors)

cor_matrix$rho

corrplot(cor_matrix$rho, method = "circle", tl.cex = 0.7, number.cex = 0.7, 
         tl.col = "black", tl.srt = 45, tl.names = colnames(medical_predictors))
names(medical_predictors)

# most significant correlation between
# dyspnoea and respdistress
# oxygensat and respdistress
# dyspnoea and oxygensat
# diarrhea and vomit

# control model with no interactions
control_model <- glm(
  covidDeath ~ dateBirth + sex + vaccine + dateHosp + fever + cough + 
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
    icu + respdistress + cardio + hepatic + immuno + renal + daysInHosp,
  data = train,
  family = binomial
)

# interaction between dyspnoea * oxygensat
inter_model1 <- glm(
  covidDeath ~ dateBirth + sex + vaccine + dateHosp + fever + cough + 
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
    icu + respdistress + cardio + hepatic + immuno + renal + daysInHosp +
    dyspnoea * oxygensat,
  data = train,
  family = binomial
)


# interaction between dyspnoea * respdistress, oxygensat * respdistress,
# dyspnoea * oxygensat
inter_model2 <- glm(
  covidDeath ~ dateBirth + sex + vaccine + dateHosp + fever + cough + 
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
    icu + respdistress + cardio + hepatic + immuno + renal + daysInHosp +
    dyspnoea * respdistress + oxygensat * respdistress + dyspnoea * oxygensat,
  data = train,
  family = binomial
)

# interaction between oxygensat * respdistress, dyspnoea * respdistress
inter_model3 <- glm(
  covidDeath ~ dateBirth + sex + vaccine + dateHosp + fever + cough + 
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
    icu + respdistress + cardio + hepatic + immuno + renal + daysInHosp +
    oxygensat * respdistress + dyspnoea * respdistress,
  data = train,
  family = binomial
)

# interaction between oxygensat * respdistress
inter_model4 <- glm(
  covidDeath ~ dateBirth + sex + vaccine + dateHosp + fever + cough + 
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
    icu + respdistress + cardio + hepatic + immuno + renal + daysInHosp +
    oxygensat * respdistress,
  data = train,
  family = binomial
)

# interaction between dyspnoea * respdistress
inter_model5 <- glm(
  covidDeath ~ dateBirth + sex + vaccine + dateHosp + fever + cough + 
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
    icu + respdistress + cardio + hepatic + immuno + renal + daysInHosp +
    dyspnoea * respdistress,
  data = train,
  family = binomial
)


# cost function using 0.5 as threshold
error_cost <- function(r, pi){
  weight1 = 1 
  weight0 = 1 
  c1 = (r==1)&(pi<0.5) 
  c0 = (r==0)&(pi>0.5) 
  return(mean(weight1*c1+weight0*c0))
}

# find 10 fold cv error for all these models
cv_control <- cv.glm(train, control_model, cost = mycost, K = 10)
cv_1 <- cv.glm(train, inter_model1, cost = error_cost, K = 10)
cv_2 <- cv.glm(train, inter_model2, cost = error_cost, K = 10)
cv_3 <- cv.glm(train, inter_model3, cost = error_cost, K = 10)
cv_4 <- cv.glm(train, inter_model4, cost = error_cost, K = 10)
cv_5 <- cv.glm(train, inter_model5, cost = error_cost, K = 10)

cv_errors <- data.frame(
  model = 0:5,
  cv = c(cv_control$delta[1], cv_1$delta[1], cv_2$delta[1], 
         cv_3$delta[1], cv_4$delta[1], cv_5$delta[1])
)

cv_errors
# model 3 is the best

# another interaction that may be good is diarrhea and vomit
# reshuffle data and repeat process
test_index <- sample(1:nrow(data), 0.2 * nrow(data))
test <- data[test_index,]
train <- data[-test_index,]

# control model with no interactions
control_model <- glm(
  covidDeath ~ dateBirth + sex + vaccine + dateHosp + fever + cough + 
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
    icu + respdistress + cardio + hepatic + immuno + renal + daysInHosp,
  data = train,
  family = binomial
)

# interaction between oxygensat * respdistress, dyspnoea * respdistress
inter_model1 <- glm(
  covidDeath ~ dateBirth + sex + vaccine + dateHosp + fever + cough + 
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
    icu + respdistress + cardio + hepatic + immuno + renal + daysInHosp +
    oxygensat * respdistress + dyspnoea * respdistress,
  data = train,
  family = binomial
)

# interaction between oxygensat * respdistress, dyspnoea * respdistress
# diarrhea * vomit
inter_model2 <- glm(
  covidDeath ~ dateBirth + sex + vaccine + dateHosp + fever + cough + 
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
    icu + respdistress + cardio + hepatic + immuno + renal + daysInHosp +
    oxygensat * respdistress + dyspnoea * respdistress + diarrhea * vomit,
  data = train,
  family = binomial
)

# interaction between diarrhea * vomit
inter_model3 <- glm(
  covidDeath ~ dateBirth + sex + vaccine + dateHosp + fever + cough + 
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
    icu + respdistress + cardio + hepatic + immuno + renal + daysInHosp +
    diarrhea * vomit,
  data = train,
  family = binomial
)

cv_control <- cv.glm(train, control_model, cost = error_cost, K = 10)
cv_1 <- cv.glm(train, inter_model1, cost = error_cost, K = 10)
cv_2 <- cv.glm(train, inter_model2, cost = error_cost, K = 10)
cv_3 <- cv.glm(train, inter_model3, cost = error_cost, K = 10)

cv_errors <- data.frame(
  model = 0:3,
  cv = c(cv_control$delta[1], cv_1$delta[1], cv_2$delta[1], cv_3$delta[1])
)

cv_errors

# again, model with interactions outperforms control
# improvement in adding diarrhea * vomit


####### Model 1
####### glm with all predictors and interactions and lasso regression with lambda.1se

data_1 <- model.matrix(
  covidDeath ~ dateBirth + sex + vaccine + dateHosp + fever + cough + 
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
    icu + respdistress + cardio + hepatic + immuno + renal + daysInHosp +
    oxygensat * respdistress + dyspnoea * respdistress + diarrhea * vomit,
  data = train,
  family = binomial()
)

x <- data.matrix(data_1)

y <- train$covidDeath

lambda_1 <- cv.glmnet(x, y, alpha = 1, family = binomial())

plot(lambda_1)

model1 <- glmnet(
  x, y, alpha = 1, lambda = lambda_1$lambda.1se, family = binomial
)

coef(model1, s = lambda_1$lambda.1se)



####### Model 2
####### gam with some polynomial regression on numerical predictors
####### and analysis of deviance

# first, find which polynomials work best for each numerical predictor

# age

cv_error_vec_age <- c()

for (dof in 1:8) {
  model_temp <- glm(covidDeath ~ ns(age, df = dof), data = train, family = binomial())
  cv_error_vec_age <- c(cv_error_vec_age, cv.glm(train, model_temp, cost = error_cost, K = 10)$delta[1])
}

cv_error_age <- data.frame(
  DoF = 1:8,
  CV_Error = cv_error_vec_age
)

age_plot <- ggplot(cv_error_age, aes(x = DoF, y = CV_Error)) +
  geom_line() +
  geom_point() +
  labs(title = "10-fold Cross-Validation for age",
       x = "Degrees of Freedom",
       y = "Cross-Validation Error") +
  theme_minimal() +
  scale_x_continuous(breaks = 1:8)


# dateHosp

cv_error_vec_date <- c()

for (dof in 1:8) {
  model_temp <- glm(covidDeath ~ ns(dateHosp, df = dof), data = train, family = binomial())
  cv_error_vec_date <- c(cv_error_vec_date, cv.glm(train, model_temp, cost = error_cost, K = 10)$delta[1])
}

cv_error_date <- data.frame(
  DoF = 1:8,
  CV_Error = cv_error_vec_date
)

date_plot <- ggplot(cv_error_date, aes(x = DoF, y = CV_Error)) +
  geom_line() +
  geom_point() +
  labs(title = "10-fold Cross-Validation for dateHosp",
       x = "Degrees of Freedom",
       y = "Cross-Validation Error") +
  theme_minimal() +
  scale_x_continuous(breaks = 1:8)


# daysInHosp

cv_error_vec_days <- c()

for (dof in 1:6) {
  model_temp <- glm(covidDeath ~ ns(daysInHosp, df = dof), data = train, family = binomial())
  cv_error_vec_days <- c(cv_error_vec_days, cv.glm(train, model_temp, cost = error_cost, K = 10)$delta[1])
}

cv_error_days <- data.frame(
  DoF = 1:6,
  CV_Error = cv_error_vec_days
)

days_plot <- ggplot(cv_error_days, aes(x = DoF, y = CV_Error)) +
  geom_line() +
  geom_point() +
  labs(title = "10-fold Cross-Validation for daysInHosp",
       x = "Degrees of Freedom",
       y = "Cross-Validation Error") +
  theme_minimal() +
  scale_x_continuous(breaks = 1:8)

grid.arrange(age_plot, date_plot, days_plot, nrow = 1, ncol = 3)


# from the findings
# age - dof = 2
# dateHosp - dof = 3
# daysInHosp - degree = 4

model2 <- glm(
  covidDeath ~ ns(age, df = 2) + sex + vaccine + ns(dateHosp, df = 3) + fever + cough + 
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
    icu + respdistress + cardio + hepatic + immuno + renal + ns(daysInHosp, df = 4) +
    oxygensat * respdistress + dyspnoea * respdistress + diarrhea * vomit,
  data = train,
  family = binomial()
)

summary(model2)

cv.glm(train, model2, cost = error_cost, K = 10)$delta[1]

print(anova(model2, test = "Chisq"))

# remove cardio

model2 <- glm(
  covidDeath ~ ns(age, df = 2) + sex + vaccine + ns(dateHosp, df = 3) + fever + cough + 
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
    icu + respdistress + hepatic + immuno + renal + ns(daysInHosp, df = 4) +
    oxygensat * respdistress + dyspnoea * respdistress + diarrhea * vomit,
  data = train,
  family = binomial()
)

print(anova(model2, test = "Chisq"))

# remove respdistress * dyspnoea

model2 <- glm(
  covidDeath ~ ns(age, df = 2) + sex + vaccine + ns(dateHosp, df = 3) + fever + cough + 
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
    icu + respdistress + hepatic + immuno + renal + ns(daysInHosp, df = 4) +
    oxygensat * respdistress + diarrhea * vomit,
  data = train,
  family = binomial()
)

print(anova(model2, test = "Chisq"))

summary(model2)

cv.glm(train, model2, cost = error_cost, K = 10)$delta[1]



##### Model 3
##### tree

model3 <- rpart(
  covidDeath ~ age + sex + vaccine + dateHosp + fever + cough + 
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
    icu + respdistress + cardio + hepatic + immuno + renal + daysInHosp,
  data = data,
  method = "class",
  control = rpart.control(cp = 0)
)

## gets a very crazy looking tree....

model3 <- rpart(
  covidDeath ~ age + sex + vaccine + dateHosp + fever + cough + 
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
    icu + respdistress + cardio + hepatic + immuno + renal + daysInHosp,
  data = data,
  method = "class",
  control = rpart.control(cp = 0.0001, minsplit = 5, maxdepth = 30)
)


summary(model3)

rpart.plot(model4, main = "Regression Tree for COVID Deaths", extra = 106)

printcp(model3)
plotcp(model3)


model3 <- rpart(
  covidDeath ~ age + sex + vaccine + dateHosp + fever + cough + 
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
    icu + respdistress + cardio + hepatic + immuno + renal + daysInHosp,
  data = data,
  method = "class",
  control = rpart.control(cp = 0.0001)
)

rpart.plot(model3, main = "Regression Tree for COVID Deaths")

printcp(model3)
plotcp(model3)

model3_pruned <- prune(model3, cp = 0.00078841)

rpart.plot(model3_pruned, main = "Trimmed Tree")


model3 <- rpart(
  covidDeath ~ age + sex + vaccine + dateHosp + fever + cough + 
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
    icu + respdistress + cardio + hepatic + immuno + renal + daysInHosp,
  data = data,
  method = "class",
  control = rpart.control(cp = 0.00078841)
)

rpart.plot(model3, main = "Regression Tree for COVID Deaths")


###### reshuffle data and do cv to find the best model out of the 3

# define a common set of 10 folds for cv

cv_folds <- createFolds(data$covidDeath, k = 10, list = TRUE, returnTrain = TRUE)

# functions to perform a validation set on all the models

cv_model1 <- function(train, test, threshold) {
  x <- model.matrix(
    covidDeath ~ dateBirth + sex + vaccine + dateHosp + fever + cough + 
      sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
      downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
      icu + respdistress + cardio + hepatic + immuno + renal + daysInHosp +
      oxygensat * respdistress + dyspnoea * respdistress + diarrhea * vomit,
    train,
    family = binomial()
  )
  y <- train$covidDeath
  
  model <- cv.glmnet(x, y, alpha = 1, family = binomial())
  
  x_test <- model.matrix(
    covidDeath ~ dateBirth + sex + vaccine + dateHosp + fever + cough + 
      sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
      downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
      icu + respdistress + cardio + hepatic + immuno + renal + daysInHosp +
      oxygensat * respdistress + dyspnoea * respdistress + diarrhea * vomit,
    test,
    family = binomial()
  )
  
  predictions <- predict(model, s = "lambda.1se", newx = x_test, type = "response")
  
  return(ifelse(predictions > threshold, TRUE, FALSE))
}


cv_model2 <- function(train, test,threshold) {
  model <- glm(
    covidDeath ~ ns(age, df = 2) + sex + vaccine + ns(dateHosp, df = 3) + fever + cough + 
      sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
      downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
      icu + respdistress + hepatic + immuno + renal + ns(daysInHosp, df = 4) +
      oxygensat * respdistress + diarrhea * vomit,
    data = train,
    family = binomial()
  )
  
  predictions <- predict(model, test, type = "response")
  
  return(ifelse(predictions > threshold, TRUE, FALSE))
}


cv_model3 <- function(train, test) {
  model <- rpart(
    covidDeath ~ age + sex + vaccine + dateHosp + fever + cough + 
      sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
      downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
      icu + respdistress + cardio + hepatic + immuno + renal + daysInHosp,
    data = train,
    method = "class",
    control = rpart.control(cp = 0.00078841)
  )
  
  predictions <- predict(model, test, type = "class")
  
  return(predictions)
}

error_rates_model1 <- numeric(length(cv_folds))
error_rates_model2 <- numeric(length(cv_folds))
error_rates_model3 <- numeric(length(cv_folds))


for (i in 1:length(cv_folds)) {
  train_indices <- cv_folds[[i]]
  test_indices <- setdiff(seq_len(nrow(data)), train_indices)
  
  train_data <- data[train_indices,]
  test_data <- data[test_indices,]
  
  # model 1
  predictions_model1 <- cv_model1(train_data, test_data, 0.5)
  error_rates_model1[i] <- mean(predictions_model1 != test_data$covidDeath)
  
  # model 2
  predictions_model2 <- cv_model2(train_data, test_data, 0.5)
  error_rates_model2[i] <- mean(predictions_model2 != test_data$covidDeath)
  
  # model 3
  predictions_model3 <- cv_model3(train_data, test_data)
  error_rates_model3[i] <- mean(predictions_model3 != test_data$covidDeath)
}

# aggregate findings
cv_error <- data.frame(
  Model = c(1,2,3),
  CV = c(
    mean(error_rates_model1), mean(error_rates_model2), mean(error_rates_model3)
  )
)

cv_error

# repeat this but with f1 score

f1_model1 <- numeric(length(cv_folds))
f1_model2 <- numeric(length(cv_folds))
f1_model3 <- numeric(length(cv_folds))

for (i in 1:length(cv_folds)) {
  train_indices <- cv_folds[[i]]
  test_indices <- setdiff(seq_len(nrow(data)), train_indices)
  
  train_data <- data[train_indices,]
  test_data <- data[test_indices,]
  
  # model 1
  predictions_model1 <- cv_model1(train_data, test_data, 0.3)
  
  confusion_matrix <- table(Predicted = predictions_model1, Actual = test_data$covidDeath)
  
  precision <- confusion_matrix[2,2] / (confusion_matrix[2,2] + confusion_matrix[2,1])
  recall <- confusion_matrix[2,2] / (confusion_matrix[2,2] + confusion_matrix[1,2])
  
  f1_model1[i] <- 2 * precision * recall / (precision + recall)
  
  # model 2
  predictions_model2 <- cv_model2(train_data, test_data, 0.3)
  
  confusion_matrix <- table(Predicted = predictions_model2, Actual = test_data$covidDeath)
  
  precision <- confusion_matrix[2,2] / (confusion_matrix[2,2] + confusion_matrix[2,1])
  recall <- confusion_matrix[2,2] / (confusion_matrix[2,2] + confusion_matrix[1,2])
  
  f1_model2[i] <- 2 * precision * recall / (precision + recall)
  
  # model 3
  predictions_model3 <- cv_model3(train_data, test_data)
  
  confusion_matrix <- table(Predicted = predictions_model3, Actual = test_data$covidDeath)
  
  precision <- confusion_matrix[2,2] / (confusion_matrix[2,2] + confusion_matrix[2,1])
  recall <- confusion_matrix[2,2] / (confusion_matrix[2,2] + confusion_matrix[1,2])
  
  f1_model3[i] <- 2 * precision * recall / (precision + recall)
}

# aggregate findings
cv_f1 <- data.frame(
  Model = c(1,2,3),
  CV = c(
    mean(f1_model1), mean(f1_model2), mean(f1_model3)
  )
)

cv_f1




# plot a roc curve and find auc

x_test <- model.matrix(
  covidDeath ~ dateBirth + sex + vaccine + dateHosp + fever + cough + 
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
    icu + respdistress + cardio + hepatic + immuno + renal + daysInHosp +
    oxygensat * respdistress + dyspnoea * respdistress + diarrhea * vomit,
  test,
  family = binomial()
)

model1_pred <- predict(model1, s = "lambda.1se", newx = x_test, type = "response")

model2_pred <- predict(model2, test, type = "response")

model3_pred <- predict(model3, test, type = "prob")[,2] #get Pr(TRUE)

roc_model1 <- roc(test$covidDeath, model1_pred[, "s1"])
auc_model1 <- auc(roc_model1)

# calculate ROC and AUC for GAM with LASSO
roc_model2 <- roc(test$covidDeath, model2_pred)
auc_model2 <- auc(roc_model2)

# calculate ROC and AUC for Decision Tree
roc_model3 <- roc(train$covidDeath, model3_pred)
auc_model3 <- auc(roc_model3)

plot(roc_model1, col = "blue", main = "ROC Curves")
lines(roc_model2, col = "red")
legend("bottomright", 
       legend = c(
         paste("Model 1, AUC = ", round(auc_model1,3)), 
         paste("Model 2, AUC = ", round(auc_model2,3))),
       col = c("blue", "red"), lwd = 2)

# use the test set to measure performance
x_test <- model.matrix(
  covidDeath ~ dateBirth + sex + vaccine + dateHosp + fever + cough + 
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
    icu + respdistress + cardio + hepatic + immuno + renal + daysInHosp +
    oxygensat * respdistress + dyspnoea * respdistress + diarrhea * vomit,
  test,
  family = binomial()
)

model1_pred <- predict(model1, x_test, type = "response")
model1_pred <- ifelse(model1_pred > 0.5, TRUE, FALSE)

model2_pred <- predict(model2, test, type = "response")
model2_pred <- ifelse(model2_pred > 0.3, TRUE, FALSE)

model3_pred <- predict(model3, test, type = "class")

model1_test_err <- mean(model1_pred != test$covidDeath)
model2_test_err <- mean(model2_pred != test$covidDeath)
model3_test_err <- mean(model3_pred != test$covidDeath)

test_errors <- data.frame(
  Model = c(1,2,3),
  Test_Err = c(
    model1_test_err, model2_test_err, model3_test_err
  )
)

test_errors


# assess model fit

filtered_data <- subset(data, daysInHosp <= 20)

best_model <- glm(
  covidDeath ~ ns(age, df = 2) + sex + vaccine + ns(dateHosp, df = 3) + fever + cough + 
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
    icu + respdistress + hepatic + immuno + renal + ns(daysInHosp, df = 4) +
    oxygensat * respdistress + diarrhea * vomit,
  data = filtered_data,
  family = binomial()
)

summary(best_model)

# assess appropriateness of fit

residuals <- resid(best_model)
fitted_values <- fitted(best_model)

residual_data <- data.frame(
  age = filtered_data$age, 
  dateHosp = filtered_data$dateHosp,
  daysInHosp = filtered_data$daysInHosp,
  residuals = residuals
)

age_residuals <- ggplot(residual_data, aes(x = age, y = residuals)) +
  geom_point() +
  geom_smooth(method = "gam", col = "blue") + 
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  labs(title = "Residuals vs Age", x = "Age", y = "Residuals") +
  theme_minimal()

date_residuals <- ggplot(residual_data, aes(x = dateHosp, y = residuals)) +
  geom_point() +
  geom_smooth(method = "gam", col = "blue") + 
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  labs(title = "Residuals vs Date Hospitalised", x = "Date Hospitalised", y = "Residuals") +
  theme_minimal()

days_residuals <- ggplot(residual_data, aes(x = daysInHosp, y = residuals)) +
  geom_point() +
  geom_smooth(method = "gam", col = "blue") + 
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") + 
  labs(title = "Residuals vs Days Spent in Hospital", x = "Days Spent in Hospital", y = "Residuals") +
  theme_minimal()

grid.arrange(age_residuals, date_residuals, days_residuals, nrow = 1, ncol = 3)


# check VIF for multicollinearity

car::vif(best_model, type = "predictor")

###### Inferences

# plot each splines

# for dateHosp

dateHosp_seq <- seq(min(data$dateHosp), max(data$dateHosp), length.out = 100)

spline_basis_dateHosp <- ns(dateHosp_seq, df = 3)

dateHosp_coef <- c(-0.35147, 0.69158, -0.62602)

log_odds_dateHosp <- spline_basis_dateHosp %*% dateHosp_coef

plot_data_dateHosp <- data.frame(dateHosp = dateHosp_seq, log_odds = log_odds_dateHosp)

dateHosp_spline <- ggplot(plot_data_dateHosp, aes(x = dateHosp, y = log_odds)) +
  geom_line(color = "red", size = 1) +
  labs(title = "Spline for Date of Hospitalization",
       x = "Date of Hospitalization",
       y = "Log Odds Contribution") +
  theme_minimal()

# show hospital admissions and deaths over time

non_survivors <- data %>%
  filter(covidDeath) 

survivors <- data %>%
  filter(!covidDeath)

admissions_by_date <- data %>%
  group_by(dateHosp) %>%
  summarize(Admissions = n()) %>%
  rename(Date = dateHosp)

deaths_by_date <- non_survivors %>%
  group_by(dateEndObs) %>%
  summarize(Deaths = n()) %>%
  rename(Date = dateEndObs)

combined_data <- full_join(admissions_by_date, deaths_by_date, by = "Date") %>%
  replace_na(list(Admissions = 0, Deaths = 0))

combined_data_long <- combined_data %>%
  pivot_longer(cols = c("Admissions", "Deaths"), names_to = "Metric", values_to = "Count")

admissions <- ggplot(combined_data_long, aes(x = Date, y = Count, color = Metric)) +
  geom_line(size = 1) +
  labs(
    title = "Hospital Admissions and Deaths Over Time",
    x = "Date",
    y = "Count",
    color = "Metric"
  ) +
  scale_color_manual(values = c("Deaths" = "red", "Admissions" = "black"))
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


grid.arrange(dateHosp_spline, admissions, nrow = 1, ncol = 2)

# for daysInHosp

quantile(data$daysInHosp, probs = 0.95)
# 35
daysInHosp_seq <- seq(min(data$daysInHosp), 35, length.out = 100)

spline_basis_daysInHosp <- ns(daysInHosp_seq, df = 4)

daysInHosp_coef <- c(-1.864990, -0.795268, -3.311231, -0.040740)

log_odds_daysInHosp <- spline_basis_daysInHosp %*% daysInHosp_coef

plot_data_daysInHosp <- data.frame(daysInHosp = daysInHosp_seq, log_odds = log_odds_daysInHosp)

daysInHosp_spline <- ggplot(plot_data_daysInHosp, aes(x = daysInHosp, y = log_odds)) +
  geom_line(color = "green", size = 1) +
  labs(title = "Days in Hospital Spline",
       x = "Days in Hospital",
       y = "Log Odds Contribution") +
  theme_minimal()

# for age

age_seq <- seq(min(data$age), max(data$age), length.out = 100)

spline_basis_age <- ns(age_seq, df = 4)

age_coef <- c(1.62077, 2.45040, 4.31564, 4.08318)

log_odds_age <- spline_basis_age %*% age_coef

plot_data_age <- data.frame(age = age_seq, log_odds = log_odds_age)

age_spline <- ggplot(plot_data_age, aes(x = age, y = log_odds)) +
  geom_line(color = "blue", size = 1) +
  labs(title = "Age Spline",
       x = "Age",
       y = "Log Odds Contribution") +
  theme_minimal()


grid.arrange(daysInHosp_spline, age_spline, nrow = 1, ncol = 2)

mean(survivors$daysInHosp)
mean(non_survivors$daysInHosp)

median(survivors$daysInHosp)
median(non_survivors$daysInHosp)

# retrain model on scaled data and order coefficients

data_scaled <- filtered_data %>%
  mutate_if(is.numeric, scale)

scaled_model <- glm(
  covidDeath ~ ns(age, df = 2) + sex + vaccine + ns(dateHosp, df = 3) + fever + cough + 
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity +
    icu + respdistress + hepatic + immuno + renal + ns(daysInHosp, df = 4) +
    oxygensat * respdistress + diarrhea * vomit,
  data = data_scaled,
  family = binomial()
)

sort(scaled_model$coefficients)


############################### PART 2 ########################################

# get just the known information
known_data <- data[,-(21:29)]

names(known_data)

f1_cost <- function(r, pi) {
  predictions <- ifelse(pi > 0.3, 1, 0)
  
  tp <- sum((predictions == 1) & (r == 1))
  fp <- sum((predictions == 1) & (r == 0))
  fn <- sum((predictions == 0) & (r == 1))
  
  precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  
  f1 <- ifelse((precision + recall) == 0, 0, 2 * (precision * recall) / (precision + recall))
  
  print(f1)
  return(f1)
}


# model 4
# build upon model 2 by adding elastic net

data_4 <- model.matrix(
  covidDeath ~ ns(age, df = 2) + sex + vaccine + ns(dateHosp, df = 3) + fever + cough + 
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity,
  data = known_data,
  family = binomial()
)

x <- data.matrix(data_4)

y <- known_data$covidDeath

model4_tuning <- cva.glmnet(x, y, cost = f1_cost, family = binomial())


plot(model4_tuning)

model4_tuning

# make predictions on test set

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

training.x <- model.matrix(
  covidDeath ~ ns(age, df = 2) + sex + vaccine + ns(dateHosp, df = 3) + fever + cough + 
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity,
  data = known_data
)

final.x <- model.matrix(
  ~ ns(age, df = 2) + sex + vaccine + ns(dateHosp, df = 3) + fever + cough + 
    sorethroat + dyspnoea + oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + diabetes + neurological + pneumopathy + obesity,
  data = competition_data
)

x <- data.matrix(training.x)

y <- known_data$covidDeath

model4 <- glmnet(x, y, alpha = 0.008, lambda = 0.004367, family = binomial())

model4_predictions <- predict(model4, data.matrix(final.x), type = "response")
model4_predictions <- ifelse(model4_predictions > 0.3, TRUE, FALSE)

submission <- data.frame(
  Patient_id = competition_data$Patient_id,
  covidDeath = model4_predictions
)
write.csv(submission, file = "submission_model4.csv", row.names = FALSE)


# model 5
# local regression on dateHosp

# find optimal span

f1_vec <- c()
for (span in seq(0.05, 0.5, 0.05)) {
  temp_model <- gam(
    covidDeath ~ ns(age, df = 2) + sex + vaccine + lo(dateHosp, span = span) +
      fever + cough + sorethroat + dyspnoea + oxygensat + diarrhea + vomit +
      hematologic + downsyn + asthma + diabetes + neurological + pneumopathy +
      obesity,
    data = known_data,
    family = binomial
  )
  
  f1_vec <- c(f1_vec, cv.glm(known_data, temp_model, cost = f1_cost, K = 10)$delta[1])
}

span_f1 <- data.frame(
  span = seq(0.05, 0.5, 0.05),
  f1 = f1_vec
)

span_f1

span_plot <- ggplot(span_f1, aes(x = span, y = f1)) +
  geom_point() +         
  geom_line(colour = "blue") +          
  labs(
    title = "CV F1 Score vs Span for LOESS on dateHosp",
    x = "Span",
    y = "F1 Score"
  ) +
  theme_minimal()

# find optimal sp

f1_vec <- c()
for (sp in seq(0, 2, 0.1)) {
  temp_model <- gam(
    covidDeath ~ s(age, sp = sp) + sex + vaccine + ns(dateHosp, df = 3) +
      fever + cough + sorethroat + dyspnoea + oxygensat + diarrhea + vomit +
      hematologic + downsyn + asthma + diabetes + neurological + pneumopathy +
      obesity,
    data = known_data,
    family = binomial
  )
  
  f1_vec <- c(f1_vec, cv.glm(known_data, temp_model, cost = f1_cost, K = 10)$delta[1])
}

sp_f1 <- data.frame(
  sp = seq(0, 2, 0.1),
  f1 = f1_vec
)

sp_plot <- ggplot(span_f1, aes(x = sp, y = f1)) +
  geom_point() +         
  geom_line(colour = "blue") +          
  labs(
    title = "CV F1 Score vs Smoothing Parameter for age",
    x = "Smoothing Parameter",
    y = "F1 Score"
  ) +
  theme_minimal()

grid.arrange(span_plot, sp_plot, ncol = 2, nrow = 1)


# from the results = span = 0.15, sp = 0.2
# apply elastic net

data_5 <- model.matrix(
  covidDeath ~ s(age, sp = 0.2) + sex + vaccine + lo(dateHosp, df = 0.15) + fever + cough + 
    sorethroat + dyspnoea * oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + neurological + pneumopathy + obesity * diabetes,
  data = known_data,
  family = binomial()
)

x <- data.matrix(data_5)

y <- known_data$covidDeath

model5_tuning <- cva.glmnet(x, y, cost = f1_cost, family = binomial())


plot(model5_tuning)

model5_tuning

# make predictions on the test set

training.x <- model.matrix(
  covidDeath ~ s(age, sp = 0.2) + sex + vaccine + lo(dateHosp, df = 0.15) + fever + cough + 
    sorethroat + dyspnoea * oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + neurological + pneumopathy + obesity * diabetes,
  data = known_data
)

final.x <- model.matrix(
  ~ s(age, sp = 0.2) + sex + vaccine + lo(dateHosp, df = 0.15) + fever + cough + 
    sorethroat + dyspnoea * oxygensat + diarrhea + vomit + hematologic + 
    downsyn + asthma + neurological + pneumopathy + obesity * diabetes,
  data = competition_data
)

x <- data.matrix(training.x)

y <- known_data$covidDeath

model5 <- glmnet(x, y, alpha = 0.001, lambda = 0.01161, family = binomial())

model5_predictions <- predict(model5, data.matrix(final.x), type = "response")
model5_predictions <- ifelse(model5_predictions > 0.3, TRUE, FALSE)

submission <- data.frame(
  Patient_id = competition_data$Patient_id,
  covidDeath = model5_predictions
)
write.csv(submission, file = "submission_model5.csv", row.names = FALSE)


# model 6
# xgboost trained with race

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

xgb_grid <- grid_latin_hypercube(
  trees(range = c(300, 500)),        
  min_n(range = c(1, 20)),         
  mtry(range = c(1, 10)),             
  tree_depth(range = c(3, 10)),       
  learn_rate(range = c(0.01, 0.3)),      
  loss_reduction(range = c(0, 10)),      
  sample_size = sample_prop(range = c(0.5, 1)),
  stop_iter(range = c(10, 50)),       
  size = 150                              
)

xgb_rs <- tune_race_anova(
  xgb_wf,
  resamples = folds,
  grid = xgb_grid,
  metrics = metric_set(f_meas),
  control = control_race(verbose_elim = TRUE)
)

print(select_best(xgb_rs))

best_mtry <- 5 
best_trees <- 408
best_min_n <- 18
best_treedepth <- 8
best_loss_red <- 3890
best_sample_size <- 0.681
best_stop_iter <- 23
best_learn_rate <- 1.42

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
  loss_reduction = best_loss_red,
  sample_size = best_sample_size,
  stop_iter = best_stop_iter,
  learn_rate = best_learn_rate,
  max_depth = best_treedepth
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
competition_data$covidDeath <- ifelse(competition_data$covidDeath == 1, "DIED", "SURVIVED")
competition_data$covidDeath <- as.factor(competition_data$covidDeath)

xgb_predictions <- predict(final_fit, new_data = competition_data)
xgb_predictions <- ifelse(predictions == "DIED", TRUE, FALSE)

submission <- data.frame(
  Patient_id = competition_data$Patient_id,
  covidDeath = predictions
)
write.csv(submission, file = "submission.csv", row.names = FALSE)

# model 7
# random forest, tuned with race

library(ranger)


prep(rec)
rf_spec <- rand_forest(
  mtry = tune(),
  trees = tune(),
  min_n = tune()
) %>%
  set_engine("ranger") %>%
  set_mode("classification")

rf_wf <- workflow(rec, rf_spec)

split_category <- function(x) {
  x %>%
    stringr::str_split(", ") %>%
    purrr::map(stringr::str_remove_all, "[:punct:]") %>%
    purrr::map(stringr::str_squish) %>%
    purrr::map(stringr::str_to_lower) %>%
    purrr::map(stringr::str_replace_all, " ", "_")
}

doParallel::registerDoParallel()

rf_grid <- grid_latin_hypercube(
  finalize(mtry(), train),
  trees(range = c(300, 1000)),
  min_n(range = c(2, 10)),
  size = 100                         
)

rf_rs <- tune_race_anova(
  rf_wf,
  resamples = folds,
  grid = rf_grid,
  metrics = metric_set(f_meas),
  control = control_race(verbose_elim = TRUE)
)

plot_race(rf_rs)

print(select_best(rf_rs))

best_mtry <- 15    
best_trees <- 1554
best_treedepth <- 5
best_min_n <- 34    

# make TRUE and FALSE cols into 0's and 1's for ranger
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
  learn_rate = 0.01
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

# make TRUE and FALSE cols into 0's and 1's for ranger
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

####### analysing predictions

predictions <- model4_predictions

test_data <- read.csv("CovidHospDataBrasil_test.csv", header = TRUE)

test_data$sex <- as.factor(test_data$sex)

for (i in 1:length(test_data)) {
  if (grepl("date", colnames(test_data)[i])) {
    test_data[,i] = as.Date(test_data[,i],
                                   format = "%d/%m/%Y")
  } else {
    if (i != 1 && i !=3 && i != 4) {
      test_data[,i] = as.logical(test_data[,i])
    }
  }
}

test_data$covidDeath <- predictions

View(test_data)

test_data$covidDeath

non_survivors <- test_data %>%
  filter(covidDeath) 

survivors <- test_data %>%
  filter(!covidDeath)

survivors$group <- "Survivors"
non_survivors$group <- "Non-Survivors"

combined_data <- rbind(survivors, non_survivors)

# visualise age

ggplot(combined_data, aes(x = group, y = age, fill = group)) +
  geom_boxplot() +
  labs(
    title = "Distribution of Age in Non-Survivors",
    x = "Group",
    y = "Age"
  ) +
  theme_minimal() +  
  theme(legend.position = "none")

summary(survivors$age)

summary(non_survivors$age)

# visualise admissions and deaths

aggregated_data <- non_survivors %>%
  group_by(dateHosp, group) %>%
  summarise(admissions = n(), .groups = 'drop')

ggplot(aggregated_data, aes(x = dateHosp, y = admissions, color = group)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Hospital Admissions Over Time for Non-Survivors",
    x = "Date of Hospitalization",
    y = "Number of Admissions"
  ) +
  theme_minimal()


test_data <- test_data %>%
  mutate(month_year = floor_date(dateHosp, "month"))

non_survivors <- non_survivors %>%
  mutate(month_year = floor_date(dateHosp, "month"))

total_patients_monthly <- test_data %>%
  group_by(month_year) %>%
  summarise(total_hospitalized = n())

# calculate the number of non-survivors each month
non_survivors_monthly <- non_survivors %>%
  group_by(month_year) %>%
  summarise(non_survivors = n())

# calculate the mortality rate
mortality_rate <- left_join(total_patients_monthly, non_survivors_monthly, by = "month_year") %>%
  mutate(non_survivors = replace_na(non_survivors, 0),
         mortality_rate = non_survivors / total_hospitalized)


ggplot(mortality_rate, aes(x = month_year, y = mortality_rate)) +
  geom_line(color = "red") +
  geom_point(color = "red") +
  labs(
    title = "Monthly Mortality Rate",
    x = "Month",
    y = "Mortality Rate"
  ) +
  scale_y_continuous(labels = scales::percent, limits = c(0.2,1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# adjust proportion of symptoms for total admissions, and make an ordered plot

total_counts <- combined_data %>%
  summarise(
    fever = sum(fever),
    cough = sum(cough),
    dyspnoea = sum(dyspnoea),
    sorethroat = sum(sorethroat),
    oxygensat = sum(oxygensat),
    diarrhea = sum(diarrhea),
    vomit = sum(vomit),
    hematologic = sum(hematologic),
    downsyn = sum(downsyn),
    asthma = sum(asthma),
    diabetes = sum(diabetes),
    neurological = sum(neurological),
    pneumopathy = sum(pneumopathy),
    obesity = sum(obesity),
    .groups = "drop"
  )

non_survivor_counts <- non_survivors %>%
  summarise(
    fever = sum(fever),
    cough = sum(cough),
    dyspnoea = sum(dyspnoea),
    sorethroat = sum(sorethroat),
    oxygensat = sum(oxygensat),
    diarrhea = sum(diarrhea),
    vomit = sum(vomit),
    hematologic = sum(hematologic),
    downsyn = sum(downsyn),
    asthma = sum(asthma),
    diabetes = sum(diabetes),
    neurological = sum(neurological),
    pneumopathy = sum(pneumopathy),
    obesity = sum(obesity)
  )

proportions <- non_survivor_counts / total_counts

symptom_proportions_long <- proportions %>%
  pivot_longer(cols = everything(), names_to = "symptom", values_to = "proportion")

symptom_proportions_long <- symptom_proportions_long %>%
  arrange(proportion) %>%
  mutate(symptom = factor(symptom, levels = symptom))


ggplot(symptom_proportions_long, aes(x = symptom, y = proportion)) +
  geom_bar(stat = "identity", fill = "#E86D6D") +
  labs(
    title = "Proportion of Non-Survivors with Each Symptom (Adjusted for Total Admissions)",
    x = "Symptom",
    y = "Proportion"
  ) +
  theme_minimal()  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 2 proportion z test on vaccinated proportions

vaccinated_survivors <- sum(survivors$vaccine)
total_survivors <- nrow(survivors)

vaccinated_non_survivors <- sum(non_survivors$vaccine)
total_non_survivors <- nrow(non_survivors)

prop.test(
  x = c(vaccinated_survivors, vaccinated_non_survivors),
  n = c(total_survivors, total_non_survivors),
  alternative = "greater"
)

vaccinated_patient_prop <- sum(test_data$vaccine) / nrow(test_data)

vaccinated_patient_prop

# show optimal threshold for f1

predictions <- predict(model2, test, type = "response")
f1_vec <- c()

for (i in seq(0.1, 0.8, 0.01)) {
  binary_predictions <- ifelse(predictions > i, TRUE, FALSE)
  
  confusion_matrix <- table(Predicted = binary_predictions, Actual = test$covidDeath)
  
  precision <- confusion_matrix[2,2] / (confusion_matrix[2,2] + confusion_matrix[2,1])
  recall <- confusion_matrix[2,2] / (confusion_matrix[2,2] + confusion_matrix[1,2])
  
  f1_vec <- c(f1_vec, 2 * precision * recall / (precision + recall))
}

glm_f1 <- data.frame(
  threshold = seq(0.1, 0.8, 0.01),
  f1 = f1_vec
)

plot(glm_f1)

glm_f1[which.max(glm_f1$f1),]

max_f1_index <- which.max(glm_f1$f1)
max_threshold <- glm_f1$threshold[max_f1_index]

# Create the ggplot
ggplot(glm_f1, aes(x = threshold, y = f1)) +
  geom_point(color = "blue") +
  geom_vline(xintercept = max_threshold, color = "red", linetype = "dashed") +
  labs(
    title = "F1 Score vs Threshold",
    x = "Threshold",
    y = "F1 Score"
  ) +
  theme_minimal() +
  annotate("text", x = max_threshold, y = min(glm_f1$f1), label = paste("Max threshold:", round(max_threshold, 2)), vjust = -1, color = "red", angle = 90, hjust = 1.5)
