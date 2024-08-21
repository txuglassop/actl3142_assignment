library(dplyr)
library(pROC)
library(ggplot2)
library(AMR)
library(bestglm)
library(gridExtra)
library(MASS)
library(broom)
library(car)
library(caret)

# read data
data <- read.csv("CovidHospDataBrasil.csv", header = TRUE)

# clean data and get a vector of all Symptoms
data$sex <- as.factor(data$sex)

Symptoms <- c()

for (i in 1:length(data)) {
    if (grepl("date", colnames(data)[i])) {
        data[,i] = as.Date(data[,i])
    } else {
        if (i != 1 && i !=3 && i != 4) {
            data[,i] = as.logical(data[,i])
            Symptoms <- c(Symptoms, names(data)[i])
        }
    }
}

# the above method counted "vaccine" "icu" "covidDeath" as Symptoms
# remove them
Symptoms <- Symptoms[-c(1,16,22)]

# add more variables of interest
data <- transform(
    data, 
    "daysInIcu" = as.numeric(dateDisIcu - dateAdmIcu),
    "daysInHosp" = as.numeric(dateEndObs - dateHosp),
    "daysBeforeIcu" = as.numeric(dateAdmIcu - dateHosp),
    "icuAndDeath" = ifelse(icu & covidDeath, TRUE, FALSE),
    "vacAndDeath" = ifelse(vaccine & covidDeath, TRUE, FALSE),
    "noVacAndDeath" = ifelse(!vaccine & covidDeath, TRUE, FALSE),
    "icuAndVac" = ifelse(vaccine & icu, TRUE, FALSE),
    "icuAndNoVac" = ifelse(!vaccine & icu, TRUE, FALSE),
    "numSymptoms" = rowSums(data[, Symptoms], na.rm = TRUE)
)

##### EDA

### Vaccinated / Unvaccinated Profile
### - Number of admissions of both types
### - icu / mortality rate of both types

# count admissions, deaths and icu admissions of vaccinated and not
num_vac <- data.frame(
    Status = c("Vaccinated", "Unvaccinated"),
    Count = c(
        sum(data$vaccine == TRUE),
        sum(data$vaccine == FALSE)
    )
)

vac_death <- data.frame(
    Status = c("Vaccinated", "Unvaccinated"),
    Count = c(
        sum(data$vacAndDeath == TRUE),
        sum(data$noVacAndDeath == TRUE)
    )
)

vac_icu <- data.frame(
    Status = c("Vaccinated", "Unvaccinated"),
    Count = c(
        sum(data$icuAndVac == TRUE),
        sum(data$icuAndNoVac == TRUE)
    )
)

# plot vaccinated vs unvaccinated covid deaths and icu admissions

neutral_colors <- c("Vaccinated" = "#6baed6", "Unvaccinated" = "#fd8d3c")
y_limit <- c(0, 100000)

p1 <- ggplot(num_vac, aes(x = Status, y = Count, fill = Status)) +
  geom_bar(stat = "identity", width = 0.3) +
  labs(
    title = "Total Admissions",
    x = "Status",
    y = "Number of Admissions"
  ) +
  scale_fill_manual(values = neutral_colors) +
  coord_cartesian(ylim = y_limit) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.position = "none"
  )

p2 <- ggplot(vac_icu, aes(x = Status, y = Count, fill = Status)) +
  geom_bar(stat = "identity", width = 0.3) +
  labs(
    title = "ICU Admissions",
    x = "Status",
    y = "Number of ICU Admissions"
  ) +
  scale_fill_manual(values = neutral_colors) +
  coord_cartesian(ylim = y_limit) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.position = "none"
  )

p3 <- ggplot(vac_death, aes(x = Status, y = Count, fill = Status)) +
  geom_bar(stat = "identity", width = 0.3) +
  labs(
    title = "Deaths from Covid",
    x = "Status",
    y = "Number of Deaths"
  ) +
  scale_fill_manual(values = neutral_colors) +
  coord_cartesian(ylim = y_limit) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.position = "none"
  )

# Arrange the plots side by side
grid.arrange(p1, p2, p3, ncol = 3)

# Perform the two-proportion z-test on the 
vac_death_test <- prop.test(
  x = c(vac_death[1,2], vac_death[2,2]),
  n = c(num_vac[1,2], num_vac[2,2]), 
  alternative = "less", 
  correct = FALSE 
)

vac_icu_test <- prop.test(
  x = c(vac_icu[1,2], vac_icu[2,2]),
  n = c(num_vac[1,2], num_vac[2,2]), 
  alternative = "less", 
  correct = FALSE 
)

vac_death_test
vac_icu_test


### Age profile of survivors and non-survivors
### - Show the age distribution of all admissions
### - Investigate the mortality rate of each age

# Find age distribution and show survivors and non-survivors

data$survivalStatus <- ifelse(data$covidDeath == TRUE, "Non-Survivor", "Survivor")

ggplot(data, aes(x = age, fill = survivalStatus)) +
    geom_histogram(binwidth = 2, color = "black", position = "stack") +
    labs(
        title = "Age Distribution of COVID-19 Admissions with Mortality",
        x = "Age",
        y = "Count",
        fill = "Survival Status"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(size =15, hjust = 0.5)
    )

# Investigate mortality rate of each age group
data$ageGroup <- age_groups(data$age, split_at = seq(14, 90, by = 2))

# Calculate mortality rate for each age group
mortality_rate <- data %>%
    group_by(ageGroup) %>%
    summarize(
        Total = n(),
        Deaths = sum(covidDeath),
        MortalityRate = Deaths / Total * 100
    )


# Plot the mortality rate for each age group

ggplot(mortality_rate, aes(x = ageGroup, y = MortalityRate, group = 1)) +
    geom_point(colour = "black", size = 1.5) +
	geom_smooth(method = "lm", colour = "red", size = 1) + 
    labs(
        title = "Mortality Rate by Age Group",
        x = "Age Group",
        y = "Mortality Rate (%)"
    ) +
    theme_minimal() +
	scale_x_discrete(
        breaks = levels(mortality_rate$ageGroup)[seq(1, length(levels(mortality_rate$ageGroup)), by = 3)]
    ) +
    theme(
        plot.title = element_text(size = 12, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
    )

### Symptoms Profile of Patientes
### - Number of Symptoms of survivors and non-survivors

symptoms_summary_data <- data %>%
    group_by(covidDeath) %>%
    summarise(
		count = n(),
		mean_symptoms = mean(numSymptoms, na.rm = TRUE),
		median_symptoms = median(numSymptoms, na.rm = TRUE),
		min_symptoms = min(numSymptoms, na.rm = TRUE),
		max_symptoms = max(numSymptoms, na.rm = TRUE),
		sd_symptoms = sd(numSymptoms, na.rm = TRUE)
    )

symptoms_summary_data

# perform t test for test of significant difference
non_survivors <- data %>%
	filter(covidDeath) 

survivors <- data %>%
	filter(!covidDeath)

t.test(
	x = survivors$numSymptoms, y = non_survivors$numSymptoms,
	alternative = "less",
	mu =  0
)

# plot the distribution of the number of symptoms of survivors and non-survivors

survivors$group <- "survivors"
non_survivors$group <- "nonsurvivors"

combined_data <- rbind(survivors, non_survivors)

ggplot(combined_data, aes(x = group, y = numSymptoms)) +
	geom_boxplot(alpha = 0.5) +
	labs(title = "Number of Symptoms of Survivors and Non-Survivors",
		x = "Number of Symptoms",
		y = "Density") +
	theme_minimal()


### Duration of Hospital Stay

stay_summary_data <- data %>%
    group_by(covidDeath) %>%
    summarise(
		count = n(),
		mean_stay = mean(daysInHosp, na.rm = TRUE),
		median_stay = median(daysInHosp, na.rm = TRUE),
		min_stay = min(daysInHosp, na.rm = TRUE),
		max_stay = max(daysInHosp, na.rm = TRUE),
		sd_stay = sd(daysInHosp, na.rm = TRUE)
    )

stay_summary_data


# perform t test for test of significant difference
t.test(
	x = survivors$daysInHosp, y = non_survivors$daysInHosp,
	alternative = "less",
	mu =  0
)

# plot hospital stay of survivors and nonsurvivors
combined_data$daysInHospCapped <- ifelse(combined_data$daysInHosp > 100, 100, combined_data$daysInHosp)

ggplot(combined_data, aes(x = daysInHospCapped, fill = group)) +
	geom_histogram(aes(y = ..density.. ), binwidth = 0.5, alpha = 0.5, position = "identity") +
	labs(title = "Days Spent in Hospital of Survivors and Non-Survivors",
		fill = "Group",
		x = "Days",
		y = "Density") +
	theme_minimal() +
	scale_fill_manual(
		values = c("red", "blue"),
		labels = c("Non-Survivors", "Survivors")
	)



#######################
##### Fitting GLM #####
#######################

### Fit naive model

naive_model <- glm(
  covidDeath ~ age + sex + vaccine + fever + cough + sorethroat + hematologic +
  dyspnoea + oxygensat + diarrhea + vomit + asthma + diabetes + downsyn +
  neurological + pneumopathy + obesity + icu + respdistress + cardio +
  hepatic + immuno + renal,
  data = data,
  family = binomial()
)

summary(naive_model)

### Fit improved model, removing insignificant variables and adding
### extra variables that can be inferred from the data

scope <- list(
	lower = covidDeath ~ 1,
	upper = covidDeath ~ age + sex + vaccine + fever + cough + sorethroat + hematologic +
  dyspnoea + oxygensat + diarrhea + vomit + asthma + diabetes + downsyn +
  neurological + pneumopathy + obesity + icu + respdistress + cardio +
  hepatic + immuno + renal + daysInHosp
)

improved_model <- stepAIC(
	object = naive_model,
	scope = scope,
	direction = "both"	
)

# change out age for date of birth and compare the models
model_with_dob <- glm(
	covidDeath ~ dateBirth + sex + vaccine + fever + cough + dyspnoea + oxygensat +
	diarrhea + asthma + diabetes + neurological + pneumopathy + obesity + icu +
	downsyn + hematologic + respdistress + hepatic + immuno + renal + daysInHosp,
	data = data,
	family = binomial()
)

glance(improved_model)
glance(model_with_dob)

# model with dob is better.
improved_model <- model_with_dob

summary(improved_model)

# plot ROC curves and find AUC
improved_predicted_odds <- predict(improved_model, type = "response")

improved_roc_curve <- roc(data$covidDeath, improved_predicted_odds)

improved_tpr <- improved_roc_curve$sensitivities
improved_fpr <- 1 - improved_roc_curve$specificities
naive_predicted_odds <- predict(naive_model, type = "response")

naive_roc_curve <- roc(data$covidDeath, naive_predicted_odds)

naive_tpr <- naive_roc_curve$sensitivities
naive_fpr <- 1 - naive_roc_curve$specificities

naive_auc_value <- auc(naive_roc_curve)
improved_auc_value <- auc(improved_roc_curve)

# Plot the ROC curves 
plot(
  naive_fpr, naive_tpr, type = "l",
  main = "ROC Curves of Naive and Improved Models", col = "blue", lwd = 2,
  xlab = "False Positive Rate",
  ylab = "True Positive Rate"
)

lines(improved_fpr, improved_tpr, col = "green", lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "red")
legend("bottomright", 
       legend = c(paste("Naive AUC =", round(naive_auc_value, 2)), 
                  paste("Improved AUC =", round(improved_auc_value, 2))), 
       col = c("blue", "green"), lwd = 2)



### Order coefficients by magnitude
sort(abs(improved_model$coefficients))


### Model comparison

# compared the models
glance(naive_model)
glance(improved_model)

car::vif(naive_model)
car::vif(improved_model)

# create confusion matrices for both models
# split the data 80 - 20 for training and testing
index <- createDataPartition(data$covidDeath, p = 0.8, list = FALSE)
train_data <- data[index,]
test_data <- data[-index,]

naive_model_2 <- glm(
  covidDeath ~ age + sex + vaccine + fever + cough + sorethroat +
  dyspnoea + oxygensat + diarrhea + vomit + asthma + diabetes + 
  neurological + pneumopathy + obesity + icu + respdistress + cardio +
  hepatic + immuno + renal + downsyn + hematologic,
  data = train_data,
  family = binomial()
)

improved_model_2 <- glm(
	covidDeath ~ dateBirth + sex + vaccine + fever + cough + dyspnoea + oxygensat +
	diarrhea + asthma + diabetes + neurological + pneumopathy + obesity + icu +
	downsyn + hematologic + respdistress + hepatic + immuno + renal + daysInHosp,
	data = train_data,
	family = binomial()
)

threshold <- 0.5

naive_predictions <- predict(
	naive_model_2, test_data, type = "response"
)

naive_final <- ifelse(
	naive_predictions > threshold,
	1, 0
)

improved_predictions <- predict(
	improved_model_2, test_data, type = "response"
)

improved_final <- ifelse(
	improved_predictions > threshold,
	1, 0
)

naive_prediction_matrix <- table(naive_final, test_data$covidDeath)
naive_tpr <- naive_prediction_matrix[2,2] / sum(naive_prediction_matrix[,2])
naive_fpr <- naive_prediction_matrix[2,1] / sum(naive_prediction_matrix[,1])

improved_prediction_matrix <- table(improved_final, test_data$covidDeath)
improved_tpr <- improved_prediction_matrix[2,2] / sum(improved_prediction_matrix[,2])
improved_fpr <- improved_prediction_matrix[2,1] / sum(improved_prediction_matrix[,1])

data.frame(
	Model = c("Naive", "Improved"),
	TPR = c(naive_tpr, improved_tpr),
	FPR = c(naive_fpr, improved_fpr)
)

# Aggregate the number of deaths per day
deaths_per_day <- non_survivors %>%
  group_by(dateEndObs) %>%
  summarise(deaths = n())

# Plot deaths over time
ggplot(deaths_per_day, aes(x = dateEndObs, y = deaths)) +
  geom_line(color = "red") +
  geom_point(color = "red") +
  labs(
    title = "Daily COVID Deaths",
    x = "Date",
    y = "Number of Deaths"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Aggregate hospital admissions by date
admissions_by_date <- data %>%
  group_by(dateHosp) %>%
  summarize(Admissions = n()) %>%
  rename(Date = dateHosp)

# Aggregate deaths by date
deaths_by_date <- non_survivors %>%
  group_by(dateEndObs) %>%
  summarize(Deaths = n()) %>%
  rename(Date = dateEndObs)

# Combine both datasets
combined_data <- full_join(admissions_by_date, deaths_by_date, by = "Date") %>%
  replace_na(list(Admissions = 0, Deaths = 0))

# Pivot the data for easier plotting
combined_data_long <- combined_data %>%
  pivot_longer(cols = c("Admissions", "Deaths"), names_to = "Metric", values_to = "Count")

# Plot the data
ggplot(combined_data_long, aes(x = Date, y = Count, color = Metric)) +
  geom_line(size = 1) +
  labs(
    title = "Hospital Admissions and Deaths Over Time",
    x = "Date",
    y = "Count",
    color = "Metric"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Aggregate hospital admissions by date
admissions_by_date <- data %>%
  group_by(dateHosp) %>%
  summarize(Admissions = n()) %>%
  rename(Date = dateHosp)

# Aggregate deaths by date
deaths_by_date <- non_survivors %>%
  group_by(dateEndObs) %>%
  summarize(Deaths = n()) %>%
  rename(Date = dateEndObs)

# Combine both datasets
combined_data <- full_join(admissions_by_date, deaths_by_date, by = "Date") %>%
  replace_na(list(Admissions = 0, Deaths = 0))

# Calculate mortality rate
combined_data <- combined_data %>%
  mutate(MortalityRate = ifelse(Admissions > 0, Deaths / Admissions * 100, 0))

# Determine the date to exclude the last month
last_month_start <- floor_date(max(combined_data$Date), unit = "month")

# Filter out the last month's data
filtered_data <- combined_data %>%
  filter(Date < last_month_start)

# Pivot the data for easier plotting
filtered_data_long <- filtered_data %>%
  pivot_longer(cols = c("Admissions", "Deaths", "MortalityRate"), names_to = "Metric", values_to = "Count")

# Plot the data
ggplot(filtered_data_long, aes(x = Date, y = Count, color = Metric)) +
  geom_line(size = 1) +
  labs(
    title = "Hospital Admissions, Deaths, and Mortality Rate Over Time",
    x = "Date",
    y = "Count / Mortality Rate (%)",
    color = "Metric"
  ) +
  scale_y_continuous(
    sec.axis = sec_axis(~ ., name = "Mortality Rate (%)")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
