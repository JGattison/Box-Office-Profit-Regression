# Load required libraries
library(dplyr)
library(DataExplorer)
library(ggplot2)
library(stringr)
library(lubridate)
library(tidyr)
library(lmtest)
library(MASS)
library(car)
library(sandwich)  # For robust GLM

options(scipen = 999)

# Data Cleaning and Merging
# --------------------------------------------------

# Merge IMDB, earnings, and genre datasets
box_office <- IMDB %>%
  inner_join(earning, by = "Movie_id") %>%
  left_join(genre, by = "Movie_id") %>%
  left_join(Movie_release_dates, by = "Movie_id")

# Handling genres with small counts
box_office <- box_office %>%
  group_by(genre) %>%
  mutate(genre = ifelse(n() < 10, NA, genre)) %>%
  ungroup()

# Create additional columns for date-related analysis
box_office$year <- str_sub(box_office$Title.x, -6) %>%
  gsub("[()]", "", .)
box_office$profit <- box_office$Domestic + box_office$Worldwide - box_office$Budget
box_office$Revenue <- box_office$Domestic + box_office$Worldwide
box_office$Profit_perc <- (box_office$profit / box_office$Budget) * 100

# Convert runtime to numeric
box_office$Runtime <- as.integer(substr(box_office$Runtime, 1, nchar(box_office$Runtime) - 4))

# Ensure proper date formatting
box_office$`Release Date_2` <- mdy(box_office$`Release Date`)

# Weekend Summary Cleaning
# --------------------------------------------------

# Function to clean and process weekend summary data for each year
clean_weekend_summary <- function(df, year) {
  df %>%
    mutate(date_2 = mdy(paste0(sub("-.*", "", date), ", ", year))) %>%
    group_by(date_2) %>%
    slice_max(order_by = overall_gross) %>%
    ungroup()
}

# Apply cleaning function to each year
weekend_summaries <- list(
  weekend_summary_2010 = clean_weekend_summary(weekend_summary_2010, 2010),
  weekend_summary_2011 = clean_weekend_summary(weekend_summary_2011, 2011),
  weekend_summary_2012 = clean_weekend_summary(weekend_summary_2012, 2012),
  # Repeat for each year...
  weekend_summary_2020 = clean_weekend_summary(weekend_summary_2020, 2020)
)

# Combine all weekend summaries
weekend_summary <- bind_rows(weekend_summaries)

# Adjust Release Dates to Fridays
# --------------------------------------------------
box_office <- box_office %>%
  mutate(
    `Release Date_3` = case_when(
      wday(`Release Date_2`) == 6 ~ `Release Date_2`,
      wday(`Release Date_2`) < 6 ~ `Release Date_2` + (6 - wday(`Release Date_2`)),
      TRUE ~ `Release Date_2` - (wday(`Release Date_2`) - 6)
    )
  )

# Expanding Data by Adding Rows for 4 Weeks Post Release
days_to_add <- c(0, 7, 14, 21)

bo_4weeks <- box_office %>%
  rowwise() %>%
  mutate(new_dates = list(`Release Date_3` + days_to_add)) %>%
  unnest(new_dates) %>%
  distinct(Movie_id, `Release Date_3`, new_dates)

# Merge expanded data with box office data
box_office <- box_office %>%
  left_join(bo_4weeks, by = c("Movie_id", "Release Date_3"))

# Summarize Weekend Box Office Data
# --------------------------------------------------
box_office_weekend <- box_office %>%
  left_join(
    weekend_summary %>%
      select(date_2, occasion, top10_gross, overall_gross, num_releases, top_release),
    by = c("new_dates" = "date_2")
  ) %>%
  distinct()

box_office_weekend_agg <- box_office_weekend %>%
  group_by(Movie_id, Title.y, `Release Date_3`) %>%
  summarise(
    top_bo_4weeks = ifelse(any(top_release %in% Title.y), 1, 0),
    debut_top = ifelse(first(top_release == Title.y), 1, 0),
    top10_gross_4weeks = sum(top10_gross, na.rm = TRUE),
    overall_gross_4weeks = sum(overall_gross, na.rm = TRUE),
    num_releases = first(num_releases),
    occasion = ifelse(any(!is.na(occasion)), 1, 0),
    release_month = first(month(`Release Date_3`))
  )

# Merge final aggregated weekend data
box_office_FINAL <- box_office %>%
  left_join(box_office_weekend_agg, by = c("Movie_id", "Title.y")) %>%
  distinct() %>%
  filter(!is.na(genre))

# Exploratory Data Analysis (EDA) 
# --------------------------------------------------

# Define a helper function to generate histograms for each variable
plot_histogram <- function(data, feature, title, binwidth = NULL, x_label = feature, y_label = "Count") {
  ggplot(data, aes_string(x = feature)) +
    geom_histogram(binwidth = binwidth, color = 'black', fill = 'skyblue') +
    labs(title = title, x = x_label, y = y_label) +
    theme_minimal()
}

# Generate histograms for key variables
plot_histogram(box_office_FINAL, "profit", "Profit Distribution", binwidth = 100000000)
plot_histogram(box_office_FINAL, 'VotesM', 'Male Rating Distribution', binwidth = 0.1)
plot_histogram(box_office_FINAL, 'VotesF', 'Female Rating Distribution', binwidth = 0.1)

# Model Building and Evaluation
# --------------------------------------------------

# GLM with Gamma distribution model 
colnames(box_office_FINAL_binary_encoded)

# 1. Define a Function to Fit a Model and Get Robust Standard Errors
fit_glm_with_robust_se <- function(formula, data) {
  model <- glm(formula, family = gaussian(link = "identity"), data = data)
  robust_se <- coeftest(model, vcov = vcovHC(model, type = "HC3"))
  list(model = model, robust_se = robust_se)
}

# 2. Function to Perform Assumption Checks
assumption_checks <- function(model) {
  # Linearity check
  plot(fitted(model), residuals(model), 
       main = "Residuals vs Fitted Values",
       xlab = "Fitted Values", ylab = "Residuals")
  abline(h = 0, col = "red")
  
  # Autocorrelation test
  print(dwtest(model))
  
  # Multicollinearity test
  print(vif(model))
  
  # Heteroscedasticity check
  print(bptest(model))
  
  # Cook's Distance for influential observations
  cooks_dist <- cooks.distance(model)
  plot(cooks_dist, main = "Cook's Distance")
}

# 3. Define your main formula and data
formula_main <- profit ~ Rating + TotalVotes + Budget + Runtime + MetaCritic + 
  year + release_month + top_bo_4weeks + debut_top + 
  top10_gross_4weeks + overall_gross_4weeks + num_releases + 
  occasion + Action + Adventure + Animation + Biography + 
  Comedy + Crime + Drama + Romance + `Sci-Fi` + Thriller + VotesM + VotesF

data <- box_office_FINAL_binary_encoded_clean

# 4. Fit the main GLM and get robust SEs
result_main <- fit_glm_with_robust_se(formula_main, data)
summary(result_main$model)

# 5. Perform assumption checks on the main model
assumption_checks(result_main$model)

# 6. Stepwise AIC Model Selection
stepwise_glm <- stepAIC(result_main$model, direction = "both")

# 7. Check for Robust SEs and Assumptions on the Stepwise Model
robust_se_stepwise <- coeftest(stepwise_glm, vcov = vcovHC(stepwise_glm, type = "HC3"))
summary(stepwise_glm)
robust_se_stepwise

# Perform assumption checks on stepwise model
assumption_checks(stepwise_glm)

# Given the heteroscedasticity, we need to identify the problem
plot(fitted(stepwise_glm), residuals(model, type = "pearson"),
     xlab = "Fitted values", ylab = "Pearson residuals",
     main = "Residuals vs Fitted")

# Predicted values at the lowest profit are more clustered together and at highest profit have higher variance of residuals44

# Given this we use Robust Standard Errors. 

# Final Robust Model with Selected Features
final_formula <- profit ~ TotalVotes + Budget + Runtime + MetaCritic + 
  overall_gross_4weeks + num_releases + VotesF
final_result <- fit_glm_with_robust_se(final_formula, data)
summary(final_result$model)
final_result$robust_se

# 10. Perform Assumption Checks on Final Model
assumption_checks(final_result$model)

# Export data

export <- box_office_FINAL_binary_encoded %>% 
  dplyr::select(profit, TotalVotes, Budget, MetaCritic, overall_gross_4weeks, VotesM, VotesF, Rating)

write.csv(export, file = "export.csv")
