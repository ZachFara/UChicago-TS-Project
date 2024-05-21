### ADSP 31006 IP02 Time Series Analysis and Forecasting Final Project ###
### Daichi Ishikawa ###

# Install library
library(ggplot2)
library(readxl)
library(lubridate)
library(hms)
library(dplyr)
library(zoo)
library(corrplot)
library(lmtest)
library(stats)
library(tseries)
library(Metrics)
library(forecast)
library(xts)
library(TSA)

#######################
#### Preprocessing ####
#######################

# Load data
data <- read_excel("data/AirQualityUCI.xlsx")
data <- data[, c("Date", "Time", "CO(GT)", "T", "RH", "AH")]

# Date time Format Change
data$Date <- as.Date(data$Date)
data$Time <- hms::as_hms(data$Time)
data$DateTime <- ymd_hms(paste(data$Date, data$Time))
data <- data[, c("DateTime", "CO(GT)", "T", "RH", "AH")]

# Function to replace -200 with NA, perform linear interpolation, and set negative values to 0
replace_anomalies_and_fill <- function(column) {
  # Replace -200 with NA
  column[column == -200] <- NA
  # Perform linear interpolation to replace NA values
  column <- na.approx(column, na.rm = FALSE, rule = 2)
  # Replace negative values with 0
  column <- pmax(column, 0)
  return(column)
}

# Apply the function to the relevant columns
data <- data %>%
  mutate(
    `CO(GT)` = replace_anomalies_and_fill(`CO(GT)`),
    T = replace_anomalies_and_fill(T),
    RH = replace_anomalies_and_fill(RH),
    AH = replace_anomalies_and_fill(AH)
  )

head(data)

#######################
######## EDA ##########
#######################

## Plot histograms for each variable
# Function to create a histogram for a given variable
create_histogram <- function(data, variable, binwidth, fill_color, title, x_label) {
  ggplot(data, aes(x = .data[[variable]])) +
    geom_histogram(binwidth = binwidth, fill = fill_color, color = "black", alpha = 0.7) +
    labs(title = title, x = x_label, y = "Frequency")
}
# Create histograms
create_histogram(data, "CO(GT)", 0.1, "blue", "Histogram of CO(GT)", "CO(GT)")
create_histogram(data, "T", 1, "red", "Histogram of Temperature (T)", "Temperature (T)")
create_histogram(data, "RH", 1, "green", "Histogram of Relative Humidity (RH)", "Relative Humidity (RH)")
create_histogram(data, "AH", 0.01, "purple", "Histogram of Absolute Humidity (AH)", "Absolute Humidity (AH)")

## Time Series Plot
# Function to create a time series plot for a given variable
create_time_series <- function(data, datetime_var, variable, line_color, title, y_label) {
  ggplot(data, aes(x = .data[[datetime_var]], y = .data[[variable]])) +
    geom_line(color = line_color) +
    labs(title = title, x = "DateTime", y = y_label)
}
# Create time series plots
create_time_series(data, "DateTime", "CO(GT)", "blue", "Time Series of CO(GT)", "CO(GT)")
create_time_series(data, "DateTime", "T", "red", "Time Series of Temperature (T)", "Temperature (T)")
create_time_series(data, "DateTime", "RH", "green", "Time Series of Relative Humidity (RH)", "Relative Humidity (RH)")
create_time_series(data, "DateTime", "AH", "purple", "Time Series of Absolute Humidity (AH)", "Absolute Humidity (AH)")

## Correlation
# Select relevant columns for correlation matrix
correlation_data <- data %>% select(`CO(GT)`, T, RH, AH)
# Calculate the correlation matrix
correlation_matrix <- cor(correlation_data, use = "complete.obs")
# Display the correlation matrix
print(correlation_matrix)
# Visualize the correlation matrix
corrplot(correlation_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, addCoef.col = "black")

## Stationality
# ADF test for each variable
adf_test_CO <- adf.test(data$`CO(GT)`, alternative = "stationary")
# adf_test_T <- adf.test(data$T, alternative = "stationary")
# adf_test_RH <- adf.test(data$RH, alternative = "stationary")
# adf_test_AH <- adf.test(data$AH, alternative = "stationary")

# Display the results of ADF test
print(adf_test_CO)
# print(adf_test_T)
# print(adf_test_RH)
# print(adf_test_AH)

# KPSS test for each variable
kpss_test_CO <- kpss.test(data$`CO(GT)`)
# kpss_test_T <- kpss.test(data$T)
# kpss_test_RH <- kpss.test(data$RH)
# kpss_test_AH <- kpss.test(data$AH)

# Display the results of KPSS test
print(kpss_test_CO)
# print(kpss_test_T)
# print(kpss_test_RH)
# print(kpss_test_AH)

## Autocorrelation
# Plot ACF
acf(data$`CO(GT)`, main="ACF of CO(GT)")
# Plot PACF
pacf(data$`CO(GT)`, main="PACF of CO(GT)")


## HoltWinters
co_ts <- ts(data$`CO(GT)`, frequency = 24)
hw_model <- HoltWinters(co_ts)
print(hw_model)
plot(hw_model)

## STL
co_ts <- ts(data$`CO(GT)`, frequency = 24)
stl_decomposition <- stl(co_ts, s.window = "periodic")
plot(stl_decomposition)


#######################
### Train/Test Split ##
#######################

# Convert DateTime column to POSIXct type if not already
data$DateTime <- as.POSIXct(data$DateTime)

# Find the maximum date in the data
max_date <- max(data$DateTime)

# Define the cutoff date as 1 month before the maximum date
cutoff_date <- max_date - months(1)

# Split the data into training and test sets
train_data <- data %>% filter(DateTime < cutoff_date)
test_data <- data %>% filter(DateTime >= cutoff_date)

# Display the first few rows of each dataset to verify the split
head(train_data)
head(test_data)

########################
## Model and Eveluate ##
########################

## Periodogram
co_ts <- ts(train_data$`CO(GT)`, frequency = 24)
# Calculate the periodogram
p <- periodogram(co_ts)
# Display the periodogram result
print(p)

# Number of top peaks to consider
top_n <- 10
# Find the frequencies with the highest spectral densities
top_freq_indices <- order(p$spec, decreasing = TRUE)[1:top_n]
top_freqs <- p$freq[top_freq_indices]
# Calculate the corresponding seasonal periods
seasonality_periods <- 1 / top_freqs
# Print the detected seasonal periods
print("Detected seasonality periods:")
print(seasonality_periods)


## TBATS

# Creating MSTS objects
train_msts <- msts(train_data$`CO(GT)`, seasonal.periods = c(24, 24*7))

# Creating and fitting TBATS model
tbats_model <- tbats(train_msts)
# Plot the model
plot(tbats_model)
# Displaying model summary
print(as.character(tbats_model))

# fit train
train_fitted <- fitted(tbats_model)

# Prediction
h <- nrow(test_data)  # Number of steps to forecast
test_msts <- msts(test_data$`CO(GT)`, seasonal.periods = c(24, 24*7))
forecast_tbats <- forecast(tbats_model, h = h)

# Plot forcast
plot(forecast(tbats_model, h=h))

# Calculate residuals for training data
train_residuals <- train_data$`CO(GT)` - train_fitted
# Plot residuals for training period
plot(train_data$DateTime, train_residuals, type = "l", col = "purple", xlab = "Time", ylab = "Residuals", main = "Residuals for Training Period")
abline(h = 0, col = "red", lty = 2)

# Plotting for train period
plot(train_data$DateTime, train_data$`CO(GT)`, type = "l", col = "black", xlab = "Time", ylab = "CO(GT)", main = "Model Fit for Training Period")
lines(train_data$DateTime, train_fitted, col = "blue", lty = 1)
legend("topright", legend = c("Actual", "Fitted"), col = c("black", "blue"), lty = c(1, 1))

# Plotting for test period
plot(test_data$DateTime, test_data$`CO(GT)`, type = "l", col = "black", xlab = "Time", ylab = "CO(GT)", main = "Forecast vs Actual for Test Period")
lines(test_data$DateTime, forecast_tbats$mean, col = "red", lty = 1)
legend("topright", legend = c("Actual", "Forecast"), col = c("black", "red"), lty = c(1, 1))

# Calculating RMSE and MAE
rmse_value <- rmse(test_data$`CO(GT)`, forecast_tbats$mean)
mae_value <- mae(test_data$`CO(GT)`, forecast_tbats$mean)
print(paste("RMSE:", rmse_value))
print(paste("MAE:", mae_value))


# Show Parameters
lambda <- tbats_model$lambda
alpha <- tbats_model$alpha
beta <- tbats_model$beta
damping_parameter <- tbats_model$damping.parameter
gamma_one_values <- tbats_model$gamma.one.values
gamma_two_values <- tbats_model$gamma.two.values
ar_coefficients <- tbats_model$ar.coefficients
ma_coefficients <- tbats_model$ma.coefficients
likelihood <- tbats_model$likelihood
optim_return_code <- tbats_model$optim.return.code
variance <- tbats_model$variance
AIC <- tbats_model$AIC
parameters <- tbats_model$parameters
seed_states <- tbats_model$seed.states
fitted_values <- tbats_model$fitted.values
errors <- tbats_model$errors
x <- tbats_model$x
seasonal_periods <- tbats_model$seasonal.periods
k_vector <- tbats_model$k.vector
y <- tbats_model$y
p <- tbats_model$p
q <- tbats_model$q
call <- tbats_model$call
series <- tbats_model$series
method <- tbats_model$method

print(paste("lambda:", lambda))
print(paste("alpha:", alpha))
print(paste("beta:", beta))
print(paste("damping.parameter:", damping_parameter))
print(paste("gamma.one.values:", paste(gamma_one_values, collapse = ", ")))
print(paste("gamma.two.values:", paste(gamma_two_values, collapse = ", ")))
print(paste("ar.coefficients:", ifelse(length(ar_coefficients) == 0, "NULL", paste(ar_coefficients, collapse = ", "))))
print(paste("ma.coefficients:", ifelse(length(ma_coefficients) == 0, "NULL", paste(ma_coefficients, collapse = ", "))))
print(paste("likelihood:", likelihood))
print(paste("variance:", variance))
print(paste("AIC:", AIC))
print(paste("parameters:", paste(parameters, collapse = ", ")))
print(paste("seasonal.periods:", paste(seasonal_periods, collapse = ", ")))
print(paste("k.vector:", paste(k_vector, collapse = ", ")))
print(paste("p:", p))
print(paste("q:", q))
