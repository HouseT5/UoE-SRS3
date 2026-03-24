# ==============================================================================
# MATH11188 Statistical Research Skills - Assignment 2
# University of Edinburgh, MSc Statistics with Data Science
# ==============================================================================

# SECTION 1: Load Libraries
library(ncdf4)      # Reading NetCDF files
library(dplyr)      # Data manipulation
library(ggplot2)    # Plotting
library(MASS)       # Negative Binomial GLM (glm.nb)
library(AER)        # Formal overdispersion test (dispersiontest)
library(DHARMa)     # Simulation-based absolute goodness-of-fit diagnostics
library(tibble)     # Modern data frames
library(tidyr)      # Data reshaping
library(lubridate)  # Date handling
library(patchwork)  # Combining ggplot figures

# SECTION 2: Load NetCDF Data
 setwd("C:/Users/toddh/SRS/temp_data_SRS.nc") # Ensure this points to your file location
#nc <- nc_open("temp_data_SRS.nc")
#lon <- ncvar_get(nc, "lon")
#lat <- ncvar_get(nc, "lat")
#time <- ncvar_get(nc, "time")
#tmp <- ncvar_get(nc, "tmp")
#nc_close(nc)
df_raw <- read.csv("uk_temperature_data.csv")
df_raw$date <- as.Date(df_raw$time, format = "%m/%d/%Y")
#dates <- as.Date(time, origin = "1900-01-01")

# SECTION 3: Spatial Subsetting - UK Bounding Box
#lon_uk_idx <- which(lon >= -8.5 & lon <= 2.0)
#lat_uk_idx <- which(lat >= 49.5 & lat <= 61.0)
#tmp_uk <- tmp[lon_uk_idx, lat_uk_idx, ]
df_raw <- df_raw %>% rename(lon = lon, lat = lat, tmp = tmp)

# SECTION 4: Compute UK Monthly Mean Temperature
#n_months <- dim(tmp_uk)[3]
#uk_mean_tmp <- sapply(1:n_months, function(t) mean(tmp_uk[,,t], na.rm = TRUE))
#df_monthly <- tibble(date = dates, year = year(dates), month = month(dates), temperature = uk_mean_tmp)

df_monthly <- df_raw %>%
  group_by(date) %>%
  summarise(temperature = mean(tmp, na.rm = TRUE), .groups = "drop") %>%
  mutate(year = year(date), month = month(date))


# SECTION 5: Define Extreme Events - WMO 1961-1990 Baseline
baseline <- df_monthly %>% filter(year >= 1961 & year <= 1990)
thresh_hot <- quantile(baseline$temperature, 0.95, na.rm = TRUE)
thresh_cold <- quantile(baseline$temperature, 0.05, na.rm = TRUE)

df_monthly <- df_monthly %>%
  mutate(
    extreme_warm = as.integer(temperature > thresh_hot),
    extreme_cold = as.integer(temperature < thresh_cold),
    extreme_any = as.integer(extreme_warm == 1 | extreme_cold == 1)
  )

# SECTION 6: Annual Aggregation
df_yearly <- df_monthly %>%
  group_by(year) %>%
  summarise(
    count_warm = sum(extreme_warm, na.rm = TRUE),
    count_cold = sum(extreme_cold, na.rm = TRUE),
    count_any = sum(extreme_any, na.rm = TRUE),
    mean_temp = mean(temperature, na.rm = TRUE),
    n_months = n(),
    .groups = "drop"
  ) %>%
  filter(n_months == 12)

# SECTION 7: Exploratory Data Visualisation
df_plot <- df_yearly %>%
  dplyr::select(year, count_warm, count_cold) %>%
  pivot_longer(cols = c(count_warm, count_cold), names_to = "type", values_to = "count") %>%
  mutate(type = ifelse(type == "count_warm", "Extreme Warm", "Extreme Cold"))

fig1 <- ggplot(df_plot, aes(x = year, y = count, fill = type)) +
  geom_col(position = "stack", width = 0.9) +
  scale_fill_manual(values = c("Extreme Warm" = "#d73027", "Extreme Cold" = "#4575b4")) +
  geom_smooth(data = df_yearly, aes(x = year, y = count_any), inherit.aes = FALSE, method = "loess", span = 0.3, linetype = "dashed", se = FALSE) +
  theme_bw()

fig2 <- ggplot(df_monthly, aes(x = date, y = temperature)) +
  geom_line(colour = "grey60", linewidth = 0.3, alpha = 0.6) +
  geom_smooth(method = "loess", span = 0.1, colour = "black", linewidth = 1.2, se = FALSE) +
  geom_hline(yintercept = thresh_hot, linetype = "dashed", colour = "#d73027", linewidth = 0.8) +
  geom_hline(yintercept = thresh_cold, linetype = "dashed", colour = "#4575b4", linewidth = 0.8) +
  annotate("text", x = min(df_monthly$date), y = thresh_hot + 0.6,
           label = paste0("Warm threshold: ", round(thresh_hot, 1), "°C"),
           colour = "#d73027", hjust = 0, size = 3.2, fontface = "italic") +
  annotate("text", x = min(df_monthly$date), y = thresh_cold + 0.6,
           label = paste0("Cold threshold: ", round(thresh_cold, 1), "°C"),
           colour = "#4575b4", hjust = 0, size = 3.2, fontface = "italic") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_x_date(date_breaks = "20 years", date_labels = "%Y") +
  labs(
    x = "Date",
    y = "Temperature (°C)",
    title = "UK Monthly Mean Temperature, 1901–2012",
    subtitle = "Grey lines = monthly values; black curve = LOESS trend; dashed lines = WMO baseline thresholds"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, colour = "grey40")
  )

fig2
# SECTION 8: Preliminary Overdispersion Assessment
cat("Var/Mean ratio:", var(df_yearly$count_any) / mean(df_yearly$count_any), "\n")

# SECTION 9: Model Fitting
fit_poisson <- glm(count_any ~ year + offset(log(n_months)), family = poisson(link = "log"), data = df_yearly)
fit_quasi   <- glm(count_any ~ year + offset(log(n_months)), family = quasipoisson(link = "log"), data = df_yearly)
fit_nb      <- MASS::glm.nb(count_any ~ year + offset(log(n_months)), data = df_yearly)

# SECTION 10: Relative Goodness-of-Fit
AIC(fit_poisson, fit_nb)
dispersiontest(fit_poisson)

# SECTION 11: Absolute Goodness-of-Fit
pearson_p <- sum(residuals(fit_poisson, type = "pearson")^2)
pchisq(pearson_p, df.residual(fit_poisson), lower.tail = FALSE)

set.seed(42)
sim_pois <- simulateResiduals(fit_poisson, n = 1000)
testUniformity(sim_pois)

# SECTION 12: Fitted Values Plot
df_yearly$fit_poisson <- fitted(fit_poisson)

df_plot <- df_yearly %>%
  dplyr::select(year, count_warm, count_cold, fit_poisson) %>%
  pivot_longer(cols = c(count_warm, count_cold), names_to = "type", values_to = "count") %>%
  mutate(type = ifelse(type == "count_warm", "Extreme Warm", "Extreme Cold"))

ggplot(df_plot, aes(x = year, y = count, fill = type)) +
  geom_col(position = "stack", width = 0.9) +
  geom_line(
    data = df_yearly,
    aes(x = year, y = fit_poisson),
    inherit.aes = FALSE,
    colour = "black", linewidth = 1, linetype = "dashed"
  ) +
  scale_fill_manual(
    values = c("Extreme Warm" = "#d73027", "Extreme Cold" = "#4575b4"),
    name = "Event Type"
  ) +
  scale_x_continuous(breaks = seq(1900, 2010, by = 25)) +
  scale_y_continuous(breaks = 0:4, expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Year",
    y = "Number of Extreme Months",
    title = "Annual Frequency of Extreme Temperature Months in the UK",
    subtitle = paste0(
      "Extreme warm: >", round(thresh_hot, 1), "°C; ",
      "Extreme cold: <", round(thresh_cold, 1), "°C",
      " (WMO 1961–1990 baseline)"
    )
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(size = 9, colour = "grey40")
  )


# SECTION 13: Final Summary Table
# (This produces the consolidated stats shown in your report's Table 1)
summary_stats <- data.frame(
  Model = c("Poisson", "NegBin"),
  AIC = c(AIC(fit_poisson), AIC(fit_nb)),
  BIC = c(BIC(fit_poisson), BIC(fit_nb))
)
print(summary_stats)

# SECTION 14: Interpretation of Best Model
coef_year <- coef(fit_poisson)["year"]
ci_year <- confint(fit_poisson)["year", ]
cat("Rate change per decade:", round((exp(coef_year * 10) - 1) * 100, 2), "%\n")

# For the Poisson model (the one chosen as best in your report)
summary(fit_poisson)

# For the Negative Binomial model
summary(fit_nb)

# For the Quasi-Poisson model
summary(fit_quasi)

