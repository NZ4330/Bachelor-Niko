# package
library("tidyverse")
library("fredr")
library("dplyr")
library("readxl")
library("readr")
library("lubridate")
library("zoo")
library("tseries")
library(nlsur) #devtools::install_github("JanMarvin/nlsur")
library(systemfit)
library(dplyr)
fredr_set_key("ea90d424181505bd4afbd1f54bece638")




# Data
Shadow_rate <- read_csv("/Users/nikolaizwisler/Desktop/Bachelor/Rstudio/Bachelor so far/Formatted_Shadow_Rate_Data.csv") %>%
  rename(Rate = Value) %>%  # Brug de korrekte kolonnenavne
  mutate(Rate = as.numeric(gsub(",", ".", Rate)))

Crude_Oil <- fredr(series_id = "DCOILBRENTEU", 
                   observation_start = as.Date("2004-09-01"), 
                   observation_end = as.Date("2022-08-01")) %>% 
  select(date, Price = value) %>% 
  mutate(date = as.Date(date)) %>%  # Konverter til datoformat
  mutate(Month = format(date, "%Y-%m")) %>%  # Lav måneds-søjle i formatet "ÅÅÅÅ-MM"
  group_by(Month) %>%
  summarise(Price = last(Price)) # Beholder sidste værdi pr. måned
Crude_Oil$Month <- as.Date(paste0(Crude_Oil$Month, "-01"))

#Transformering
Crude_Oil$Diff_Price <- c(NA, diff(Crude_Oil$Price))  # Tilføj en NA først

Shadow_rate$Diff_rate <- c(NA, diff(Shadow_rate$Rate))

# plot
plot(Shadow_rate$Date, Shadow_rate$Diff_rate,
     xlab = "Date", ylab = "Shadow_rate$diff_rate",
     main = "Shadow rate Diff")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted", lwd = par("lwd"))

plot(Shadow_rate$Date, Shadow_rate$Rate, xlab = "Date", ylab = "Rate",
     main = "Shadow rate" )
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted", lwd = par("lwd"))

plot(Crude_Oil$Month, Crude_Oil$Diff_Price,
     xlab = "Date", ylab = "Crude_Oil$Diff_Price",
     main = "Crude oile Diff")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted", lwd = par("lwd"))

plot(Crude_Oil$Month, Crude_Oil$Price,
     xlab = "Date", ylab = "Crude_Oil$Price",
     main = "Crude oile")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted", lwd = par("lwd"))


Shadow_rate$regime <- ifelse(
  Shadow_rate$Date >= as.Date("2012-01-01") & Shadow_rate$Date <= as.Date("2019-12-31"),
  "Regime 1",
  ifelse(Shadow_rate$Date >= as.Date("2020-01-01"), "Regime 2", NA)
)

ggplot(Shadow_rate, aes(x = Date, y = Diff_rate, color = regime)) +
  geom_point(na.rm = TRUE) +
  
  # Vertical lines marking start/end of Regime 1, and start of Regime 2
  geom_vline(xintercept = as.Date("2012-01-01"), linetype = "dashed", color = "blue", size = 1) +
  geom_vline(xintercept = as.Date("2019-12-31"), linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = as.Date("2020-01-01"), linetype = "dashed", color = "red", size = 1) +
  
  # Labels
  annotate("text", x = as.Date("2012-01-01"), y = max(Shadow_rate$Diff_rate, na.rm = TRUE),
           label = "Start Regime 1", angle = 90, vjust = -0.5, hjust = 0, size = 3) +
  annotate("text", x = as.Date("2019-12-31"), y = max(Shadow_rate$Diff_rate, na.rm = TRUE),
           label = "End Regime 1", angle = 90, vjust = -0.5, hjust = 0, size = 3) +
  annotate("text", x = as.Date("2020-01-01"), y = max(Shadow_rate$Diff_rate, na.rm = TRUE),
           label = "Start Regime 2", angle = 90, vjust = -0.5, hjust = 0, size = 3) +
  
  labs(
    title = "Scatterplot with Regime Differentiation",
    x = "Date",
    y = "Change in Shadow Rate",
    color = "Regime"
  ) +
  theme_minimal()

# Analyse samt Bootstrap fra Tamás med lidt variation

# Regime 1
Regime_1_Rate <- subset(Shadow_rate, Date >= as.Date("2012-01-01") & Date <= as.Date("2019-12-31"))
Regime_1_Olie <- subset(Crude_Oil, Month >= as.Date("2012-01-01") & Month <= as.Date("2019-12-31"))

# Regime 2
Regime_2_Rate <- subset(Shadow_rate, Date >= as.Date("2020-01-01"))
Regime_2_Olie <- subset(Crude_Oil, Month >= as.Date("2020-01-01") )

# Regime 1
Regime_1 <- merge(Regime_1_Rate[, c("Date", "Diff_rate")],
                  Regime_1_Olie[, c("Month", "Diff_Price")],
                  by.x = "Date", by.y = "Month")

# Regime 2
Regime_2 <- merge(Regime_2_Rate[, c("Date", "Diff_rate")],
                  Regime_2_Olie[, c("Month", "Diff_Price")],
                  by.x = "Date", by.y = "Month")

# Covariance
cov_regime1 <- cov(na.omit(Regime_1[, c("Diff_Price", "Diff_rate")]))
cov_regime2 <- cov(na.omit(Regime_2[, c("Diff_Price", "Diff_rate")]))

# Brug dine tal her:
var_p1 <- cov_regime1[1,1]
var_r1 <- cov_regime1[2,2]
cov_pr1 <- cov_regime1[1,2]

var_p2 <- cov_regime2[1,1]
var_r2 <- cov_regime2[2,2]
cov_pr2 <- cov_regime2[1,2]

# Rigobon & Sack analytisk løsning
a <- cov_pr1 - cov_pr2
b <- var_p1 - var_p2 - var_r1 + var_r2
c <- -(cov_pr1 - cov_pr2)

alpha_roots <- polyroot(c(a, b, c))
alpha <- Re(alpha_roots[1])  # Brug realdelen

beta <- (cov_pr1 - alpha * var_r1) / var_p1

cat("Alpha =", alpha, "\n")
cat("Beta  =", beta, "\n")

### bootstap generation 

bootIter <- 1000
bootRes <- data.frame(alpha = rep(NA, bootIter), beta = rep(NA, bootIter))

set.seed(123)  # for reproducibility

for (i in 1:bootIter) {
  
  # Træk tilfældige observationer (med replacement) fra hvert regime
  boot1 <- Regime_1[sample(1:nrow(Regime_1), replace = TRUE), c("Diff_Price", "Diff_rate")]
  boot2 <- Regime_2[sample(1:nrow(Regime_2), replace = TRUE), c("Diff_Price", "Diff_rate")]
  
  # Beregn kovariansmatricer
  cov1 <- cov(boot1, use = "complete.obs")
  cov2 <- cov(boot2, use = "complete.obs")
  
  # Hent varians og kovarians
  var_p1 <- cov1[1,1]
  var_r1 <- cov1[2,2]
  cov_pr1 <- cov1[1,2]
  
  var_p2 <- cov2[1,1]
  var_r2 <- cov2[2,2]
  cov_pr2 <- cov2[1,2]
  
  # Analytisk løsning: kvadratisk ligning i alpha
  a <- cov_pr1 - cov_pr2
  b <- var_p1 - var_p2 - var_r1 + var_r2
  c <- -(cov_pr1 - cov_pr2)
  
  alpha_roots <- polyroot(c(a, b, c))
  alpha <- Re(alpha_roots[1])  # vælg første løsning
  
  if (is.finite(alpha)) {
    beta <- (cov_pr1 - alpha * var_r1) / var_p1
    if (is.finite(beta)) {
      bootRes$alpha[i] <- alpha
      bootRes$beta[i] <- beta
    }
  }
}

summary(bootRes)

cat("Alpha (mean ± CI):", mean(bootRes$alpha, na.rm = TRUE), 
    quantile(bootRes$alpha, probs = c(0.05, 0.95), na.rm = TRUE), "\n")

cat("Beta (mean ± CI):", mean(bootRes$beta, na.rm = TRUE), 
    quantile(bootRes$beta, probs = c(0.05, 0.95), na.rm = TRUE), "\n")


hist(bootRes$alpha, breaks = 30, main = "Bootstrap Distribution of Alpha", xlab = "Alpha")
hist(bootRes$beta, breaks = 30, main = "Bootstrap Distribution of Beta", xlab = "Beta")
