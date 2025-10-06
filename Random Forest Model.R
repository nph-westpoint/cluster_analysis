library(tidyverse)
library(janitor)
library(randomForest)
library(iml)
# Load the dataset
data <- read.csv("real_world_synthetic.csv")

# Select relevant columns
data_rf <- data %>% select(-age_wave1,-sex_wave1)
  
library(caret)
dummies <- dummyVars(~ cluster, data = data_rf)
dummy_df <- predict(dummies, newdata = data_rf)

# Remove rows with missing values
data_rf_dummy <- data_rf %>% 
  select(-cluster) %>% 
  bind_cols(dummy_df)
data_rf_dummy <- na.omit(data_rf_dummy)

str(data_rf$cluster)
set.seed(123)
rf_model <- randomForest(auc_pred ~ ., data = data_rf_dummy, importance = TRUE, mtry = 2)
print(rf_model)
rf_model$importance
varImpPlot(rf_model) 

oob_mse  <- tail(rf_model$mse, 1)
oob_rmse <- sqrt(oob_mse)
oob_r2   <- tail(rf_model$rsq, 1)  # %Var explained ~ R^2

cat("\nRandom Forest OOB performance:\n",
    "  RMSE  :", round(oob_rmse, 3), "\n",
    "  R^2   :", round(oob_r2, 3),   "\n")
# Create a Predictor object
X <- data_rf_dummy %>% select(-auc_pred)  # predictors
y <- data_rf_dummy$auc_pred
predictor <- Predictor$new(rf_model, data = X, y = y)

# Compute SHAP values
shap <- Shapley$new(predictor, x.interest = X[1, ])

# Plot SHAP values for the first observation
plot(shap)

# Change this butt ugly graph to something nice
library(ggplot2)
library(dplyr)

# Extract SHAP results as a data.frame
shap_df <- shap$results %>%
  as.data.frame() %>%
  arrange(desc(abs(phi)))   # order by absolute importance

# Clean up variable names if needed
shap_df$feature <- factor(shap_df$feature, levels = shap_df$feature)

# --- GLOBAL SHAP SUMMARY (mean |SHAP| across all rows) ---
install.packages("fastshap")   # uncomment if needed
library(fastshap)
library(dplyr)
library(ggplot2)
library(tidyr)

# X and y should already exist from your code:
# X <- data_rf_dummy %>% dplyr::select(-auc_pred)
# y <- data_rf_dummy$auc_pred

# Wrapper for predict()
pred_fun <- function(object, newdata) predict(object, newdata)

set.seed(123)
shap_mat <- fastshap::explain(
  object       = rf_model,
  X            = as.data.frame(X),
  pred_wrapper = pred_fun,
  nsim         = 100,      # increase for more stable estimates (e.g., 200–500)
  adjust       = TRUE
)
# shap_mat is an N x P matrix of SHAP values

# Summarize to global importance = mean absolute SHAP per feature
global_shap <- shap_mat |>
  as.data.frame() |>
  mutate(.row = row_number()) |>
  pivot_longer(cols = - .row, names_to = "feature", values_to = "shap") |>
  group_by(feature) |>
  summarise(mean_abs_shap = mean(abs(shap), na.rm = TRUE), .groups = "drop") |>
  arrange(desc(mean_abs_shap))

print(global_shap)

# Pretty bar chart
ggplot(global_shap, aes(x = reorder(feature, mean_abs_shap), y = mean_abs_shap)) +
  geom_col(fill = "#2E86AB") +
  coord_flip() +
  labs(
    title = "Global SHAP Importance (Random Forest)",
    subtitle = "Mean absolute SHAP value across all observations",
    x = "Feature",
    y = "Mean |SHAP|"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 11)
  )


# Get model prediction and average baseline
y_hat <- predict(rf_model, X[1, , drop = FALSE])
y_mean <- mean(y)

# Nice SHAP plot
ggplot(shap_df, aes(x = feature, y = phi, fill = phi > 0)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#2E86AB", "FALSE" = "#E94F37")) +
  labs(
    title = "SHAP Explanation for One Observation",
    subtitle = paste0("Actual prediction = ", round(y_hat, 2),
                      " | Average prediction = ", round(y_mean, 2)),
    x = "Feature",
    y = "SHAP value (impact on prediction)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 11)
  )
