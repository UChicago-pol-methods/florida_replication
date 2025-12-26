#+ setup, include=FALSE
knitr::opts_chunk$set(dev = 'png', dev.args = list(type = "cairo"))

#'
library(ggplot2)
library(ggdist)
library(dplyr)
library(distributional)

load("../data/models.RData")

# List of models for each outcome
models <- list(
  `(1) Affect` = list(mod1_therm, mod2_therm, mod3_therm),
  `(2) Tolerance` = list(mod1_tlr, mod2_tlr, mod3_tlr),
  `(3) Policy` = list(mod1_law, mod2_law, mod3_law)
)

# Extract data from models
mod_data <- purrr::map_dfr(names(models), function(outcome) {
  purrr::imap_dfr(models[[outcome]], function(model, model_num) {
    # Extract treatment effect and confidence intervals
    coef <- model$coefficients["treat_indTRUE"]
    std_error <- model$std.error["treat_indTRUE"]
    conf_low <- model$conf.low["treat_indTRUE"]
    conf_high <- model$conf.high["treat_indTRUE"]
    
    # Create a data frame for this model
    tibble(
      outcome = outcome,
      Model = paste0("Model ", model_num),
      estimate = coef,
      std_error = std_error,
      conf_low = conf_low,
      conf_high = conf_high
    )
  })
})

# Update model labels for legend
mod_data <- mod_data |> 
  mutate(Model = factor(Model, 
                        levels = c("Model 1", "Model 2", "Model 3"),
                        labels = c("Difference-in-means", "Covariate adjusted", "Covariate adjusted\nand reweighted")))

# Create the plot
ggplot(mod_data, aes(ydist = distributional::dist_normal(estimate, std_error), 
                     fill = Model)) +
  ggdist::stat_gradientinterval(position = "dodge", .width = c(.9, .95)) +
  facet_wrap(~outcome, scales = "free_y") +
  scale_fill_manual(
    values = c("#D55E00", "#0072B2", "#CC79A7"),  # Colorblind-friendly palette
    name = "Model Type"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "lightgray") + 
  labs(
    x = "Model",
    y = "Treatment effect estimates"
  ) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), 
    legend.title = element_blank()
  )

# Save as high-resolution PNG for LaTeX
ggsave("../figures/Figure1_treatment_effects.png", width = 10, height = 5, dpi = 600)
# Save as TIFF for journal submission
ggsave("../figures/Figure1_treatment_effects.tiff", width = 10, height = 5, dpi = 600, compression = "lzw")

# by subgroups

load("../data/models_dems.RData")   # Democrats
load("../data/models_reps.RData")   # Republicans
load("../data/models_inds.RData")   # Independents

# Define groups and their models
models_list <- list(
  "Overall" = list(
    `(1) Affect` = list(mod1_therm, mod2_therm, mod3_therm),
    `(2) Tolerance` = list(mod1_tlr, mod2_tlr, mod3_tlr),
    `(3) Policy` = list(mod1_law, mod2_law, mod3_law)
  ),
  "Democrats" = list(
    `(1) Affect` = list(mod1_thermd, mod2_thermd, mod3_thermd),
    `(2) Tolerance` = list(mod1_tlrd, mod2_tlrd, mod3_tlrd),
    `(3) Policy` = list(mod1_lawd, mod2_lawd, mod3_lawd)
  ),
  "Republicans" = list(
    `(1) Affect` = list(mod1_thermr, mod2_thermr, mod3_thermr),
    `(2) Tolerance` = list(mod1_tlrr, mod2_tlrr, mod3_tlrr),
    `(3) Policy` = list(mod1_lawr, mod2_lawr, mod3_lawr)
  ),
  "Independents" = list(
    `(1) Affect` = list(mod1_thermi, mod2_thermi, mod3_thermi),
    `(2) Tolerance` = list(mod1_tlri, mod2_tlri, mod3_tlri),
    `(3) Policy` = list(mod1_lawi, mod2_lawi, mod3_lawi)
  )
)

# Extract data from models
mod_data <- purrr::imap_dfr(models_list, function(models, group) {
  purrr::map_dfr(names(models), function(outcome) {
    purrr::imap_dfr(models[[outcome]], function(model, model_num) {
      tibble(
        group = group,
        outcome = outcome,
        Model = paste0("Model ", model_num),
        estimate = model$coefficients["treat_indTRUE"],
        std_error = model$std.error["treat_indTRUE"],
        conf_low = model$conf.low["treat_indTRUE"],
        conf_high = model$conf.high["treat_indTRUE"]
      )
    })
  })
})

# Update factor levels for legend ordering
mod_data <- mod_data |> 
  mutate(
    Model = factor(Model, levels = c("Model 1", "Model 2", "Model 3"),
                   labels = c("Difference-in-means", "Covariate adjusted", "Covariate adjusted\nand reweighted")),
    group = factor(group, levels = c("Overall", "Democrats", "Republicans", "Independents"))
  )

# Define colors for subgroups & overall effects
color_palette <- c(
  "Overall" = "#000000",        # Black for overall effect
  "Democrats" = "#D55E00",      # Orange
  "Republicans" = "#0072B2",    # Blue
  "Independents" = "#CC79A7"    # Pink
)

# Compute y-axis positions for group labels
bracket_positions <- mod_data |> 
  group_by(outcome) |> 
  summarize(max_est = max(conf_high), min_est = min(conf_low)) |> 
  mutate(label_y = min_est - (max_est - min_est) * 0.1,
         label_ybar = min_est - (max_est - min_est) * 0.05)  # Adjust spacing

# Create the plot
ggplot(mod_data, aes(ydist = distributional::dist_normal(estimate, std_error), 
                     fill = Model, x = group)) +
  ggdist::stat_gradientinterval(position = position_dodge(width = 0.6),
                                point_size = 0.5, interval_size_range = c(0.25, 0.75), .width = c(.9, .95)) +
  facet_wrap(~outcome, scales = "free_y") +
  scale_fill_manual(
    values = c("#D55E00", "#0072B2", "#CC79A7"),  # Colorblind-friendly palette
    name = "Model Type"
  ) +
  geom_text(data = bracket_positions, aes(x = 1, y = label_y, label = "Overall"), 
            inherit.aes = FALSE, size = 2) +
  geom_text(data = bracket_positions, aes(x = 2, y = label_y, label = "Democrats"), 
            inherit.aes = FALSE, size = 2) +
  geom_text(data = bracket_positions, aes(x = 3, y = label_y, label = "Republicans"), 
            inherit.aes = FALSE, size = 2) +
  geom_text(data = bracket_positions, aes(x = 4, y = label_y, label = "Independents"),
            inherit.aes = FALSE, size = 2) +
  geom_segment(data = bracket_positions, 
               aes(x = 0.8, xend = 1.2, y = label_ybar, yend = label_ybar), 
               inherit.aes = FALSE, color = "darkgray") +
  geom_segment(data = bracket_positions, 
               aes(x = 1.8, xend = 2.2, y = label_ybar, yend = label_ybar), 
               inherit.aes = FALSE, color = "darkgray") +
  geom_segment(data = bracket_positions, 
               aes(x = 2.8, xend = 3.2, y = label_ybar, yend = label_ybar), 
               inherit.aes = FALSE, color = "darkgray") +
  geom_segment(data = bracket_positions, 
               aes(x = 3.8, xend = 4.2, y = label_ybar, yend = label_ybar), 
               inherit.aes = FALSE, color = "darkgray") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "lightgray") + 
  labs(
    x = "Model",
    y = "Treatment effect estimates"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank()
  )



# Save as high-resolution PNG for LaTeX
ggsave("../figures/FigureSI2_treatment_effects_by_subgroup.png", width = 10, height = 5, dpi = 600)
# Save as TIFF for journal submission
ggsave("../figures/FigureSI2_treatment_effects_by_subgroup.tiff", width = 10, height = 5, dpi = 600, compression = "lzw")



