library(tidyverse)
library(patchwork)
library(metafor)
library(brms)
library(rstan)
library(tidybayes)
library(marginaleffects)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores() - 1)

meta_analysis_data <- read_csv("data/meta_analysis_data.csv") |>
  janitor::clean_names() |>
  mutate(condition_code = case_when(
    condition_code == "control" ~ "Control",
    condition_code == "ammonia" ~ "Ammonia",
    condition_code == "placebo" ~ "Placebo"
  )) |>
  mutate(condition_code = factor(condition_code, levels = c("Control", "Ammonia", "Placebo")))

# Check mean-sd scaling
meta_analysis_data |>
  ggplot(aes(x = log(post_m), y = log(post_sd))) +
  geom_point() +
  geom_smooth(method = "lm")

# Calculate standardised single group means
meta_analysis_data <-  escalc(
  data = meta_analysis_data,
    measure = "SMN",
    mi = post_m,
    sdi = post_sd,
    ni = n
  )

# Fit meta-analytic network model 

meta_model <- brm(yi | se(sqrt(vi)) ~ condition_code + 
                    (condition_code | study_no) + 
                    (condition_code | group_no) +
                    (1 | effect_no),
                  data = meta_analysis_data,
                  chains = 4,
                  cores = 4,
                  seed = 1988,
                  # warmup = 4000,
                  # iter = 40000,
                  control = list(adapt_delta = 0.99, max_treedepth = 11)
                  )

plot(meta_model)

pp_check(meta_model)

# plot meta-analysis

meta_pred <- predictions(
  meta_model,
  re_formula = NA,
  newdata = datagrid(condition_code = c("Control", "Ammonia", "Placebo"))
) |>
  posterior_draws()

meta_pred_labels <- meta_pred |>
  group_by(condition_code) |>
  mean_qi(draw)

study_pred <- predictions(
  meta_model,
  re_formula = NULL,
) |>
  posterior_draws()

study_labels <- meta_analysis_data |>
  select(1,2) |>
  group_by(study_no) |>
  slice(1)

study_pred <- left_join(study_pred, study_labels, by = "study_no") |>
  group_by(authors) |>
  mutate(mean_draw = mean(draw))

study_pred_labels <- study_pred |>
  group_by(authors, condition_code) |>
  mean_qi(draw)


meta_slopes <- avg_slopes(
  meta_model,
  re_formula = NA,
  variables = "condition_code"
  ) |>
  posterior_draws()

meta_slopes_labels <- meta_slopes |>
  group_by(contrast) |>
  mean_qi(draw)


meta_pred_plot <- ggplot(meta_pred, aes(x = draw, fill = condition_code)) +
  geom_vline(xintercept = 0, lty = "dashed") +
  stat_halfeye(slab_alpha = .5) +
  facet_grid(condition_code~.) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00", "#009E73")) +
  geom_text(
    data = mutate_if(meta_pred_labels,
                     is.numeric, round, 2),
    aes(
      label = glue::glue("{draw} [{.lower}, {.upper}]"),
      x = mean(draw), y = 0.25
    ),
    size = 3
  ) +
  labs(
    x = "Standardised Mean Score",
    fill = "Condition",
    title = "Global Grand Mean Estimates for Condition"
  ) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(title = element_text(size=8))



study_pred_plot <- ggplot(study_pred, aes(x = draw, 
                                                 y = reorder(authors, mean_draw), 
                                                 fill = condition_code)) +
  geom_vline(xintercept = 0, lty = "dashed") +
  # Add individual study data
  geom_point(
    data = meta_analysis_data,
    aes(x = yi, y = authors, color = condition_code),
    position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.1),
    alpha = 0.5
    ) +
  stat_halfeye(slab_alpha = .5, position = position_dodge(width = 0.5), size = 0.25) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00", "#009E73")) +
  scale_color_manual(values = c("#56B4E9", "#E69F00", "#009E73")) +
  geom_text(
    data = mutate_if(study_pred_labels,
                     is.numeric, round, 2),
    aes(
      label = glue::glue("{condition_code}: {draw} [{.lower}, {.upper}]"),
      x = 40, y =reorder(authors, draw), group = condition_code
    ),
    size = 2, position = position_dodge(width = 0.75),
    hjust = "inward"
  ) +
  labs(
    x = "Standardised Mean Score",
    title = "Conditional Estimates for Condition by Study"
  ) +
  guides(
    fill = "none",
    color = "none"
  ) +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(title = element_text(size=8))


contrast_plot <- ggplot(meta_slopes, aes(x = draw)) +
  geom_vline(xintercept = 0, lty = "dashed") +
  stat_halfeye(slab_alpha = .5, fill = "black") +
  facet_grid(contrast~.) +
  geom_text(
    data = mutate_if(meta_slopes_labels,
                     is.numeric, round, 2),
    aes(
      label = glue::glue("{draw} [{.lower}, {.upper}]"),
      x = -6, y = 0.1
    ),
    size = 3, hjust = "inward"
  ) +
  labs(
    x = "Standardised Mean Score",
    title = "Contrasts Between Conditions"
  ) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(title = element_text(size=8))

meta_plots <- (study_pred_plot + meta_pred_plot + contrast_plot) + 
  plot_annotation(title = "Meta-Analysis of Prior Studies Examining the Effects of Ammonia Inhalants",
                  caption = "Point estimates and 95% quantile intervals reported") +
  plot_layout(guides = "collect", axis_titles = "collect",
              widths = c(2,1,1))  &
  theme(axis.title.x = element_text(size=10),
        legend.position = "bottom")

meta_plots

ggsave("meta_plots.tiff", dpi = 300, w=12.5, h=5)


# Get posterior to use as prior
draws <- meta_model |>
  spread_draws(b_Intercept, b_condition_codeAmmonia, b_condition_codePlacebo)

control_prior <- MASS::fitdistr(draws$b_Intercept, "t")$estimate

ammonia_prior <- MASS::fitdistr(draws$b_condition_codeAmmonia, "t")$estimate

placebo_prior <- MASS::fitdistr(draws$b_condition_codePlacebo, "t")$estimate


empirical_data <- read_csv("data/empirical_data.csv") |>
  janitor::clean_names() |>
  mutate(condition = factor(condition, levels = c("Control", "Ammonia", "Placebo"))) |>
  mutate(peak_force_z = peak_force_n/sd(peak_force_n),
         rfd_z = rfd_200ms_n_s/sd(rfd_200ms_n_s))

prior <- c(
  set_prior(paste("student_t(",control_prior[3],",", control_prior[1],",", control_prior[2],")"),  class = "Intercept"),
  set_prior(paste("student_t(",ammonia_prior[3],",", ammonia_prior[1],",", ammonia_prior[2],")"),  class = "b", coef = "conditionAmmonia"),
  set_prior(paste("student_t(",placebo_prior[3],",", placebo_prior[1],",", placebo_prior[2],")"),  class = "b", coef = "conditionPlacebo")
)

prior_only_model <-
  brm(
    peak_force_z ~ condition + (1 | name),
    data = empirical_data,
    chains = 4,
    cores = 4,
    seed = 1988,
    prior = prior,
    sample_prior = "only",
  )


prior_pred <- predictions(
  prior_only_model,
  re_formula = NA,
  newdata = datagrid(condition = c("Control", "Ammonia", "Placebo"))
) |>
  posterior_draws()

prior_slopes <- avg_slopes(
  prior_only_model,
  re_formula = NA,
  variables = "condition"
) |>
  posterior_draws()

# Fit empirical models

peak_force_model <- brm(peak_force_z ~ condition + (1 | name),
                  data = empirical_data,
                  chains = 4,
                  cores = 4,
                  seed = 1988,
                  prior = prior
                  # warmup = 4000,
                  # iter = 40000,
                  # control = list(adapt_delta = 0.99, max_treedepth = 11)
)

plot(peak_force_model)

pp_check(peak_force_model)

rfd_model <- brm(rfd_z ~ condition + (1 | name),
                        data = empirical_data,
                        chains = 4,
                        cores = 4,
                        seed = 1988,
                        prior = prior
                        # warmup = 4000,
                        # iter = 40000,
                        # control = list(adapt_delta = 0.99, max_treedepth = 11)
)

plot(rfd_model)

pp_check(rfd_model)


# Plot empirical models 

peak_force_pred <- predictions(
  peak_force_model,
  re_formula = NA,
  newdata = datagrid(condition = c("Control", "Ammonia", "Placebo"))
) |>
  posterior_draws() 

peak_force_pred_labels <- peak_force_pred |>
  group_by(condition) |>
  mean_qi(draw)

peak_force_slopes <- avg_slopes(
  peak_force_model,
  re_formula = NA,
  variables = "condition"
) |>
  posterior_draws() 

peak_force_slopes_labels <- peak_force_slopes |>
  group_by(contrast) |>
  mean_qi(draw)

peak_force_pred_plot <- ggplot(peak_force_pred, aes(x = draw, y = condition, fill = condition)) +
  geom_vline(xintercept = 0, lty = "dashed") +
  geom_point(
    data = empirical_data,
    aes(x = peak_force_z, y = condition, color = condition),
    position = position_nudge(y=-0.05),
    alpha = 0.5
  ) +
  stat_slab(data = prior_pred,
            aes(x = draw, y = condition),
            slab_alpha = .5, fill = "grey") +
  stat_halfeye(slab_alpha = .5,position = position_dodge()) +
  scale_y_discrete(limits = rev) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00", "#009E73")) +
  scale_color_manual(values = c("#56B4E9", "#E69F00", "#009E73")) +
  geom_text(
    data = mutate_if(peak_force_pred_labels,
                     is.numeric, round, 2),
    aes(
      label = glue::glue("{draw} [{.lower}, {.upper}]"),
      x = mean(draw), y = condition
    ),
    size = 3, position = position_nudge(y = 0.25)
  ) +
  labs(
    x = "Peak Force (Standardised Score)",
    fill = "Condition",
    title = "Global Grand Mean Estimates for Condition"
  ) +
  guides(
    color = "none"
  ) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(title = element_text(size=8))

peak_force_contrast_plot <- ggplot(peak_force_slopes, aes(x = draw)) +
  stat_slab(data = prior_slopes,
            slab_alpha = .5, fill = "grey") +
  geom_vline(xintercept = 0, lty = "dashed") +
  stat_halfeye(slab_alpha = .5, fill = "black") +
  facet_grid(contrast~.) +
  geom_text(
    data = mutate_if(peak_force_slopes_labels,
                     is.numeric, round, 2),
    aes(
      label = glue::glue("{draw} [{.lower}, {.upper}]"),
      x = -3, y = 0.1
    ),
    size = 3, hjust = "inward"
  ) +
  labs(
    x = "Peak Force (Standardised Score)",
    title = "Contrasts Between Conditions"
  ) +
  scale_x_continuous(limits = c(-3,3)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(title = element_text(size=8))

empirical_plots <- (peak_force_pred_plot + peak_force_contrast_plot) + 
  plot_annotation(title = "Posterior Estimates for Current Study (Prior Distributions in Lighter Grey)",
                  caption = "Point estimates and 95% quantile intervals reported") +
  plot_layout(guides = "collect", axis_titles = "collect")  &
  theme(axis.title.x = element_text(size=10))

empirical_plots

ggsave("empirical_plots.tiff", dpi = 300, w=10, h=5)
