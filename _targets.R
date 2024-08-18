
source("_targets_packages.R")

tar_option_set(
  controller = crew_controller_local(workers = 50)
)

# Load and preprocess data
tar_data_pbc <- tar_target(
  data_pbc,
  {
    pbc %>%
      drop_na(time, status) %>%
      mutate(
        status = case_when(
          status <= 1 ~ 0, # code transplant as censored
          status == 2 ~ 1  # convert status to a 0/1 column
        ),
        across(
          .cols = c(sex, ascites, hepato, spiders, edema, trt, stage),
          .fns = factor
        ),
        surv = Surv(time, status)
      ) %>%
      relocate(surv, .before = 1) %>%
      select(-id, -time, -status)
  }
)

# Setup resampling for model validation
tar_data_resamples <- tar_target(
  data_resamples,
  vfold_cv(data_pbc)
)

# Generate recipes based on a grid of parameters
tar_recipe <- tar_target(
  recipe,
  {
    recipes_init <- expand_grid(
      boxcox = c(TRUE, FALSE),
      pca = c(0, 3, 5),
      spline_b = c(0, 5)
    )

    recipes_init %>%
      mutate(
        recipe = pmap(
          .l = list(boxcox, pca, spline_b),
          .f = function(use_boxcox, use_pca, use_spline_b){
            rec_init <- recipe(surv ~ ., data = data_pbc) %>%
              step_impute_mean(all_numeric_predictors()) %>%
              step_impute_mode(all_nominal_predictors())

            if (use_boxcox) {
              rec_init <- rec_init %>%
                step_BoxCox(all_numeric_predictors())
            }

            if (use_pca > 0) {
              rec_init <- rec_init %>%
                step_pca(all_numeric_predictors(),
                         num_comp = use_pca)
            }

            if(use_spline_b > 0) {
              rec_init <- rec_init %>%
                step_spline_b(all_numeric_predictors(),
                              deg_free = use_spline_b)
            }

            rec_init %>%
              step_dummy(all_nominal_predictors())
          }
        )
      ) %>%
      mutate(
        boxcox = factor(boxcox, labels = c("box_no", "box_all")),
        pca = factor(pca, labels = c("pca_no", "pca_3", "pca_5")),
        spline_b = factor(spline_b, labels = c("spline_no", "spline_5"))
      ) %>%
      unite(col = 'name', boxcox, pca, spline_b) %>%
      mutate(name = paste0(name, "..")) %>%
      deframe()

  }
)

tar_model_spec <- tar_target(
  model_spec,
  {
    list(
      ph_glmnet = proportional_hazards(penalty = tune(),
                                       mixture = 1) %>%
        set_engine("glmnet") %>%
        set_mode("censored regression"),

      survreg = survival_reg(dist = tune()) %>%
        set_engine("survival") %>%
        set_mode("censored regression"),

      aorsf = rand_forest(min_n = tune()) %>%
        set_engine("aorsf") %>%
        set_mode("censored regression"),

      mboost = boost_tree() %>%
        set_engine("mboost") %>%
        set_mode("censored regression")
    )
  }
)

tar_wflow <- tar_target(
  wflow,
  {
    wf <- workflow_set(preproc = recipe, models = model_spec)

    survival_metrics <- metric_set(brier_survival,
                                   roc_auc_survival,
                                   brier_survival_integrated)

    evaluation_time_points <- seq(500, 3000, length.out = 100)

    workflow_map(
      object = wf,
      fn = 'tune_grid',
      resamples = data_resamples,
      metrics = survival_metrics,
      eval_time = evaluation_time_points,
      control = control_resamples(save_pred = TRUE)
    )

  },
  pattern = cross(recipe, model_spec)
)

tar_fig <- tar_target(
  fig,
  {

    data_gg <- wflow %>%
      rank_results("brier_survival_integrated", select_best = TRUE) %>%
      filter(.metric == "brier_survival_integrated") %>%
      separate(col = wflow_id,
               into = c(".preproc", ".model"),
               sep = '\\.\\._',
               remove = FALSE)

    data_gg <- data_gg %>%
      arrange(.preproc, desc(rank)) %>%
      mutate(text = NA_character_) %>%
      split(.$.preproc) %>%
      map2_dfr(names(.), ~add_row(.x, text = .y, .before = 1)) %>%
      mutate(x = rev(seq(n())))

    min_val <- min(data_gg$mean, na.rm = TRUE)

    fig <- ggplot(data_gg, aes(x = x,
                        y = mean,
                        ymin = mean - 1.96 * std_err,
                        ymax = mean + 1.96 * std_err,
                        fill = .preproc)) +
      geom_col(show.legend = FALSE) +
      coord_flip() +
      geom_text(aes(label = .model, y = mean, hjust = -0.1)) +
      geom_text(aes(label = text, y = 0, hjust = 0)) +
      geom_hline(yintercept = min_val, linetype = 2) +
      theme_minimal() +
      scale_fill_viridis_d() +
      theme(axis.text.y = element_blank(),
            panel.grid = element_blank(),
            axis.line.x = element_line()) +
      labs(y = "Integrated Brier score",
           x = "",
           fill = "Pre-processing")

    # ggsave('tmp.png', fig, width = 6, height = 12, dpi = 300)

  }
)

# Define the list of targets for the pipeline
list(
  tar_data_pbc,
  tar_data_resamples,
  tar_recipe,
  tar_model_spec,
  tar_wflow,
  tar_fig
)
