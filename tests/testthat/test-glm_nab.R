test_that("glm_nab_binomial runs", {
  expect_type(glm_nab(y = 0,
                      x_standardized = data.frame(x1 = 0, x2 = 1),
                      family = "binomial",
                      alpha_prior_mean = 0,
                      alpha_prior_cov = matrix(1),
                      beta_orig_scale = 1,
                      beta_aug_scale = 1,
                      beta_aug_scale_tilde = 1,
                      mc_warmup = 50,
                      mc_iter_after_warmup = 50,
                      only_prior = TRUE), "list")
})

test_that("glm_nab_gaussian runs", {
  expect_type(glm_nab(y = 0,
                      x_standardized = data.frame(x1 = 0, x2 = 1),
                      family = "gaussian",
                      alpha_prior_mean = 0,
                      alpha_prior_cov = matrix(1),
                      beta_orig_scale = 1,
                      beta_aug_scale = 1,
                      beta_aug_scale_tilde = 1,
                      mc_warmup = 50,
                      mc_iter_after_warmup = 50,
                      only_prior = TRUE), "list")
})

test_that("glm_nab_simple_binomial runs", {
  expect_type(glm_nab(y = 0,
                      x_standardized = data.frame(x1 = 0, x2 = 1),
                      family = "binomial",
                      alpha_prior_mean = 0,
                      alpha_prior_cov = matrix(1),
                      beta_orig_scale = 1,
                      beta_aug_scale = 1,
                      beta_aug_scale_tilde = 1,
                      phi_mean = 1,
                      phi_sd = 0,
                      eta_param = Inf,
                      mc_warmup = 50,
                      mc_iter_after_warmup = 50,
                      only_prior = TRUE), "list")
})

test_that("glm_nab_simple_gaussian runs", {
  expect_type(glm_nab(y = 0,
                      x_standardized = data.frame(x1 = 0, x2 = 1),
                      family = "gaussian",
                      alpha_prior_mean = 0,
                      alpha_prior_cov = matrix(1),
                      beta_orig_scale = 1,
                      beta_aug_scale = 1,
                      beta_aug_scale_tilde = 1,
                      phi_mean = 1,
                      phi_sd = 0,
                      eta_param = Inf,
                      mc_warmup = 50,
                      mc_iter_after_warmup = 50,
                      only_prior = TRUE), "list")
})
