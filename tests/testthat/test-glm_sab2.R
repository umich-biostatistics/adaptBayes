test_that("glm_sab2_binomial runs", {
  expect_type(glm_sab2(y = 0,
                       x_standardized = data.frame(x1 = 0, x2 = 1),
                       family = "binomial",
                       alpha_prior_mean = 0,
                       alpha_prior_cov = matrix(1),
                       aug_projection = matrix(1),
                       beta_orig_scale = 1,
                       beta_aug_scale = 1,
                       only_prior = TRUE), "list")
})

test_that("glm_sab2_gaussian runs", {
  expect_type(glm_sab2(y = 0,
                       x_standardized = data.frame(x1 = 0, x2 = 1),
                       family = "gaussian",
                       alpha_prior_mean = 0,
                       alpha_prior_cov = matrix(1),
                       aug_projection = matrix(1),
                       beta_orig_scale = 1,
                       beta_aug_scale = 1,
                       only_prior = TRUE), "list")
})

test_that("glm_sab2_simple_binomial runs", {
  expect_type(glm_sab2(y = 0,
                       x_standardized = data.frame(x1 = 0, x2 = 1),
                       family = "binomial",
                       alpha_prior_mean = 0,
                       alpha_prior_cov = matrix(1),
                       aug_projection = matrix(1),
                       beta_orig_scale = 1,
                       beta_aug_scale = 1,
                       phi_mean = 1,
                       phi_sd = 0,
                       psi_mean = 0,
                       psi_sd = 0,
                       only_prior = TRUE), "list")
})

test_that("glm_sab2_simple_gaussian runs", {
  expect_type(glm_sab2(y = 0,
                       x_standardized = data.frame(x1 = 0, x2 = 1),
                       family = "gaussian",
                       alpha_prior_mean = 0,
                       alpha_prior_cov = matrix(1),
                       aug_projection = matrix(1),
                       beta_orig_scale = 1,
                       beta_aug_scale = 1,
                       phi_mean = 1,
                       phi_sd = 0,
                       psi_mean = 0,
                       psi_sd = 0,
                       only_prior = TRUE), "list")
})
