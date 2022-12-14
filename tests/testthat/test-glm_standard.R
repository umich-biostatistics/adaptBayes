test_that("glm_standard_binomial runs, q = 0", {
  expect_type(glm_standard(y = 0,
                           x_standardized = data.frame(x1 = 0, x2 = 1),
                           family = "binomial",
                           p = 2,
                           q = 0,
                           beta_orig_scale = 1,
                           beta_aug_scale = 1,
                           mc_warmup = 50,
                           mc_iter_after_warmup = 50,
                           only_prior = TRUE), "list")
})

test_that("glm_standard_gaussian runs, q = 0", {
  expect_type(glm_standard(y = 0,
                           x_standardized = data.frame(x1 = 0, x2 = 1),
                           family = "gaussian",
                           p = 2,
                           q = 0,
                           beta_orig_scale = 1,
                           beta_aug_scale = 1,
                           mc_warmup = 50,
                           mc_iter_after_warmup = 50,
                           only_prior = TRUE), "list")
})

test_that("glm_standard_binomial runs, q = 1", {
  expect_type(glm_standard(y = 0,
                           x_standardized = data.frame(x1 = 0, x2 = 1),
                           family = "binomial",
                           p = 1,
                           q = 1,
                           beta_orig_scale = 1,
                           beta_aug_scale = 1,
                           mc_warmup = 50,
                           mc_iter_after_warmup = 50,
                           only_prior = TRUE), "list")
})

test_that("glm_standard_gaussian runs, q = 1", {
  expect_type(glm_standard(y = 0,
                           x_standardized = data.frame(x1 = 0, x2 = 1),
                           family = "gaussian",
                           p = 1,
                           q = 1,
                           beta_orig_scale = 1,
                           beta_aug_scale = 1,
                           mc_warmup = 50,
                           mc_iter_after_warmup = 50,
                           only_prior = TRUE), "list")
})
