test_that("glm_studt_binomial runs", {
  expect_type(glm_studt(y = 0,
                       x_standardized = data.frame(x1 = 0, x2 = 1),
                       family = "binomial",
                       beta_scale = 1,
                       mc_warmup = 50,
                       mc_iter_after_warmup = 50,
                       only_prior = TRUE), "list")
})

test_that("glm_studt_gaussian runs", {
  expect_type(glm_studt(y = 0,
                       x_standardized = data.frame(x1 = 0, x2 = 1),
                       family = "gaussian",
                       beta_scale = 1,
                       mc_warmup = 50,
                       mc_iter_after_warmup = 50,
                       only_prior = TRUE), "list")
})
