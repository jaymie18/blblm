test_that("multiplication works", {
  a <- lm(mpg ~ hp + cyl + cyl*hp , mtcars)
  b <- blblm(mpg ~ hp + cyl + cyl*hp , mtcars, m = 3, B = 500)
  expect_equivalent(formula(a), formula(b))
})
