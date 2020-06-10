test_that("blblm-lenght of coefficent matches", {
  a <- blblm(mpg ~ hp + cyl + cyl*hp , mtcars, m = 2 , B=1000)
  x <- model.matrix(mpg ~ hp + cyl + cyl * hp , mtcars)
  y <- model.response(model.frame(mpg ~ hp + cyl + cyl * hp , mtcars))
  w <- rmultinom(1, 10, rep(1, nrow(mtcars)))
  b <- lm.wfit(x,y,w)
  expect_equal(length(coef(a)), length(coef(b)))
})

test_that("coefficient name matches", {
  a <- blblm(mpg ~ hp + cyl + cyl*hp , mtcars, m = 2 , B=1000)
  x <- model.matrix(mpg ~ hp + cyl + cyl * hp , mtcars)
  y <- model.response(model.frame(mpg ~ hp + cyl + cyl * hp , mtcars))
  w <- rmultinom(1, 10, rep(1, nrow(mtcars)))
  b <- lm.wfit(x,y,w)
  expect_equivalent(names(coef(a)), names(coef(b)))
})