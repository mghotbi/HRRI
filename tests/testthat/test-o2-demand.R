test_that("worked stoichiometric demand is exact", {
  x <- data.frame(Fe2=50, FeS=5, Mn2=2, NH4=2, acetate=2, CH4=.5, O2=1.3)
  z <- rri_o2_demand(x, fe2_col="Fe2", fes_col="FeS", mn2_col="Mn2",
    nh4_col="NH4", acetate_col="acetate", ch4_col="CH4", o2_supply_col="O2")
  expect_equal(z$o2_demand, 33.75)
  expect_equal(z$o2_deficit_ratio, 33.75/1.3)
  expect_equal(z$n_species_observed, 6L)
})
test_that("units, zero stock and missing values are explicit", {
  a <- rri_o2_demand(data.frame(CH4=1000), ch4_col="CH4", ch4_unit="umol_kg")
  b <- rri_o2_demand(data.frame(CH4=1), ch4_col="CH4")
  expect_equal(a$o2_demand, b$o2_demand)
  z <- rri_o2_demand(data.frame(Fe2=c(0,NA), O2=c(0,1)),
                     fe2_col="Fe2", o2_supply_col="O2")
  expect_true(is.na(z$o2_deficit_ratio[1]))
  expect_true(is.na(z$o2_demand[2]))
})
test_that("reservoir columns cannot be counted twice", {
  expect_error(rri_o2_demand(data.frame(x=1), fe2_col="x", mn2_col="x"),
               "counted twice")
  expect_error(rri_o2_demand(data.frame(x=factor("1")), fe2_col="x"), "numeric")
})
