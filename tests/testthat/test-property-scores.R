test_that("capacity and root composites enforce numeric observations", {
  c <- rri_capacity_index(data.frame(EAC=c(2,4),EDC=c(3,1)))
  expect_length(c$capacity_score,2)
  expect_true(all(c$coverage==1))
  expect_error(rri_capacity_index(data.frame(EAC=factor(c(2,4))),edc_col=NULL),"numeric")
  r <- rri_root_physio(data.frame(ROL=c(.2,.4)), rol_col="ROL",
                       directions=c(ROL=1))
  expect_length(r$root_physio_score,2)
})

test_that("root trait direction is explicit and reversible", {
  x <- data.frame(ROL = c(1, 2, 3))
  hi <- rri_root_physio(x, rol_col = "ROL", directions = c(ROL = 1),
                        scaling = "minmax")$root_physio_score
  lo <- rri_root_physio(x, rol_col = "ROL", directions = c(ROL = -1),
                        scaling = "minmax")$root_physio_score
  expect_equal(lo, 1 - hi)
})
test_that("single-row memory score preserves its dimension", {
  x <- data.frame(H_hysteresis=.2,H_axis="forcing",incomplete_return_frac=-.3)
  z <- rri_memory_index(x)
  expect_equal(nrow(z),1L)
  expect_true(is.finite(z$memory_index))
})
test_that("property summary keeps unavailable diagnostics missing", {
  rs <- data.frame(Physio=c(.2,.4,.6,.8),Soil=c(.3,.4,.7,.7),
                   Micro=c(.2,.5,.5,.9),RRI=c(.23,.43,.61,.8))
  res <- list(row_scores=rs,meta=list(graph=NULL))
  z <- rri_property_scores(res,soil_df=data.frame(EAC=1:4,EDC=4:1))
  expect_true(is.finite(z$property_scores["Capacity"]))
  expect_true(is.finite(z$property_scores["Connectivity"]))
  expect_true(is.na(z$property_scores["Kinetics"]))
  expect_true(is.na(z$property_scores["Memory"]))
})
