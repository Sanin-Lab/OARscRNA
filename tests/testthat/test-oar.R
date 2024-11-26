test_that("oar works with defaults", {
  skip_on_cran()
  
  oar_input <- readRDS(test_path("fixtures", "pdcs.rds"))
  oar_output <- readRDS(test_path("fixtures", "oar_pdcs_out.rds"))
  check_oar <- oar(oar_input)
  check_oar_meta <- check_oar@meta.data
  
  expect_snapshot(
    waldo::compare(oar_output, check_oar_meta)
  )
})

test_that("oar works with no mismatch", {
  skip_on_cran()
  
  oar_input <- readRDS(test_path("fixtures", "pdcs.rds"))
  oar_output <- readRDS(test_path("fixtures", "oar_pdcs_no_mismatch_output.rds"))
  check_oar_no_mis <- oar(oar_input, mismatch = F)
  check_oar_no_mis_meta <- check_oar_no_mis@meta.data
  
  expect_snapshot(
    waldo::compare(oar_output, check_oar_no_mis_meta)
  )
}) 

test_that("oar produces error when object is not seuratv5", {
  pdcs_v3 <- readRDS(test_path("fixtures", "pdcs_v3.rds"))
  
  expect_error(oar(pdcs_v3))
})