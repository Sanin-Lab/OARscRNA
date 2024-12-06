test_that("oar works with defaults - starting from Seurat", {
  skip_on_cran()
  
  oar_input <- readRDS(test_path("fixtures", "oar_data.rds"))
  oar_output <- readRDS(test_path("fixtures", "oar_out.rds"))
  check_oar <- oar(oar_input, seurat_v5 = T, cores = parallelly::availableCores()-2)
  check_oar <- check_oar@meta.data
  
  expect_snapshot(
    waldo::compare(oar_output, check_oar)
  )
})

test_that("oar works with no mismatch - starting from Seurat", {
  skip_on_cran()
  
  oar_input <- readRDS(test_path("fixtures", "oar_data.rds"))
  oar_output <- readRDS(test_path("fixtures", "oar_out_nm.rds"))
  check_oar <- oar(oar_input, mismatch = F, seurat_v5 = T, cores = parallelly::availableCores()-2)
  check_oar <- check_oar@meta.data
  
  expect_snapshot(
    waldo::compare(oar_output, check_oar)
  )
})

test_that("oar works with defaults - starting from Matrix", {
  skip_on_cran()
  
  oar_input <- readRDS(test_path("fixtures", "oar_data_m.rds"))
  oar_output <- readRDS(test_path("fixtures", "oar_out_m.rds"))
  check_oar <- oar(oar_input, seurat_v5 = F, cores = parallelly::availableCores()-2)
  
  expect_snapshot(
    waldo::compare(oar_output, check_oar)
  )
})

test_that("oar works with no mismatch - starting from Matrix", {
  skip_on_cran()
  
  oar_input <- readRDS(test_path("fixtures", "oar_data_m.rds"))
  oar_output <- readRDS(test_path("fixtures", "oar_out_m_nm.rds"))
  check_oar <- oar(oar_input, mismatch = F, seurat_v5 = F, cores = parallelly::availableCores()-2)
  
  expect_snapshot(
    waldo::compare(oar_output, check_oar)
  )
})
