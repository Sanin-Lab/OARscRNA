test_that("oar_base works when searching for tolerance", {
  oar_base_input <- readRDS(test_path("fixtures", "oar_data.rds"))
  mdp <- readRDS(test_path("fixtures", "mdp_1.rds"))
  oar_base_search_tol <- readRDS(test_path("fixtures", "oar_base_search_tol.rds"))
  check_oar_base_search <- oar_base(oar_base_input, mdp)
  
  expect_snapshot(
    waldo::compare(oar_base_search_tol, check_oar_base_search) 
  )
})

test_that("oar_base works, user input tolerance", {
  oar_base_input <- readRDS(test_path("fixtures", "oar_data.rds"))
  mdp <- readRDS(test_path("fixtures", "mdp_2.rds"))
  oar_base_give_tol <- readRDS(test_path("fixtures", "oar_base_give_tol.rds"))
  check_oar_base_give <- oar_base(oar_base_input, mdp)
  
  expect_snapshot(
    waldo::compare(oar_base_give_tol, check_oar_base_give) 
  )
})