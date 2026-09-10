# Extracted from test-cpp-implem-helpers.R:10

# test -------------------------------------------------------------------------
expect_equal(object = ilrInvcpp(Pindep, norm = "orthogonal"), expected = ilrInv(Pindep, norm = "orthogonal"), tolerance = 1e-6)
