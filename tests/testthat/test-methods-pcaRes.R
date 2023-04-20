test_that("DModX is invariant to the number of selected observations", {
    data("mtcars")
    ind <- sample(1:nrow(mtcars), 5)
    test <- mtcars[ind,]
    train <- mtcars[-ind,]
    pcaModel <- pca(train, scale="uv", center=TRUE, nPcs=3)
    fiveObs <- DModX(pcaModel, newdata=test, type="normalized")
    threeObs <- DModX(pcaModel, newdata=test[1:3,], type="normalized")
    expect_equal(fiveObs[1], threeObs[1])
})
