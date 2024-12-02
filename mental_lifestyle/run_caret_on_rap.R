install.packages('caret')
install.packages('doParallel')
install.packages("glmnet")
library(glmnet)
library(doParallel)
library(data.table)
library(caret)

cl <- makePSOCKcluster(96)
registerDoParallel(cl)

system(paste('dx download /zillur/bip_bmi_smd/pheno/*fin2.csv'))

df1 = fread('proteomic_pheno_caucasian_imputed_filtered_fin2.csv',data.table = F)


# Step 1: Make Predictions
get_predictions = function(pheno){
    predictors1 = df1[,3:2912]
    rownames(predictors1) = df1$FID
    
    #Near zero variance
    nzv <- nearZeroVar(predictors1, saveMetrics= TRUE)
    
    #Correlated variables
    descrCor <-  cor(predictors1)
    highCorr <- sum(abs(descrCor[upper.tri(descrCor)]) > .999)
    highlyCorDescr <- findCorrelation(descrCor, cutoff = .75)
    predictors2 <- predictors1[,-highlyCorDescr]
    
    #Linear dependency 
    comboInfo <- findLinearCombos(predictors2)
    
    #Near zero variance
    nzv <- nearZeroVar(predictors1, saveMetrics= TRUE)
    
    # Create a partition
    set.seed(123)  # Setting the seed for reproducibility
    #trainIndex <- createDataPartition(predictors2, p = 0.6, list = FALSE)
    # Split the data into training and testing sets
    ##Question; split only predictors not the outcome?
    ##Question; also transform the outcomes?
    # Calculate the number of rows to sample
    n <- nrow(predictors2)
    half_n <- floor(n / 2)  # You can use floor() if you prefer rounding down
    # Sample half the rows
    random_rows <- sample(seq_len(n), size = half_n)
    # Subset the data frame
    #half_df <- df[random_rows, ]

    
    trainData <- predictors2[random_rows, ]
    testData <- predictors2[-random_rows, ]
    preProcValues <- preProcess(trainData, method = c("center", "scale"))
    trainTransformed <- predict(preProcValues, trainData)
    testTransformed <- predict(preProcValues, testData)

    trainTransformed$MD=df1$MD[match(rownames(trainTransformed),df1$FID)]
    testTransformed$MD=df1$MD[match(rownames(testTransformed),df1$FID)]
    trainTransformed$Obesity=df1$Obesity[match(rownames(trainTransformed),df1$FID)]
    testTransformed$Obesity=df1$Obesity[match(rownames(testTransformed),df1$FID)]
    trainTransformed$Age=df1$Age[match(rownames(trainTransformed),df1$FID)]
    testTransformed$Age=df1$Age[match(rownames(testTransformed),df1$FID)]
    trainTransformed$Sex=df1$Sex[match(rownames(trainTransformed),df1$FID)]
    testTransformed$Sex=df1$Sex[match(rownames(testTransformed),df1$FID)]
    trainTransformed$PhysAct=df1$PhysAct[match(rownames(trainTransformed),df1$FID)]
    testTransformed$PhysAct=df1$PhysAct[match(rownames(testTransformed),df1$FID)]
    
    
    #Corss validation
    fitControl <- trainControl(## 10-fold CV
                           method = "repeatedcv",
                           number = 10,
                           ## repeated ten times
                           repeats = 10)
    
     glmnetFit1 <- train(as.formula(paste(pheno, "~ .")), data = trainTransformed, 
                      method = "glmnet", 
                      trControl = fitControl,
                      verbose = TRUE)
  
    predictions <- predict(glmnetFit1, newdata = testTransformed)
    actual <- testTransformed[[pheno]]
    rsquared <- caret::R2(predictions, actual)
    rmse <- caret::RMSE(predictions, actual)
    mae <- caret::MAE(predictions, actual)
    # Print performance metrics
    print(paste("R-squared:", rsquared))
    print(paste("RMSE:", rmse))
    print(paste("MAE:", mae))
    
    # Ensure valid lambda value for coefficients extraction
    best_lambda <- glmnetFit1$finalModel$lambda.min
    if (is.null(best_lambda) || is.na(best_lambda)) {
        best_lambda <- glmnetFit1$bestTune$lambda
    }
    
    # Step 3: Inspect Model Coefficients
    final_model <- glmnetFit1$finalModel
    coeff_matrix <- coef(final_model, s = best_lambda)
    coefficients <- as.data.frame(as.matrix(coeff_matrix))
    
    # The names of the predictors
    predictor_names <- rownames(coefficients)
    
    # Create the dataframe
    results <- data.frame(
        Predictor = predictor_names,
        Coefficient = coefficients,
        R_squared = rsquared)
    fwrite(results,file=paste0(pheno,'_glmfit_result_summary1.csv'),sep='\t',quote = F,na = 'NA')
    return(results)
    
}

md1 = get_predictions('MD')
system(paste('dx upload *summary1.csv --path=/zillur/bip_bmi_smd/caret/'))
obesity1 = get_predictions('Obesity')
system(paste('dx upload *summary1.csv --path=/zillur/bip_bmi_smd/caret/'))
physact1 = get_predictions('PhysAct')
system(paste('dx upload *summary1.csv --path=/zillur/bip_bmi_smd/caret/'))


