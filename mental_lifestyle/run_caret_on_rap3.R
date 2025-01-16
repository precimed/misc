install.packages(c('caret', 'doParallel', 'glmnet', 'recipes', 'dplyr', 'QSARdata',
                   'e1071','kernlab','rsample','randomForest','ggplot2'))
my_libraries = c('caret', 'doParallel', 'glmnet', 'recipes', 'dplyr', 'QSARdata', 
                 'e1071','data.table','kernlab','rsample','randomForest','ggplot2')
lapply(my_libraries, library, character.only = TRUE)

#system(paste('dx download /zillur/bip_bmi_smd/pheno/*fin2.csv'))
args = commandArgs(trailingOnly=TRUE)
pheno_file = args[1]
df1 = fread(pheno_file, data.table = FALSE)
cl <- makePSOCKcluster(48)
registerDoParallel(cl)

get_predictions = function(pheno, model_method = 'glmnet') {
  predictors1 = df1[, c(3:2912, 2915, 2916, which(colnames(df1) == pheno))]
  rownames(predictors1) = df1$FID
  
  predictors1[[pheno]] = as.numeric(as.character(predictors1[[pheno]]))
  predictors1$Sex = as.factor(predictors1$Sex)
  
  recipe_obj = recipe(as.formula(paste(pheno, "~ .")), data = predictors1) %>%
    step_center(all_numeric_predictors()) %>%
    step_scale(all_numeric_predictors()) %>%
    step_nzv(all_predictors()) %>%
    step_naomit(all_predictors(), all_outcomes())
  
  set.seed(1234)
  train_indices = createDataPartition(predictors1[[pheno]], p = 0.8, list = FALSE, times = 1)
  train_data = predictors1[train_indices, ]
  test_data = predictors1[-train_indices, ]
  
  prep_obj = prep(recipe_obj, training = train_data)
  train_transformed = juice(prep_obj)
  test_transformed = bake(prep_obj, new_data = test_data)
  
  fit_control = trainControl(method = "repeatedcv", number = 10, repeats = 100,
                             savePredictions = TRUE, verboseIter = FALSE)
  
  fit_model = train(as.formula(paste(pheno, "~ .")), data = train_transformed,
                    method = model_method, trControl = fit_control)
  
  predictions = predict(fit_model, newdata = test_transformed)
  actual = test_transformed[[pheno]]
  
  post_resam = postResample(pred = predictions, obs = actual)
  
  coefficients = NULL
  if (model_method == "glmnet" && !is.null(fit_model$finalModel)) {
    best_lambda = fit_model$bestTune$lambda
    coefficients = as.data.frame(as.matrix(coef(fit_model$finalModel, s = best_lambda)))
    predictor_names = rownames(coefficients)
    coefficients = data.frame(Predictor = predictor_names, Coefficient = coefficients[, 1])
  }
  
  summary = data.frame(
    Model = model_method,
    R_squared = post_resam["Rsquared"],
    RMSE = post_resam["RMSE"],
    MAE = post_resam["MAE"]
  )
  
  fwrite(summary, file = paste0(pheno, '_ ', model_method, '_summary_fin2.csv'), sep = '\t', quote = FALSE, na = 'NA')
  if (!is.null(coefficients)) {
    fwrite(coefficients, file = paste0(pheno, '_', model_method, '_coefficients_fin2.csv'), sep = '\t', quote = FALSE, na = 'NA')
  }
  
  # Residuals
  residuals = actual - predictions
  residuals_df = data.frame(Actual = actual, Predicted = predictions, Residuals = residuals)
  
  # Save Residuals to CSV
  fwrite(residuals_df, file = paste0(pheno, '_residuals_fin2.csv'), sep = '\t', quote = FALSE, na = 'NA')
  
  # Diagnostic plots
  # 1. Residuals vs Fitted Plot
  png(filename = paste0(pheno, '_residuals_vs_fitted_fin2.png'))
  ggplot(residuals_df, aes(x = Predicted, y = Residuals)) +
    geom_point() + geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
    ggtitle("Residuals vs Fitted") + ylab("Residuals") + xlab("Fitted Values")
  dev.off()
  
  # 2. Predicted vs Actual Plot
  png(filename = paste0(pheno, '_predicted_vs_actual_fin2.png'))
  ggplot(residuals_df, aes(x = Actual, y = Predicted)) +
    geom_point() + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +
    ggtitle("Predicted vs Actual") + ylab("Predicted Values") + xlab("Actual Values")
  dev.off()
  
  return(summary)
}

# Example usage
#md1 = get_predictions('MD', model_method = 'rf')
#system(paste('dx upload *.csv --path=/zillur/bip_bmi_smd/caret/'))
#obesity1 = get_predictions('Obesity', model_method = 'rf')
#system(paste('dx upload *.csv --path=/zillur/bip_bmi_smd/caret/'))
#physact1 = get_predictions('PhysAct', model_method = 'rf')
#system(paste('dx upload *.csv --path=/zillur/bip_bmi_smd/caret/'))

md2 = get_predictions('MD', model_method = 'glmnet')
system(paste('dx upload MD*.csv --path=/zillur/bip_bmi_smd/caret/'))
obesity2 = get_predictions('Obesity', model_method = 'glmnet')
system(paste('dx upload Obesity*.csv --path=/zillur/bip_bmi_smd/caret/'))
physact2 = get_predictions('PhysAct', model_method = 'glmnet')
system(paste('dx upload PhysAct*.csv --path=/zillur/bip_bmi_smd/caret/'))
system(paste('dx upload *.png --path=/zillur/bip_bmi_smd/caret/'))

stopCluster(cl)