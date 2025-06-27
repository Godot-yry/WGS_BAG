rm(list=ls())
library(dplyr)
library(data.table)
library(readr)
library(caret)
library(parallel)
library(glmnet)
library(doParallel)
library(foreach)

#list organ specific protein
organ_specific_path<-"./select_protein4_impute_onlymean/"
organ_base_filename<-c("Adrenal","Artery","Brain","Esophagus","Heart","Immune","Intestine","Kidney","Liver","Lung","Muscle","Pancreas","Pituitary","Salivary","Skin","Stomach","Thyroid")
Healthy_group<-fread("./healthy_group_eid_onlyhes.csv")
save_path<-"./elg_results/"
dir.create(save_path)


for(organ in c(1:17)){
  print(organ_base_filename[organ])
  collect_eval<-list()
  
  current_protein1 <- fread(paste0(organ_specific_path,organ_base_filename[organ],"_select_protein_new.csv"))%>%
    subset(eid %in% Healthy_group$eid)%>%
    select(-"eid")%>%
    na.omit()%>%
    as.data.frame()
  
  current_protein_nocov<-current_protein%>%
    select(-"Sex",-"Age")
  
  protein_name<-colnames(current_protein_nocov)
  rm(current_protein_nocov)
  current_protein$Sex<-as.factor(current_protein$Sex)
  
  current_protein_predage<-fread(paste0(organ_specific_path,organ_base_filename[organ],"_select_protein_new.csv"))%>%
    subset(eid %in% Healthy_group$eid)%>%
    select("eid","Age")
  current_protein_predage$predAge<-NA
  current_protein_predage$gap_resid<-NA
  
  for(gentle in c(0,1)){
    current_protein_gentle<-subset(current_protein,Sex==gentle)
    
    current_protein_predage_gentle<-fread(paste0(organ_specific_path,organ_base_filename[organ],"_select_protein_new.csv"))%>%
      subset(eid %in% Healthy_group$eid)%>%
      subset(Sex==gentle)%>%
      select("eid","Age")
    current_protein_predage_gentle$predAge<-NA
    current_protein_predage_gentle$gapresid<-NA
    
    #########training model
    folds <- createFolds(current_protein_gentle$Age, k = 5)
    list_models<-list()
    list_parameter<-list()
    for (k in c(1:5)){
      print(k)
      test_indices <- folds[[k]]
      train_indices <- setdiff(1:nrow(current_protein_gentle), test_indices)
      
      current_protein_train_data <- current_protein_gentle[train_indices,]
      current_protein_test_data <- current_protein_gentle[test_indices,]
      
      X_train<-as.matrix(current_protein_train_data[,protein_name])
      Y_train<-as.matrix(current_protein_train_data[,"Age"])
      
      numCores <- detectCores() - 4
      cl <- makeCluster(numCores)
      registerDoParallel(cl)
      
      alpha_grid <- seq(0, 1, by = 0.1)  
      lambda_grid <- exp(seq(log(0.0001), log(10), length.out = 10))  
      
      tuneGrid <- expand.grid(alpha = alpha_grid, lambda = lambda_grid)
      
      results <- foreach(i = 1:length(alpha_grid), .packages = "glmnet") %dopar% {
        set.seed(123)  
        cv_model <- cv.glmnet(X_train, Y_train, alpha = alpha_grid[i], family = "gaussian",type.measure = 'mae', nfolds = 5)
        
        best_lambda <- cv_model$lambda.min
        min_cv_error <- min(cv_model$cvm)
        
        return(list(model = cv_model, alpha = alpha_grid[i], lambda = best_lambda, cv_error = min_cv_error))
      }
      
      best_index <- which.min(sapply(results, function(x) x$cv_error))
      best_model <- results[[best_index]]$model 
      best_params <- list(alpha = results[[best_index]]$alpha, lambda = results[[best_index]]$lambda)
      
      list_models[[k]]<-best_model
      list_parameter[[k]]<-best_params
      
      predictions_train <- predict(best_model, newx = as.matrix(X_train), s = best_params$lambda)
      results_train<-postResample(pred = predictions_train,obs = Y_train)
      
      X_test<-as.matrix(current_protein_test_data[,protein_name])
      Y_test<-as.matrix(current_protein_test_data[,"Age"])
      predictions_test <- predict(best_model, newx = as.matrix(X_test), s = best_params$lambda)
      
      results_test<-postResample(pred = predictions_test,obs = Y_test)
      
      r_predic_train <- cor(Y_train, predictions_train)
      collect_eval$fold[k]<-k
      collect_eval$r_pred_train[k]<-r_predic_train
      collect_eval$RMSE_train[k]<-results_train[1]
      collect_eval$Rsquared_train[k]<-results_train[2]
      collect_eval$MAE_train[k]<-results_train[3]
      
      r_predic_test <- cor(Y_test, predictions_test)
      collect_eval$r_pred_test[k]<-r_predic_test
      collect_eval$RMSE_test[k]<-results_test[1]
      collect_eval$Rsquared_test[k]<-results_test[2]
      collect_eval$MAE_test[k]<-results_test[3]
      
      current_protein_predage_gentle$predAge[test_indices]<-predictions_test
      
    }
    
    collect_eval<-as.data.frame(collect_eval)
    rm(best_model,best_index,best_params)
    
    ###age gap correction
    
    for(k in c(1:5)){
      test_indices <- folds[[k]]
      train_indices <- setdiff(1:nrow(current_protein_gentle), test_indices)
      
      ytrain <- current_protein_predage_gentle$Age[train_indices]        
      ytrain_predic <- current_protein_predage_gentle$predAge[train_indices]  
      gap <- ytrain_predic - ytrain        
      
      b <- glm(gap ~ ytrain, family = gaussian)  
      
      ytest <- current_protein_predage_gentle$Age[test_indices]     
      ytest_predic <- current_protein_predage_gentle$predAge[test_indices] 
      
      gap <- ytest_predic - ytest      
      
      yfit <- predict(b, newdata = data.frame(ytrain = ytest), type = "response")  
      current_protein_predage_gentle$gapresid[test_indices] <- gap - yfit  
      
    }
    
    
    ###model eval total
    r <- cor(current_protein_predage_gentle$Age, current_protein_predage_gentle$predAge)
    results_eval<-postResample(pred = current_protein_predage_gentle$predAge,obs = current_protein_predage_gentle$Age)
    collect_eval_total<-data.frame(
      organ = organ_base_filename[organ],
      gentle = gentle,
      r = r,
      mae_total = results_eval[3],
      RMSE_total = results_eval[1],
      Rsquared_total = results_eval[2]
    )
    
    
    #####save results
    
    ###find the best model across 5 folds
    
    best_folds_index = which.min(collect_eval$MAE_test)
    
    best_model<-list_models[[best_folds_index]]
    best_params<-list_parameter[[best_folds_index]]
    
    if(gentle==0){
      write.csv(current_protein_predage_gentle,file = paste0(save_path,"F/",organ_base_filename[organ],"_predAge.csv"),row.names = F)
      write.csv(collect_eval,file = paste0(save_path,"F/",organ_base_filename[organ],"_collect_eval.csv"),row.names = F)
      write.csv(collect_eval_total,file = paste0(save_path,"F/",organ_base_filename[organ],"_totalcollect_eval.csv"),row.names = F)
      save(best_model,best_params,file = paste0(save_path,"F/",organ_base_filename[organ],"_bestmodel.RData"))
    }else{
      write.csv(current_protein_predage_gentle,file = paste0(save_path,"M/",organ_base_filename[organ],"_predAge.csv"),row.names = F)
      write.csv(collect_eval,file = paste0(save_path,"M/",organ_base_filename[organ],"_collect_eval.csv"),row.names = F)
      write.csv(collect_eval_total,file = paste0(save_path,"M/",organ_base_filename[organ],"_totalcollect_eval.csv"),row.names = F)
      save(best_model,best_params,file = paste0(save_path,"M/",organ_base_filename[organ],"_bestmodel.RData"))
    }
    
    
  }
  
}


##################### Apply the model trained on normal subjects to other subjects (excluding normal controls)
rm(list = ls())
organ_specific_path<-"./select_protein4_impute_onlymean//"
organ_base_filename<-c("Adrenal","Artery","Brain","Esophagus","Heart","Immune","Intestine","Kidney","Liver","Lung","Muscle","Pancreas","Pituitary","Salivary","Skin","Stomach","Thyroid")
Healthy_group<-fread("./healthy_group_eid_onlyhes.csv")
save_path<-"./elg_results/"

for(organ in c(1:17)){
  print(organ_base_filename[organ])
  collect_eval<-list()
  
  current_protein <- fread(paste0(organ_specific_path,organ_base_filename[organ],"_select_protein_new.csv"))%>%
    subset(eid %in% setdiff(eid,Healthy_group$eid))%>%
    select(-"eid")%>%
    na.omit()%>%
    as.data.frame()
  
  current_protein_nocov<-current_protein%>%
    select(-"Sex",-"Age")
  
  protein_name<-colnames(current_protein_nocov)
  rm(current_protein_nocov)
  current_protein$Sex<-as.factor(current_protein$Sex)
  
  current_protein_predage<-fread(paste0(organ_specific_path,organ_base_filename[organ],"_select_protein_new.csv"))%>%
    subset(eid %in% setdiff(eid,Healthy_group$eid))%>%
    select("eid","Age")
  current_protein_predage$predAge<-NA
  current_protein_predage$gap_resid<-NA
  
  for(gentle in c(0,1)){
    current_protein_gentle<-subset(current_protein,Sex==gentle)
    
    current_protein_predage_gentle<-fread(paste0(organ_specific_path,organ_base_filename[organ],"_select_protein_new.csv"))%>%
      subset(eid %in% setdiff(eid,Healthy_group$eid))%>%
      subset(Sex==gentle)%>%
      select("eid","Age")
    current_protein_predage_gentle$predAge<-NA
    current_protein_predage_gentle$gapresid<-NA
    
    if(gentle==0){
      load(paste0(save_path,"F/",organ_base_filename[organ],"_bestmodel.RData"))
      
      X_data<-as.matrix(current_protein_gentle[,protein_name])
      Y_data<-as.matrix(current_protein_gentle[,"Age"])
      
      predictions_results <- predict(best_model, newx = as.matrix(X_data), s = best_params$lambda)
      current_protein_predage_gentle$predAge<-predictions_results
      #age gap correction
      age_pred<-predictions_results
      age_real<-Y_data
      
      gap <- age_pred - age_real        
      
      b <- glm(gap ~ age_real, family = gaussian)  
      
      yfit <- predict(b, newdata = as.data.frame(age_real), type = "response")  
      current_protein_predage_gentle$gapresid <- gap - yfit  
      
      results_pred_eval<-postResample(pred = current_protein_predage_gentle$predAge,obs = current_protein_predage_gentle$Age)
      r = cor(current_protein_predage_gentle$predAge,current_protein_predage_gentle$Age)
      eval_outofhealth<-data.frame(
        organ_base_filename[organ],
        gentle = gentle,
        r = r,
        MAE = results_pred_eval[3],
        RMSE = results_pred_eval[1],
        Rsquared = results_pred_eval[2]
      )
      
      write.csv(current_protein_predage_gentle,file = paste0(save_path,"F/",organ_base_filename[organ],"_outofH_predage.csv"),row.names = F)
      write.csv(eval_outofhealth,file = paste0(save_path,"F/",organ_base_filename[organ],"_totaloutofH_collect_eval.csv"),row.names = F)
      
    }else{
      load(paste0(save_path,"M/",organ_base_filename[organ],"_bestmodel.RData"))
      
      X_data<-as.matrix(current_protein_gentle[,protein_name])
      Y_data<-as.matrix(current_protein_gentle[,"Age"])
      
      predictions_results <- predict(best_model, newx = as.matrix(X_data), s = best_params$lambda)
      current_protein_predage_gentle$predAge<-predictions_results
      #age gap correction
      age_pred<-predictions_results
      age_real<-Y_data
      
      gap <- age_pred - age_real        
      
      b <- glm(gap ~ age_real, family = gaussian)  
      
      yfit <- predict(b, newdata = as.data.frame(age_real), type = "response")  
      current_protein_predage_gentle$gapresid <- gap - yfit  
      
      results_pred_eval<-postResample(pred = current_protein_predage_gentle$predAge,obs = current_protein_predage_gentle$Age)
      r = cor(current_protein_predage_gentle$predAge,current_protein_predage_gentle$Age)
      eval_outofhealth<-data.frame(
        organ_base_filename[organ],
        gentle = gentle,
        r = r,
        MAE = results_pred_eval[3],
        RMSE = results_pred_eval[1],
        Rsquared = results_pred_eval[2]
      )
      
      write.csv(current_protein_predage_gentle,file = paste0(save_path,"M/",organ_base_filename[organ],"_outofH_predage.csv"),row.names = F)
      write.csv(eval_outofhealth,file = paste0(save_path,"M/",organ_base_filename[organ],"_totaloutofH_collect_eval.csv"),row.names = F)
      
    }
  }
}




################### combined F/M and health/outofH adjusted age gap ###################
rm(list = ls())
organ_specific_path<-"./select_protein4_impute_onlymean//"
organ_base_filename<-c("Adrenal","Artery","Brain","Esophagus","Heart","Immune","Intestine","Kidney","Liver","Lung","Muscle","Pancreas","Pituitary","Salivary","Skin","Stomach","Thyroid")
Healthy_group<-fread("./healthy_group_eid_onlyhes.csv")
save_path<-"./elg_results/"
newsavepath<-"./elg_results_allpheno/"
dir.create(newsavepath)
for(organ in c(1:17)){
  print(organ_base_filename[organ])
  female_health<-fread(paste0(save_path,"F/",organ_base_filename[organ],"_predAge.csv"))
  female_outofhealth<-fread(paste0(save_path,"F/",organ_base_filename[organ],"_outofH_predage.csv"))
  male_health<-fread(paste0(save_path,"M/",organ_base_filename[organ],"_predAge.csv"))
  male_outofhealth<-fread(paste0(save_path,"M/",organ_base_filename[organ],"_outofH_predage.csv"))
  
  female_all<-rbind(female_health,female_outofhealth)%>%
    mutate(Sex = 0)
  
  male_all<-rbind(male_health,male_outofhealth)%>%
    mutate(Sex = 1)
  
  allpredage<-rbind(female_all,male_all)%>%
    arrange(eid)
  
  write.csv(allpredage,paste0(newsavepath,organ_base_filename[organ],"_predage_all.csv"),row.names = F)
}

rm(list = ls())
organ_specific_path<-"./select_protein4_impute_onlymean//"
organ_base_filename<-c("Adrenal","Artery","Brain","Esophagus","Heart","Immune","Intestine","Kidney","Liver","Lung","Muscle","Pancreas","Pituitary","Salivary","Skin","Stomach","Thyroid")
Healthy_group<-fread("./healthy_group_eid_onlyhes.csv")
save_path<-"./elg_results/"
newsavepath<-"./elg_results_allpheno/"
dir.create(newsavepath)
for(organ in c(1:17)){
  print(organ_base_filename[organ])
  
  allpredage<-fread(paste0(newsavepath,organ_base_filename[organ],"_predage_all.csv"))
  
  results_pred_eval<-postResample(pred = allpredage$predAge,obs = allpredage$Age)
  r = cor(allpredage$predAge,allpredage$Age)
  eval_outofhealth<-data.frame(
    organ_base_filename[organ],
    r = r,
    MAE = results_pred_eval[3],
    RMSE = results_pred_eval[1],
    Rsquared = results_pred_eval[2]
  )
  
  write.csv(eval_outofhealth,paste0(newsavepath,organ_base_filename[organ],"_eval_all.csv"),row.names = F)
}



folder_path<-"./elg_results_allpheno/"
file_suffix <- "eval_all.csv" 
file_list <- list.files(folder_path, pattern = file_suffix, full.names = TRUE)
combined_data <- data.frame()
for (file in file_list) {
  data <- read.csv(file, header = TRUE)  
  combined_data <- rbind(combined_data, data)
}

write.csv(combined_data, paste0(folder_path,"combined_eval_all.csv"), row.names = FALSE)
