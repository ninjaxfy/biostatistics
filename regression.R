#######################
#' @Regression_task ###

library(glmnet)
##read fasta file
fasta <- as.data.frame(readBStringSet(filepath = 'x-2k_sequence.fa'))
colnames(fasta) <- 'seq'
fasta$Geneid <- row.names(fasta)
row.names(fasta) <- NULL

##read expression file
tpm <- read.csv(file = 'y.csv',header = T)
rownames(tpm) <- tpm$Geneid
tpm <- tpm['TPM']

data_name <- fasta$Geneid

result <- list()
for (k in c(1:6)) {
  ##get k-mer
  message(paste0("argument k=",k))
  x <- as.matrix(seqR::count_kmers(sequences = fasta$seq,k = k))
  colnames(x) <- gsub("0",'',colnames(x))
  colnames(x) <- gsub("[.]",'',colnames(x))
  colnames(x) <- gsub("_",'',colnames(x))
  x <- x[,-grep('N',colnames(x))]
  rownames(x) <- fasta$Geneid
  
  set.seed(1)
  train <- sample(1:length(data_name),length(data_name)/2)
  
  x.train <- x[train,]
  y.train <- tpm[train,]$TPM
  x.test <- x[-train,]
  y.test <- tpm[-train,]$TPM
  
  ##Linear Regression
  
  
  ##ridge
  Ridge.model.cv <- cv.glmnet(x.train,y.train,alpha=0)
  best.lambda <- Ridge.model.cv$lambda.min
  ##re-train
  Ridge.best.model <- glmnet(x.train,y.train,alpha = 0, lambda = best.lambda)
  y.test.pred <- predict(Ridge.best.model, newx = x.test, s=best.lambda)
  PCC <- cor.test(y.test,y.test.pred)
  
  result[[paste0('Ridge.',k,'-mer')]] <- PCC
  
  ##LASSO
  lasso.cv.fit <- cv.glmnet(data.matrix(x.train),y.train,alpha=1)
  best_lambda <- lasso.cv.fit$lambda.min
  lasso.best.fit <- glmnet(data.matrix(x.train),y.train,alpha=1,lambda = best_lambda)
  y.test.predict <- predict(lasso.best.fit,data.matrix(x.test),s=best_lambda)
  lasso.pcc <- cor(y.test,y.test.predict)
  lasso.mse <- mean((y.test.predict-y.test)^2)
  result[[paste0('lasso',k,'-mer')]] <- lasso.pcc
  
  
}

