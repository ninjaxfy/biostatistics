#######################
#' @Classfication_task ###

library(seqR)
library(caret)
library(Biostrings)
library(e1071)
library(MASS)
library(randomForest)

#inpu_file
pos_fasta <- as.data.frame(readBStringSet(filepath = 'pos.fa'))
neg_fasta <- as.data.frame(readBStringSet(filepath = 'neg.fa'))
colnames(pos_fasta) <- 'seq'
colnames(neg_fasta) <- 'seq'
pos_fasta$label <- 'pos'
neg_fasta$label <- 'neg'

result <- list()
result1 <- list()
for (k in c(1:6)) {
  ##get k-mer
  message(paste0("argument k=",k,' ',Sys.time()))
  pos_x <- as.matrix(seqR::count_kmers(sequences = pos_fasta$seq,k = k))
  neg_x <- as.matrix(seqR::count_kmers(sequences = neg_fasta$seq,k = k))
  
  colnames(pos_x) <- gsub("0",'',colnames(pos_x))
  colnames(pos_x) <- gsub("[.]",'',colnames(pos_x))
  colnames(pos_x) <- gsub("_",'',colnames(pos_x))
  pos_x <- as.data.frame(pos_x)
  row.names(pos_x) <- NULL
  pos_x$label <- 'pos'
  
  colnames(neg_x) <- gsub("0",'',colnames(neg_x))
  colnames(neg_x) <- gsub("[.]",'',colnames(neg_x))
  colnames(neg_x) <- gsub("_",'',colnames(neg_x))
  neg_x <- as.data.frame(neg_x)
  neg_x$label <- 'neg'
  neg_x <- neg_x[colnames(pos_x)]
  row.names(neg_x) <- NULL
  
  set.seed(1)
  pos_train <- sample(1:nrow(pos_x),nrow(pos_x)/2)
  set.seed(2)
  neg_train <- sample(1:nrow(neg_x),nrow(neg_x)/2)
  
  train <- rbind(pos_x[pos_train,],neg_x[neg_train,])
  train$label <- as.factor(train$label)
  test <- rbind(pos_x[-pos_train,],neg_x[-neg_train,])
  test$label <- as.factor(test$label)
  
  ###logictic regression
  lr.fit <- glm(label~.,data = train,family = 'binomial')
  lr.prod <- predict(lr.fit,test,type = 'response')
  lr.pred <- rep('neg',nrow(test))
  lr.pred[which(lr.prod>=0.5)] <- 'pos'
  ROC <- roc(test$label,lr.prod)
  result[[paste0('LR.',k,'-mer')]] <- confusionMatrix(table(test$label,lr.pred))
  result1[[paste0('LR.',k,'-mer')]] <- ROC
  
  ###SVM
  tune.out <- tune(method = svm,
                   label~.,data = train,
                   kernel='radial',decision.value=T,
                   ranges = list(cost = c(0.1,1,10,100), 
                                 gamma = c(0.001,0.1,1,5)))
  
  svm.best.model <- tune.out$best.model
  svm.pred <- predict(svm.best.model,test,decision.values = T)
  fitted <- attributes(svm.pred)$decision.value
  t <- roc(test$label,fitted)
  
  ##Random Forest
  ntree = 50
  mtry <- round(ncol(train[,-1])^(1/2))
  rf.fit <- randomForest(x.train,y.train,ntree = ntree,mtry = mtry,importance = T)
  
  rf.predict <- predict(rf.fit,x.test,decision.values = T)
  acc.res <- table(y.test,rf.predict)
  acc <- (acc.res[1,1]+acc.res[2,2])/sum(acc.res)
  result1[k] <- acc
  
  ###LDA
  lda.fit <- lda(label~.,data = train)
  lda.prod <- predict(lda.fit,test)
  lda.pred <- lda.prod$class
  lda.prosterior <- lda.prod$posterior[,2]
  result1[[paste0('LDA.',k,'-mer')]] <- confusionMatrix(table(lda.pred,test$label))
  ROC <- roc(test$label,lda.prosterior)
  result2[[paste0('LDA.',k,'-mer')]] <- ROC
  
  
}

save.image(file = 'cla.RData')
