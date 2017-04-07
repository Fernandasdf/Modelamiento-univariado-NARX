library(e1071)
inverse.step.test <- function(
  model,
  inverse.step.lags,
  lag
){
  x <- inverse.step.lags$folded.signal["MABP"]
  y <- inverse.step.lags$folded.signal["CBFV"]
  
  mean.cbfv <- mean(y[,])
  
  abp.step.signal <- data.frame(
    MABP = x,
    CBFV = round(mean.cbfv,4)
  )
  
  response.ARX.model <- NULL
  step.response <- NULL
  
  lag.abp.step <- lag.abp(abp.step.signal, lag)
  columns <- sum(unlist(lag)) + 1
  
  for (indice in 1:nrow(lag.abp.step)){
    step.response <- predict(model, lag.abp.step[indice,1:columns])
    response.ARX.model[indice] <- round(step.response,4)
    #
    if(lag[["CBFV"]]>1){
      for(b in lag[["CBFV"]]:2){
        var.name.old <- paste0("CBFV.lag",b)
        var.name.new <- paste0("CBFV.lag",b-1)
        lag.abp.step[indice+1,var.name.old] <- lag.abp.step[indice,var.name.new]
      }
      lag.abp.step[indice+1,"CBFV.lag1"] <- round(step.response,4)
    }else{
      var.name <- "CBFV.lag1"
      lag.abp.step[indice+1,var.name] <- round(step.response,4)
    }
  }
  #
  #
  start <- max(unlist(lag)) +1
  y <-  y[start:length(y$CBFV),1]
  pred <- response.ARX.model
  pred <- round(pred, 3)
  test.cor.step <<- cor(response.ARX.model,y)
  print(test.cor.step)
  ts = 0.6
  before.drop = 10.2 #seconds
  step.duration = 42 #seconds
  step.length = step.duration/ts
  #fin = 300*time.sample
  half.step = floor(length(inverse.step[,1])/2)
  ajuste = pred[half.step] - 0.8
  
  pred <- pred[(half.step-before.drop/ts):(half.step-before.drop/ts+step.length)]-ajuste
  
  time <- seq(1,length(pred),0.6)
  pred<-pred-ajuste
  df <- data.frame(pred,time[1:length(pred)])
  #colnames(df) <- c("VFSC estimado (cm/seg)","Tiempo (seg)")
  title <- paste('Step response subject: ',name.subject)
  plot<-ggplot(df, aes(x = time.1.length.pred.., y = pred)) +
    geom_line() +
    xlab("Tiempo (seg)") + ylab("VFSC estimado (cm/seg)") +
    ggtitle(title) 
  print(plot)
  ggsave(paste('C:/Users/Feffyta/Documents/Universidad/tesis/Programas/Entrenamiento en R/Modelamiento-univariado-NARX/Results/Univariado/NARX/stepResponse_',name.subject,'.jpg'), width = 8,height = 7)
  
  
}

make.inverseStep <- function(
  lag
){
  ts = 0.6 #time sample
  length.step = 150/ts #largo del escalon inverso
  inverse.step.unos <- matrix(1,floor(length.step/2),1) #parte de 1's del escalon
  inverse.step.ceros <- matrix(0,floor(length.step/2),1) #parte de 0's del escalon
  inverse.step <- rbind(inverse.step.unos,inverse.step.ceros) #escalon inverso 
  butter.filter <- butter(2,0.3)
  inverse.step <- filter(butter.filter,inverse.step)
  inverse.step = data.frame(inverse.step)
  
  
  
  
  inverse.step <- cbind(inverse.step,matrix(0,length.step,1))
  
  colnames(inverse.step) <- c("MABP","CBFV")
  
  signal.step <- retardos_multi_inverseStep(inverse.step, lag)
  
}

lag.abp <- function(
  abp.step.signal,
  lags
)
{
  max.lag <- max(unlist(lags)) + 1
  indices <- 1:nrow(abp.step.signal)
  lag.mat <- embed(indices, max.lag)
  
  col.names <- list("MABP","CBFV")
  columns <- NULL
  lagged.columns.names <- c()
  for(colname in col.names){
    lag.order <- lags[[colname]]
    if(!is.null(lag.order) && lag.order > 0)
      for(i in lag.order:1){
        new.colname <- paste(colname, paste0("lag", i), sep = ".")
        lagged.columns.names <- c(lagged.columns.names, new.colname)
        columns[[new.colname]] <- abp.step.signal[lag.mat[, i + 1], colname]
      }
    columns[[colname]] <- abp.step.signal[lag.mat[, 1], colname]
  }
  X<-data.frame(columns)
  sorting <- order(lag.mat[, 1])
  X <- X[sorting, ]
}


training <- function(
  parameter,
  lag,
  fold,
  input.var.names,
  output.var.names,
  signal.train,
  signal.test,
  inputs,
  inverse.step.lags
)
{
  fmla.str <- paste(inputs, collapse = " + ")
  fmla.str <- paste(output.var.names, "~", fmla.str)
  fmla <- formula(fmla.str)
  params <- list(formula = fmla, data = signal.train$folded.signal, scale = FALSE)
  params <- c(params, list(type = "nu-regression", kernel = "radial"))
  params <- c(params, parameter)
  
  start.time <- Sys.time()
  model <- do.call(svm, params)
  end.time <-Sys.time()
  
  stats <- eval.model(
    parameter,
    lag,
    fold,
    input.var.names,
    output.var.names,
    signal.train,
    signal.test,
    inputs,
    model,
    inverse.step.lags
  )
  results <- list(models = list(model), stats = stats) 
}
eval.model <- function(
  parameter,
  lag,
  fold,
  input.var.names,
  output.var.names,
  signal.train,
  signal.test,
  inputs,
  model,
  inverse.step.lags
){
  fitted.signal <- round(model[["fitted"]], 4)
  train.cor <- cor(fitted.signal, signal.train$folded.signal["CBFV"])
  train.cor <- round(train.cor, 4)
  
  x <- signal.test$folded.signal["MABP"]
  y <- signal.test$folded.signal["CBFV"]

  mean.cbfv <- mean(y[,])

  abp.step.signal <- data.frame(
    MABP = x,
    CBFV = round(mean.cbfv,4)
  )
  
   # print( abp.step.signal )
  #
      lag.abp.step <- lag.abp(abp.step.signal, lag)

   # print(lag.abp.step)
    lag.abp.step <- lag.abp.step[inputs]
  #
    columns <- sum(unlist(lag)) + 1
  #
    response.NARX.model <- NULL

 # pred <- predict(model, x)
 # pred <- round(pred, 3)
 # test.cor <- cor(pred, y)
 # test.cor <- round(test.cor, 3)
   for (indice in 1:nrow(lag.abp.step)){
     step.response <- predict(model, lag.abp.step[indice,1:columns])
     response.NARX.model[indice] <- round(step.response,4)
  #
     if(lag[["CBFV"]]>1){
       for(b in lag[["CBFV"]]:2){
         var.name.old <- paste0("CBFV.lag",b)
         var.name.new <- paste0("CBFV.lag",b-1)
         lag.abp.step[indice+1,var.name.old] <- lag.abp.step[indice,var.name.new]
       }
       lag.abp.step[indice+1,"CBFV.lag1"] <- round(step.response,4)
     }else{
       var.name <- "CBFV.lag1"
       lag.abp.step[indice+1,var.name] <- round(step.response,4)
     }
   }
  #
  #
   start <- max(unlist(lag)) +1
   y <-  y[start:length(y$CBFV),1]
  
   test.cor <- cor(response.NARX.model,y)

  
  
  print(test.cor)  

  inverse.step.result <- inverse.step.test(model,inverse.step.lags,lag)
  
  data.frame(
    MABP = lag["MABP"],
    MABP = lag["CBFV"],
    gamma = parameter['gamma'],
    nu = parameter['nu'],
    cost = parameter['cost'],
    train.cor = train.cor,
    test.cor = test.cor,
    fold = fold,
    test.cor.step = test.cor.step
  )
}
retardos_multi <- function(
  arch_entrada,
  lags
)
{
  signal <- read.table(arch_entrada, col.names = c("CBFV","MABP","PIC") )
  signal[['CBFV']] <- round(signal[['CBFV']],4)
  signal[['MABP']] <- round(signal[['MABP']],4)
  signal[['PIC']] <- round(signal[['PIC']],4)

  signal.uni <- data.frame(
    MABP = signal[['MABP']],
    CBFV = signal[['CBFV']]
  )
  max.lag <- max(unlist(lags)) + 1
  indices <- 1:nrow(signal.uni)
  lag.mat <- embed(indices, max.lag)
  
  col.names <- list("MABP","CBFV")
  columns <- NULL
  lagged.columns.names <- c()
  for(colname in col.names){
    lag.order <- lags[[colname]]
    columns[[colname]] <- signal.uni[lag.mat[, 1], colname]
    if(!is.null(lag.order) && lag.order > 0)
      for(i in 1:lag.order){
        new.colname <- paste(colname, paste0("lag", i), sep = ".")
        lagged.columns.names <- c(lagged.columns.names, new.colname)
        columns[[new.colname]] <- signal.uni[lag.mat[, i + 1], colname]
      }
    
  }
  folded.signal <- data.frame(columns)
  sorting <- order(lag.mat[, 1])
  folded.signal <- folded.signal[sorting, ]
  list(folded.signal = folded.signal, lagged.columns.names = lagged.columns.names)
}

retardos_multi_inverseStep <- function(
  arch_entrada,
  lags
)
{
  #signal <- read.table(arch_entrada, col.names = c("MABP","CBFV") )
  arch_entrada[['CBFV']] <- round(arch_entrada[['CBFV']],4)
  arch_entrada[['MABP']] <- round(arch_entrada[['MABP']],4)
  signal.uni <- data.frame(
    MABP = arch_entrada[['MABP']],
    CBFV = arch_entrada[['CBFV']]
  )  
  max.lag <- max(unlist(lags)) + 1
  indices <- 1:nrow(signal.uni)
  lag.mat <- embed(indices, max.lag)
  
  col.names <- list("MABP","CBFV")
  columns <- NULL
  lagged.columns.names <- c()
  for(colname in col.names){
    
    lag.order <- lags[[colname]]
    columns[[colname]] <- signal.uni[lag.mat[, 1], colname]
    if(!is.null(lag.order) && lag.order > 0)
      for(i in 1:lag.order){
        new.colname <- paste(colname, paste0("lag", i), sep = ".")
        lagged.columns.names <- c(lagged.columns.names, new.colname)
        columns[[new.colname]] <- signal.uni[lag.mat[, i+1], colname]
      }
    
    
    
  }
  folded.signal <- data.frame(columns)
  
  sorting <- order(lag.mat[, 1])
  folded.signal <- folded.signal[sorting, ]
  list(folded.signal = folded.signal, lagged.columns.names = lagged.columns.names)
}

get.results <- function(
    parametros,
    lags,
    src.basename,
    src.dir,
    src.ext,
    keep.nstats
    
  )
{
  parameter <- expand.grid(parametros)
  folds <- list(fold1 = 'p1', fold2 = 'p2')
  input.var.names <- c("MABP")
  output.var.names <- c('CBFV')
  
  
  archivo1 <- paste(src.basename,folds[['fold1']],sep="_")
  archivo1 <- paste(archivo1,src.ext,sep=".")
  
  archivo2 <- paste(src.basename,folds[['fold2']],sep="_")
  archivo2 <- paste(archivo2,src.ext,sep=".")
  
  file.fold1 <- paste(src.dir,archivo1, sep="/")
  file.fold2 <- paste(src.dir,archivo2, sep="/")
  
  lag <- list(MABP = lags[,1][["MABP"]], CBFV = lags[,1][["CBFV"]])
  fold <-lags[,1][["fold"]]
  if(fold == 1){    
    signal.train <- retardos_multi(file.fold1, lag)
    signal.test <- retardos_multi(file.fold2, lag)
  }else{
    signal.train <- retardos_multi(file.fold2, lag)
    signal.test <- retardos_multi(file.fold1, lag)
  }
  
  
  inputs <- c(signal.train[["lagged.columns.names"]],input.var.names)
  inverse.step.lags <- make.inverseStep(lag) #obtiene la matriz de retrasos a partir del escalon inverso
  output = apply(parameter, 1,function(p) training(p,lag,fold,input.var.names,output.var.names,signal.train,signal.test,inputs,inverse.step.lags))
  
  models <- sapply(output, function(l) l[[1]])
  stats.list <- lapply(output, function(l) l[[2]])
  stats <- do.call(rbind, stats.list)
  colnames(stats) <- c("MABP","CBFV","gamma","nu","cost","train.cor","test.cor","fold","test.cor.step")
  
  i <- order(
    -stats["test.cor"]
  )
  stats <- stats[i, ]
  if(nrow(stats) > keep.nstats)
  {
    i <- 1:keep.nstats
    stats <- stats[i, ]
  }
  
  rownames(stats) <- NULL
  
  stats
}