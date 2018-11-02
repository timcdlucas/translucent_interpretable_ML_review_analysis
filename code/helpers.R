
plotCV <- function(t, print = TRUE){
  stopifnot(inherits(t, 'train'))
  
  row_matches <- sapply(1:length(t$bestTune), function(x) t$pred[, names(t$bestTune)[x]] == t$bestTune[[x]])
  best_rows <- rowMeans(row_matches) == 1

  d <- t$pred[best_rows, ]

  if('weights' %in% names(d)){
    p <- ggplot(d, aes(obs, pred, size = weights))
  } else { 
    p <- ggplot(d, aes(obs, pred))
  }

  p <- p + 
        geom_point(alpha = 0.3) + 
        geom_smooth() +
        geom_abline(slope = 1, intercept = 0)
  
  if(print) print(p)

  return(invisible(p))

} 



compare_models <- function(t1, t2, print = TRUE){
  stopifnot(inherits(t1, 'train'))
  stopifnot(inherits(t2, 'train'))

  row_matches <- sapply(1:length(t1$bestTune), function(x) t1$pred[, names(t1$bestTune)[x]] == t1$bestTune[[x]])
  best_rows <- rowMeans(row_matches) == 1

  d1 <- t1$pred[best_rows, ]

  row_matches <- sapply(1:length(t2$bestTune), function(x) t2$pred[, names(t2$bestTune)[x]] == t2$bestTune[[x]])
  best_rows <- rowMeans(row_matches) == 1

  d2 <- t2$pred[best_rows, ]

  d <- left_join(d1, d2, by = 'rowIndex')

  if('weights' %in% names(d1)){
    p <- ggplot(d, aes(pred.x, pred.y, size = weights))
  } else { 
    p <- ggplot(d, aes(pred.x, pred.y))
  }

  corre <- round(cor(d$pred.x, d$pred.y), 2)
  r1 <- round(min(t1$results[, t1$metric]), 2)
  r2 <- round(min(t2$results[, t2$metric]), 2)
  p <- p + 
        geom_point(alpha = 0.3) + 
        geom_abline(slope = 1, intercept = 0) +
        labs(x = t1$method, y = t2$method) +
        ggtitle(paste0('Correlation: ', corre, '. t1 ', t1$metric, ': ', r1, '. t2 ', t1$metric, ': ', r2))
  
  if(print) print(p)

  return(invisible(p))

} 


