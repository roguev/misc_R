# creates an AB - BA scatter plot for a particular score matrix
scatter_replicates = function(scores,alpha=0.5,paired_data=NULL) {
  require(ggplot2);
  if (is.null(paired_data)) {
    paired_data = get_paired_data(scores);
  }
  
  if (is.null(paired_data)) {
    print("No paired data found");
    return(NULL);
  }
  
  r = range(paired_data[,"data1"],paired_data[,"data2"]);
  min_score = r[1];
  max_score = r[2];
  result = data.frame(replicate_1 = paired_data[,"data1"], replicate_2 = paired_data[,"data2"]);
  cc = cor(x = result$replicate_1, y = result$replicate_2);
  # do the plotting
  d = ggplot(result,aes(x=replicate_1,y=replicate_2)) +
    geom_point(color="blue", size=1,alpha=1) +
    geom_smooth(method=lm,se=F,linetype=2,color="red") +
    stat_density2d(geom="polygon",aes(fill=..level..),contour=T,alpha=alpha) +
    geom_density2d(colour="black", size = 0.5) + 
    scale_fill_gradient(low="yellow",high="red") + 
    xlim(min_score,max_score) +
    ylim(min_score,max_score) + 
    ggtitle(round(cc,digits = 2));
  print(d);
  print(sprintf("Pearson cc: %f", cc), quote = FALSE);
  return(list("df" = result, "plot" = d));
}

# creates a pretty-looking scatterplot
scatter_with_density = function(...,setX = 1, setY = 1, lX = "dataX", lY = "dataY",alpha=0.2, contour_thickness = 0.5, with_plot = TRUE, title = NULL, square_axes = F) {
  require(ggplot2);
  dots = list(...);
  dX = as.matrix(dots[[1]]);
  dY = as.matrix(dots[[2]]);
  df = data.frame(dX = dX[,setX],dY = dY[,setY]);
  rangeX = range(df$dX);
  rangeY = range(df$dY);
  
  # do the plotting
  d = ggplot(df, aes(x=dX,y=dY)) +
    geom_point(color="black", size=.75) +
    stat_density2d(aes(fill=..level..), geom="polygon",contour=T,alpha=alpha) +
    geom_density2d(colour = "blue", size = contour_thickness, alpha = 1) +
    scale_fill_gradient(low="yellow",high="red",guide="none") + 
    xlab(as.character(lX)) + ylab(as.character(lY)) + 
    ggtitle(as.character(title));
  
  if (square_axes == T) {
    data_min = min(rangeX, rangeY);
    data_max = max(rangeX,rangeY);
    d = d + xlim(data_min, data_max) + ylim(data_min, data_max);
  } else {
    d = d + xlim(rangeX[1],rangeX[2]) + ylim(rangeY[1],rangeY[2]); 
  }
  
  if(with_plot == TRUE) {
    print(d);
  }
  return(d);
}

# scatters several datasets
scatter_scores = function(..., alpha = 0.2, labels = NULL, with_plot = TRUE) {
  datasets = make_consistent_datasets(...);
  ggplots = list();
  for (i in 1:length(datasets)) {
    for (j in 1:length(datasets)) {
      ggplots[[length(ggplots) + 1]] = scatter_with_density(datasets[[i]]$data, datasets[[j]]$data, with_plot = FALSE, alpha = alpha, contour_thickness = 0.2, lX = labels[i], lY = labels[j], title = paste("CC = ",as.character(round(cor(datasets[[i]]$data, datasets[[j]]$data), digits = 2))));
    }
  }
  if (with_plot == TRUE) {
    multiplot(plotlist = ggplots, cols = length(datasets));
  }
  
  return(ggplots);
}

scatter_scores_with_facets = function(..., labels = NULL, with_plot = TRUE, alpha = 0.2, contour_thickness = 0.2) {
  # generate labels if they are not present
  scores = make_consistent_datasets(...);
  if (is.null(labels)) {
    labels = mapply(paste,rep("data_",length(scores)),c(1:length(scores)),sep="");
  }
  
  df = data.frame(scoremat1 = NULL, scoremat2 = NULL, dataX = NULL, dataY = NULL);
  for (i in 1:length(scores)) {
    for (j in 1:length(scores)) {
      if (!i == j) {
        df = rbind(df,data.frame(scoremat1 = labels[[i]], scoremat2 = labels[[j]], dataX = scores[[i]]$data, dataY = scores[[j]]$data));
      } else {
        df = rbind(df, data.frame(scoremat1 = labels[[i]], scoremat2 = labels[[j]], dataX = 0, dataY = 0));
      }
    }
  }
  
  require(ggplot2);
  d = ggplot(df, aes(x=dataX, y=dataY)) + 
    geom_point(color="black", size=.5, alpha = 0.05) +
    #stat_density2d(aes(fill=..level..), geom="polygon",contour=T,alpha=alpha, na.rm=T) +
    #geom_density2d(colour = "blue", size = contour_thickness, alpha = 1, na.rm=T) +
    #scale_fill_gradient(low="yellow",high="red", guide="none") +
    xlab(NULL) + ylab(NULL) +
    facet_grid(scoremat1 ~ scoremat2, scales="free");
  
  if (with_plot == TRUE) {
    print(d);
  }
  
  return(d);
}

# a wrapper around the scatter_scores_with_facets to save some typing
scatter_score_set = function(score_set) {
  datasets = list();
  for (i in 1:length(score_set)) {
    if ("scores" %in% names(score_set[[i]])) {
      dataset = score_set[[i]]$scores;
    } else {
      dataset = score_set[[i]];
    }
    datasets[[length(datasets)+1]] = dataset;
  }
  d = scatter_scores_with_facets(datasets, labels = names(score_set));
  return(d);
}

# draws a heatmap clustergram of a scores matrix
draw_cluster = function(scores, Qs = NULL, As = NULL, pal = colorRampPalette(c("blue", "white", "red"))(n = 256)) {
  if (is.matrix(scores)) {
    return(heatmap(scores, col = pal));
  } else {
    return(heatmap(reformat_to_matrix(scores, Qs = Qs, As = As),col = pal));
  }
}