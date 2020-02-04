# reads a data file into a dataframe for further analysis
dataset_load = function(file, field_sep = "\t", skip  = 1) {  
  data = read.table(file, sep = field_sep, skip = skip);
  # modified 20150114
  colnames(data)[1] = "Q_ID"
  colnames(data)[2] = "G1";
  colnames(data)[3] = "G2";
  D = dim(data);
  colnames(data)[4:D[2]] = sub("V","",colnames(data)[4:D[2]]);
  colnames(data)[4:D[2]] = as.numeric(colnames(data)[4:D[2]]) - 3;
  colnames(data)[4:D[2]] = paste("replicate_",colnames(data)[4:D[2]],sep="");
  data$G1 = toupper(data$G1);
  data$G2 = toupper(data$G2);
  data$G1 = gsub("[^a-zA-Z0-9_]","",data$G1); # remove any special characters
  data$G2 = gsub("[^a-zA-Z0-9_]","",data$G2);
  data$G1 = as.factor(data$G1);
  data$G2 = as.factor(data$G2);
  print(sprintf("Loaded %s, %i pairs, %i replicates",file,D[1],D[2]-3),quote = FALSE);
  return(data);
}

# computes estimated single effects as the medians of all the gene pairs containing particular gene
calc_singles = function(dataN, normalize = TRUE) {
  all_genes = union(dataN$G1,dataN$G2);
  
  # normalize data?
  if (normalize == TRUE) {
    dataN = normalize_data(dataN);
  }
  
  # compute singles
  singles = data.frame(gene=all_genes,data=as.numeric(rep(NA,length(all_genes))));
  L = length(colnames(dataN));
  print(sprintf("Calculating singles for %i genes in %i replicates",length(all_genes), L - 4), quote = FALSE);
  for (i in 1:length(singles$gene)) {
    iq = find_string(singles$gene[i], dataN$G1);
    ia = find_string(singles$gene[i],dataN$G2);
    ii = union(ia,iq);
    singles$data[i] = median(unlist(dataN[ii,5:L]));
    singles$sd[i] = sd(unlist(dataN[ii,5:L]));
  }
  return(singles);
}

# normalizes data by median cenetring the dataset. allows for shifting the center of the distribution
normalize_data = function(data, shift_center = 0) {
  D = dim(data);
  bs = unique(data$batch);
  for (i in 1:length(bs)) {
    inds = which(data$batch == bs[i]);
    print(sprintf("Normalizing batch %i, %i pairs, %i replicates, shift_center = %i",bs[i], length(inds), D[2] - 4, shift_center),quote = FALSE);
    data[inds,5:D[2]] =  data[inds,5:D[2]]/median(unlist(data[inds,5:D[2]])) + shift_center;
  }
  return(data);
}

# computes GI scores usin a multiplicative model based on single effects estimates
calc_multiplicative_scores = function(dataN, singles = NULL, normalize = TRUE, shift_center = 0, change_signs = FALSE) {
  D = dim(dataN);
  
  all_genes = union(dataN$G1, dataN$G2);
  
  # normalize ?
  if (normalize == TRUE) {
    dataN = normalize_data(dataN, shift_center = shift_center);
  }
  
  # compute singles ?
  if (is.null(singles)) {
    singles = calc_singles(dataN, normalize = FALSE);
  }
  
  print(sprintf("Computing multiplicative scores for %i genes", length(union(dataN$G1, dataN$G2))), quote = FALSE);
  scores = data.frame(G1=NULL, G2=NULL, data=NULL);
  for (i in 1:length(all_genes)) {
    df1 = scatter_singles_vs_doubles(all_genes[i], dataN, with_plot = FALSE, singles = singles)$df;
    g_ind = find_string(all_genes[i],singles$gene);
    g_single = singles$data[g_ind];
    expected = df1$singles*g_single;
    observed = df1$doubles;
    if (change_signs == TRUE) {
      df2 = data.frame(G1 = as.character(df1$G1), G2 = as.character(df1$G2), data = abs(expected-observed));
      
      # decide which scores are to become negative
      change_sign = which((expected < g_single & observed < expected) | (expected > g_single & observed > expected));
      df2$data[change_sign] = -1*df2$data[change_sign];
    } else {
      df2 = data.frame(G1 = as.character(df1$G1), G2 = as.character(df1$G2), data = expected-observed);
    }
    scores = rbind(scores,df2);
  }
  return(scores);
}

# a wrapper to score multiple datasets
# uses the regresion score by default
score_data = function(..., shift_center = 0, mode = "regression", model = "loess", data_names = NULL) {
  dots = list(...);
  if (length(shift_center) == 1) {
    shift_center = c(rep(shift_center,length(dots)));
  }
  
  if (length(mode) == 1) {
    mode = c(rep(mode, length(dots)));
  }
  
  if (length(model) == 1) {
    model = c(rep(model,length(dots)));
  }
  
  result = list();
  for (i in 1:length(dots)) { 
    if (is.character(dots[[i]])) {
      d = dataset_load(as.character(dots[[i]]));
    } else {
      d = dots[[i]];
    }
    
    if (is.null(data_names[i]) & is.character(dots[[i]])) {
      data_name = tail(strsplit(as.character(dots[i]),'/')[[1]],n=1);
    }
    
    if (is.null(data_names[i]) & !is.character(dots[[i]])) {
      data_name = paste("data_",as.character(i),"_",model[i],sep="");
    } else {
      data_name = data_names[i];
    }
    
    switch(mode[i],
           multiplicative = {
             result[[i]] = calc_multiplicative_scores(d, shift_center = shift_center[i]);
           }, regression = {
             result[[i]] = calc_regression_scores(d, shift_center = shift_center[i], model = model[i]);
           });
    print(sprintf("Dataset name: %s",data_name), quote = FALSE);
    names(result)[i] = data_name;
  }
  return(result);
}

correlate_clusters = function(...,names=NULL,m = "spearman",mode = "by_row") {
  dots = list(...);
  result = matrix(NA, nrow = length(dots), ncol = length(dots));
  for (i in 1:length(dots)) {
    for (j in 1:length(dots)) {
      if (!i == j) {
        switch(mode,
               by_row = {
                 result[i,j] = cor(dots[[i]]$rowInd, dots[[j]]$rowInd, method = m);
               }, by_col = {
                 result[i,j] = cor(dots[[i]]$colInd, dots[[j]]$colInd, method = m);
               });
      }
    }
  }
  if (is.null(names)) {
    rownames(result) = c(1:length(dots));
    colnames(result) = c(1:length(dots));
  } else {
    rownames(result) = names;
    colnames(result) = names;
  }
  return(result);
}

# extracts all the pairs containing a particular gene
extract_pairs = function(gene,data) {
  qi = find_string(gene,data$G1);
  ai = find_string(gene,data$G2);
  
  if (length(c(qi,ai)) == 0) {
    return(NULL);
  }
  
  G1 = c(as.character(data$G1[qi]),as.character(data$G2[ai]));
  G2 = c(as.character(data$G2[qi]),as.character(data$G1[ai]));
  d = calc_medians(data[c(qi,ai),5:dim(data)[2]]);
  return(data.frame(G1 = G1, G2 = G2, data = d));
}

# creates a scatter plot of measurments of single vs double knockdowns
# used in the regression score calculation. draws a scatter by default
scatter_singles_vs_doubles = function(gene, data, with_plot = TRUE, singles = NULL) {
  if (is.null(singles)) {
    singles = calc_singles(data,normalize = FALSE);
  }
  
  # extract all the pairs for the gene
  pairs = extract_pairs(gene,data);
  if (is.null(pairs)) {
    print(sprintf("No pairs for %s found in doubles data.",gene), quote = FALSE);
    return(NULL);
  }
  
  x = intersect_w_indeces(singles$gene,pairs$G2);
  if (is.null(x)) {
    print(sprintf("No singles for %s found in singles data.",gene), quote = FALSE);
    return(NULL);
  }
  
  df = data.frame(G1 = rep(gene, length(x$i)),G2 = x$i, singles = singles$data[x$ia], doubles = pairs$data[x$ib]);
  require(ggplot2);
  d = ggplot(df,aes(x = singles, y = doubles)) + geom_point() + stat_smooth(method="loess", formula = y~x) + ggtitle(as.character(gene));
  if (with_plot) {
    print(d);
  }
  return(list("df"=df,"plot"=d));
}

# makes a boxplot of the data broken by gene
plot_phenotype_estimates = function(data = data, singles = NULL, with_plot = TRUE, gene_list = NULL) {
  if (is.null(singles)) {
    singles = calc_singles(data,normalize = FALSE);
  }
  
  if (is.null(gene_list)) {
    gene_list= union(data$G1, data$G2);
  }
  
  df_G1 = NULL;
  df_G2 = NULL;
  df_data = NULL;
  for (i in 1:length(gene_list)) {
    tmp = scatter_singles_vs_doubles(gene=gene_list[i], data=data, with_plot = F, singles = singles);
    df_G1 = c(df_G1, as.character(tmp$df$G1));
    df_G2 = c(df_G2, as.character(tmp$df$G2));
    df_data = c(df_data,as.numeric(tmp$df$doubles));
  }
  df = data.frame(G1 = df_G1, G2 = df_G2, data = df_data);
  p = ggplot(df,aes(G1,data)) + geom_boxplot(varwidth=T, aes(fill=factor(G1))) + guides(fill=F);
  
  if (with_plot) {
    print(p);
  }
return(list("df"=df, "plot" = p));
}


# computes GI scores using loess regression between estimated singles and observed doubles 
calc_regression_scores = function(dataN, singles = NULL, normalize = TRUE, shift_center = 0, model = "loess", change_signs = FALSE) {
  all_genes = union(dataN$G1,dataN$G2);
  
  # normalize ?
  if (normalize == TRUE) {
    dataN = normalize_data(dataN, shift_center = shift_center);
  }
  
  # compute singles ?
  if (is.null(singles)) {
    singles = calc_singles(dataN, normalize = FALSE);
  }
  
  # gets the median of the normalized distribution, used for setting the sign to the scores
  data_median = median(as.matrix(dataN[,5:dim(dataN)[2]]));
  
  # initialize stuff
  scores = data.frame(G1 = NULL, G2 = NULL, data = NULL);
  models = list();
  intersects = data.frame(gene = NULL, data = NULL);
  
  print(sprintf("Computing regression scores for %i genes using %s stat model",length(all_genes),model),quote = FALSE);
  for (i in 1:length(all_genes)) {
    # get the singles vs doubles data
    print(sprintf("%s",all_genes[i]),quote = FALSE);
    df1 = scatter_singles_vs_doubles(all_genes[i], dataN, with_plot = FALSE, singles = singles)$df;
    
    switch(model,
           loess = {
             # loess smoothing
             l = loess(doubles ~ singles,df1);
           },lm = {
             # simple linear model
             l = lm(doubles ~ singles,df1);
           }, rlm = {
             # robust linear
             require(MASS);
             l = rlm(doubles ~ singles,df1);
           });
    
    intersect = predict(l,data.frame(singles = data_median));
    
    # for some reason the current implementation gives a really strange and uncorrelated results
    if (change_signs == TRUE) {
      df2 = data.frame(G1 = as.character(df1$G1), G2 = as.character(df1$G2), data = abs(l$residuals));
      
      # decide on the sign of the scores
      expected = l$fitted;
      observed = df1$doubles;
      
      # synergisms
      change_sign = which((expected < intersect & observed < expected) | (expected > intersect & observed > expected));
      
      df2$data[change_sign] = -1*df2$data[change_sign];
    } else {
      # do not change signs, leave the scores as they were derived from the residuals
      # 20150115 change removed -1* in data = -1*l$residuals
      #df2 = data.frame(G1 = as.character(df1$G1), G2 = as.character(df1$G2), data = -1*l$residuals);
      df2 = data.frame(G1 = as.character(df1$G1), G2 = as.character(df1$G2), data = l$residuals);
    }
    
    # add the modelto the models list
    models = c(models,list(l));
    names(models)[i] = all_genes[i];
    scores = rbind(scores, df2);
    intersects = rbind(intersects, data.frame(gene = all_genes[i], data = intersect));
  }
  scores$data = matrix(scores$data, nrow = length(scores$data), ncol = 1);
  return(list(scores = scores, models = models, singles = singles, dataN = dataN, intersects = intersects));
}

# scatter singles vs interaction degree
scatter_singles_vs_int_degree = function(data, qL = 0.05, qH = 0.95, tL = NULL, tH = NULL, with_plot = TRUE, d = NULL, point_size = 2, line_thickness = 2, average_scores = F, plot_title = "") {
  # pull the relevant data from the data structure
  singles = data$singles;
  scores = data$scores;
  
  # average scores?
  if (average_scores == T) {
    scores = average_scores(scores);
  }
  
  # calculate the thresholds for significant interactions if not defined in the function call
  if (is.null(tL)) {
    tL = quantile(as.matrix(scores$data), qL, na.rm = T);
  }
  if (is.null(tH)) {
    tH = quantile(as.matrix(scores$data),qH, na.rm = T);
  }
  
  # get the list of genes
  all_genes = union(scores$G1, scores$G2);
  
  # initialize matrix
  degrees = matrix(NA, nrow = length(all_genes), ncol = 1);
  
  # compute interaction degree for each gene and populate matrix
  for (i in 1:length(all_genes)) {
    pairs = extract_pairs(all_genes[i],scores);
    num_sig = length(which(pairs$data < tL | pairs$data > tH));
    num_valid = length(which(!is.na(pairs$data))) + length(which(!is.na(pairs$data)));
    degrees[i,1] = num_sig/num_valid;
  }
  dimnames(degrees) = list(all_genes);
  
  # reorder the degrees and singles data so they r co-linear
  i = intersect_w_indeces(singles$gene, rownames(degrees));
  
  # create data fra,me for plotting
  df = data.frame(gene = as.factor(i$i), singles = singles$data[i$ia], degree = degrees[i$ib,1]);
  
  require(ggplot2);
  d = ggplot(df,aes(x = singles, y = degree)) + geom_point(size = point_size) + stat_smooth(colour = "red", method="loess", formula = y~x, size = line_thickness) + ggtitle(sprintf("%s thresholds: %.2f %.2f",plot_title,tL,tH)) + xlab("Singles phenotype") + ylab("Interaction degree");
  # plot?
  if (with_plot) {
    print(d);
  }
  return(list("df" = df, "plot" = d));
}