# a fast implementation of searching for a matching string within a vector of strings
# returns indeces of matches
find_string = function(what,where) {
  return(which(as.character(what) == as.character(where)));
}

# a matlab style intersect also returning indeces
intersect_w_indeces = function(a,b) {
  i = intersect(a,b);
  
  if (length(i) == 0) {
    return(NULL);
  }
  
  ia = match(i,a);
  ib = match(i,b);
  ia = ia[!is.na(ia)];
  ib = ib[!is.na(ib)];
  return(list("i" = i, "ia" = ia, "ib" = ib));
}

# finds elements in common between multiple lists
get_common_elements = function(...) {
  dots = (...);
  if (length(dots) < 2) {
    return(NULL);
  }
  common = intersect(dots[[1]], dots[[2]]);
  for (i in 1:length(dots)) {
    common = intersect(common,dots[[i]])
  }
  return(common);
}

# computes medians row or columnwise in a matrix
calc_medians = function(x,mode = "by_row") {
  x = as.matrix(x);
  D = dim(x);
  
  switch(mode,
         by_row = {
           result = matrix(NA,nrow = D[1],ncol = 1);
           for (i in 1:D[1]) {
             result[i] = median(x[i,1:D[2]]);
           }
         }, by_col = {
           result = matrix(NA,nrow = 1,ncol = D[2]);
           for (i in 1:D[2]) {
             result[i] = median(x[1:D[1],i]);
           }
         });
  return(result);
}

# exports data in a format cluster3 understands
export_for_cluser = function(data, output, mode = "rect", Qs = NULL, As = NULL) {
  if (is.matrix(data)) {
    scoremat = data;
  } else {
    scoremat = reformat_to_matrix(data, mode, Qs = Qs, As = As);
  }
  
  cnames = cat(c("Gene", colnames(scoremat),"\n"), sep = "\t", file = output);
  nrows = length(rownames(scoremat));
  for (i in 1:nrows) {
    row_str = as.character(as.list(scoremat[i,]));
    row_str_na = gsub("NA","",row_str);
    cat(rownames(scoremat)[i], row_str_na,"\n",file = output, sep = "\t",append = TRUE);
  }
}

# averages scores using simple median
average_scores = function(scores, paired_data = NULL) {
  print("Averagins scores", quote = FALSE);
  if (is.null(paired_data)) {
    paired_data = get_paired_data(scores);
  }
  
  if (is.null(paired_data)) {
    print("No paired data found");
    return(NULL);
  }
  
  p = unique(paired_data$pair);
  for (i in 1:length(p)) {
    p_inds = find_string(p[i], paired_data$pair);
    av_score = median(unlist(paired_data[p_inds,c('data1','data2')]));
    #print(i);
    #av_scores = calc_medians(cbind(paired_data[,"pairs_G1"],paired_data[,"pairs_G2"]));
    scores$data[paired_data[p_inds,"ind1"]] = av_score;
    scores$data[paired_data[p_inds,"ind2"]] = av_score;
  }
  return(scores);
}

# removes queries data
remove_Qs = function(data, Qs) {
  for (i in 1:length(Qs)) {
    inds = union(find_string(Qs[i], data$Qs), find_string(Qs[i],data$Q_ID));
    if (length(inds) > 0) {
      data = data[-inds,];
    }
  }
  return(data);
}

# remove array data
remove_As = function(data, As) {
  for (i in 1:length(As)) {
    inds = find_string(As[i], data$G2);
    if (length(inds) > 0) {
      data = data[-inds,];
    }
  }
  return(data);
}

remove_data = function(data, Qs = NULL, As = NULL) {
  data = remove_Qs(data, Qs);
  data = remove_As(data, As);
  return(data);
}

# converts a score dataframe to a square or rectangular score matrix
reformat_to_matrix = function(scores, mode = "rect", Qs = NULL, As = NULL) {
  switch(mode,
         rect = {
           if (is.null(Qs)) {
             Qs = sort(unique(scores$G1));
           }
           if (is.null(As)) {
             As = sort(unique(scores$G2));
           }
         },
         square = {
           if (is.null(Qs) | is.null(As)) {
             Qs = sort(union(scores$G1, scores$G2));
             As = Qs;
           }
         });
  Qs = unique(Qs);
  As = unique(As);
  print(sprintf("Reformating to matrix of %i rows and %i columns", length(Qs), length(As)), quote = FALSE);
  result = matrix(NA,nrow = length(Qs), ncol = length(As));
  rownames(result) = Qs;
  colnames(result) = As;
  for (i in 1:dim(as.matrix(scores$data))[1]) {
    ri = find_string(scores$G1[i],Qs);
    ci = find_string(scores$G2[i],As);
    result[ri,ci] = scores$data[i];
    if (mode == 'square') {
      result[ci,ri] = scores$data[i];
    }
  }
  return(result);
}

# extracts data from synonymous gene pairs where geneA x geneB = geneB x geneA
get_paired_data = function(scores,iter=100,benchmark=F) {
  print("Extracting paired data", quote = FALSE);
  merge1 = paste(scores$G1,scores$G2);
  merge2 = paste(scores$G2,scores$G1);
  um = intersect(merge1,merge2);
  um = sort(um);
  print(length(um));
  L1 = NULL;
  
  if (benchmark==T) {
    lim = iter;
  } else {
    lim = length(um);
  }
  
  for (i in 1:lim) {
    p1 = find_string(um[i], merge1);
    p2 = find_string(um[i], merge2);
    combs = expand.grid(p1,p2);
    
    for (j in 1:dim(combs)[1]) {
      L2 = c(um[i],combs$Var1[j],scores$data[combs$Var1[j]],combs$Var2[j],scores$data[combs$Var2[j]]);
      L1 = c(L1,L2);
    }
    #print(i);
  }
  
  if (length(L1) == 0) {
    return(NULL);
  } else {
    m = matrix(L1,ncol = 5, byrow = T);
    result = data.frame(pair = m[,1], ind1 = as.numeric(m[,2]), ind2 = as.numeric(m[,4]), data1 = as.numeric(m[,3]), data2 = as.numeric(m[,5]));
    return(result);
  }
}

# extracts data from synonymous gene pairs where geneA x geneB = geneB x geneA
get_paired_data_new = function(scores,iter=100,benchmark=F) {
  print("Extracting paired data", quote = FALSE);
  merge1 = paste(scores$G1,scores$G2);
  merge2 = paste(scores$G2,scores$G1);
  um = intersect(merge1,merge2);
  um = sort(um);
  print(length(um));
  L1 = list();
  
  if (benchmark==T) {
    lim = iter;
  } else {
    lim = length(um);
  }
  
  for (i in 1:lim) {  
    p1 = find_string(um[i], merge1);
    p2 = find_string(um[i], merge2);
    combs = expand.grid(p1,p2);
    pair = um[i];
    L1[[pair]] = c(combs$Var1,combs$Var2,scores$data[combs$Var1],scores$data[combs$Var2]);
    #print(i);
  }
  
  if (length(L1) == 0) {
    return(NULL);
  } else {
    result = do.call("rbind", L1);
    return(result);
  }
}

# makes conistent datasets
make_consistent_datasets = function(...,average=F) {
  print("Creating consistent datasets", quote = FALSE);
  if (!is.list(...)) {
    dots = list(...);
  } else {
    dots = (...);
  }
  
  gene_pairs = list();
  
  # avergae all the data and generae the gene-pair strings
  for (i in 1:length(dots)) {
    if (average==T) {
      dots[[i]] = average_scores(dots[[i]]);
    }
    gene_pairs[[i]] = paste(dots[[i]]$G1, dots[[i]]$G2);
  }
  
  # get all the common gene pairs between all the datasets
  in_common = get_common_elements(gene_pairs);
  if (is.null(in_common)) {
    print ("No common pairs found", quote = FALSE);
    return(NULL);
  }
  
  # reorder all the data so its all consistent and contains the same gene pairs 
  for (i in 1:length(gene_pairs)) {
    tmp = intersect_w_indeces(in_common, gene_pairs[[i]]);
    tmp_df = data.frame(G1 = dots[[i]]$G1[tmp$ib], G2 = dots[[i]]$G2[tmp$ib], data = dots[[i]]$data[tmp$ib]);
    dots[[i]] = tmp_df;
  }
  return(dots);
}

# computes correlation coefficients between score matrice
correlate_scores = function(...,list = FALSE, names=NULL,m = "pearson") {
  if (list == FALSE) {
    dots = list(...);
  } else {
    dots = (...);
  }
  
  result = matrix(NA, nrow = length(dots), ncol = length(dots));
  for (i in 1:length(dots)) {
    for (j in 1:length(dots)) {
      if (!i == j) {
        if ("data" %in% names(dots[[i]])) {
          d1 = dots[[i]]$data;
        } else {
          d1 = dots[[i]]$scores$data;
        }
        
        if ("data" %in% names(dots[[j]])) {
          d2 = dots[[j]]$data;
        } else {
          d2 = dots[[j]]$scores$data;
        }
        
        result[i,j] = cor(d1, d2, method = m);
      }
    }
  }
  
  if (is.null(names)) {
    rownames(result) = names(dots); #c(1:length(dots));
    colnames(result) = names(dots); #c(1:length(dots));
  } else {
    rownames(result) = names;
    colnames(result) = names;
  }
  return(result);
}

# print data for control wells
get_controls_data = function(control, data, replicate_tag='replicate_') {
  c_inds = find_string(control,data$G2);
  qs = sort(unique(data$Q_ID));
  D = dim(data);
  L = list();
  for (i in 1:length(qs)) {
    q_inds = find_string(qs[i],data$Q_ID);
    q_data = unlist(data[q_inds,grep(replicate_tag,colnames(data))],use.names=F);
    q_med = median(q_data);
    
    d_inds = intersect(c_inds,q_inds);
    c_data = unlist(data[d_inds,grep(replicate_tag,colnames(data))],use.names=F);
    c_med = median(c_data);
    
    L[[qs[i]]] = c(data$G1[q_inds[1]],c_med,q_med,c_med/q_med,data$batch[q_inds[1]]);
  }
  L = do.call(rbind,L);
  result = data.frame(Q_ID = rownames(L), ctrl_med = as.numeric(L[,2]), plate_med = as.numeric(L[,3]),ratio = as.numeric(L[,4]), batch = as.numeric(L[,5]));
  return(result);
}

reformat_for_EMAP_toolbox = function(data,target_dir=test) {
  dir.create(target_dir);
  bs = unique(data$batch);
  
  for (b in 1:length(bs)) {
    b_inds = find_string(bs[b],data$batch);
    qs = unique(data$Q_ID[b_inds]);
    print(qs);
    dir_name = paste(target_dir,sprintf("batch_%i",bs[b]),sep='/');
    print(dir_name);
    dir.create(dir_name);
    screen_map = data.frame(filename=NULL, orfname = NULL, mutation = NULL);

    for (q in qs) {
      #print(q);
      q_inds = find_string(q,data$Q_ID);
      q_data = data[q_inds,];
      reps = grep('replicate_',colnames(q_data));
      for (repl in 1:length(reps)) {
        df = q_data[order(q_data$Well),c('Well',colnames(q_data)[reps[repl]])];
        df1 = data.frame(row = df$Well%/%100+1, col = df$Well%%100, size = df[,2], circ = 0.99);
        screen_map = rbind(screen_map, data.frame(filename = q, orfname = data$G1[find_string(q,data$Q_ID)[1]],mutation = 'esiRNA'));
        f_name = sprintf("%s/%s_%i.dat",dir_name,q, repl);
        print(f_name);
        write.table(df1,file = f_name,sep = "\t",row.names=F,quote=F);
      }
    }
    print(screen_map);
    write.table(screen_map,file = sprintf("%s/map.txt",dir_name),sep = "\t",row.names=F,quote=F);  
  }
}