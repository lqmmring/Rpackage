#' @title gene_set_selection
#'
#' @description Construct multiple gene subsets with some top HVGs
#' @param sc.counts Row (genes) by column (cells) log-normalized expression matrix
#' @param set_num_list A list provided the number of selected top HVGs for each gene subset
#' @return Return a list contains logNorm expression matrix and geneSet.
#' @export
#'
#'

gene_set_selection<-function(sc.counts,
                             set_num_list =NULL
){
  if (is.null(set_num_list)) {
    set_num_list<-c(50,100,250,500,800)
  }
  #
  library(Seurat)
  counts.list<-list()
  results.list<-list()
  num.counts<-length(set_num_list)
  expression_counts <- Seurat::CreateSeuratObject(counts = sc.counts)
  expression_counts <- Seurat::NormalizeData(expression_counts,
                                             normalization.method = 'LogNormalize',
                                             scale.factor = 10000)
  expression_counts <- Seurat::FindVariableFeatures(expression_counts, nfeature= 1000)
  top.genes<-expression_counts@assays[["RNA"]]@var.features

  counts.list<-parallel::mclapply(
    X = 1:num.counts,
    mc.preschedule = FALSE,
    mc.cores = num.counts,
    FUN = function(round) {
      # f<-as(sc.counts[top.genes[1:set_num_list[round]],], 'dgCMatrix')
      f<-sc.counts[top.genes[1:set_num_list[round]],]
      return(f)})

  results.list[['logNorm']]<-expression_counts@assays[["RNA"]]@data
  results.list[['geneSets']]<-counts.list

  return(results.list)
}


#' @title similarity_multi_gene
#'
#' @description Construct multiple similarity matrices with KNN graphs of multiple gene subsets
#' @param counts.list A list geneSet contains gene's name in each gene subset
#' @param neighborNum Number of nearest neighboors
#' @param cores Number of cores
#' @return Return a list contains multiple similarity matrices.
#' @export
#'
#'

similarity_multi_gene<-function(counts.list,
                                neighborNum,
                                cores = 5 # equal to the number of elemets in list
){

  num.counts<-length(counts.list)
  neighbor.list<-list()

  neighbor.list<-parallel::mclapply(
    X = 1:num.counts,
    mc.preschedule = FALSE,
    mc.cores = cores,
    FUN = function(round) {

      neighbor<-list()
      nn_network = dbscan::kNN(x = Matrix::t(counts.list[[round]]), k = neighborNum, sort = TRUE)
      neighbor_index = as.matrix(as.data.frame(nn_network['id']))
      row.names(neighbor_index)=colnames(counts.list[[round]])
      neighbor_dist = as.matrix(as.data.frame(nn_network['dist']))
      row.names(neighbor_dist)=colnames(counts.list[[round]])
      neighbor_weight = neighbor_dist/rowSums(neighbor_dist)
      neighbor[['index']]<-neighbor_index
      neighbor[['weight']]<-neighbor_weight
      return(neighbor)})

  return(neighbor.list)
}


#' @title regression_multi_gene
#'
#' @description Construct multiple similarity matrices with regression of multiple gene subsets
#' @param counts.list A list geneSet contains gene's name in each gene subset
#' @param neighborNum Number of nearest neighboors
#' @param num.sample Number of samples in regression
#' @param cores Number of cores
#' @return Return a list contains multiple similarity matrices.
#' @export
#'
#'

regression_multi_gene<-function(counts.list, #each element of list is a matrix: gene*cell
                                neighborNum,
                                num.sample,
                                cores = 5 # equal to the number of elemets in list
){

  num.counts<-length(counts.list)
  neighbor.list<-list()

  neighbor.list<-parallel::mclapply(
    X = 1:num.counts,
    mc.preschedule = FALSE,
    mc.cores = cores,
    FUN = function(round) {

      neighbor<-list()
      temp.regression<-c()
      num.gene.inset=dim(counts.list[[round]])[1]
      nn_network = dbscan::kNN(x = Matrix::t(counts.list[[round]]), k = neighborNum, sort = TRUE)
      neighbor_index = as.matrix(as.data.frame(nn_network['id']))
      row.names(neighbor_index)=colnames(counts.list[[round]])
      # neighbor_dist = as.matrix(as.data.frame(nn_network['dist']))
      # row.names(neighbor_dist)=colnames(counts.list[[round]])
      # neighbor_weight = neighbor_dist/rowSums(neighbor_dist)

      # for (i in c(1:num.sample)) {
      #   # cv<-glmnet::cv.glmnet(counts.list[[round]][,neighbor_index[i,]],
      #   #                       counts.list[[round]][,i],
      #   #                       family = 'poisson', dfmax = neighborNum)
      #   # cc<-coef(cv, s ="lambda.min")
      #   cv<-glmnet::glmnet(counts.list[[round]][,neighbor_index[i,]],
      #                      counts.list[[round]][,i],
      #                      lambda = 0.1,family = 'poisson', dfmax = neighborNum)
      #   cc<-coef(cv)
      #   temp.cc<-as(matrix(cc[-1,1], nrow = 1),'dgCMatrix')
      #   temp.regression<-rbind(temp.regression, temp.cc)
      # }
      M = counts.list[[round]]
      x<-matrix(c(1:num.sample),nrow =1)
      # temp.regression<-as(matrix(0, nrow = num.sample, ncol = num.gene.inset),'dgCMatrix')
      regression.matrix<-as(matrix(0, nrow = num.sample, ncol = num.gene.inset),'dgCMatrix')
      regression.matrix<-apply(
        X = x, MARGIN = 2,
        FUN = function(i) {
          cv<-glmnet::glmnet(M[,neighbor_index[i,]],
                             M[,i],
                             lambda = 0.1,family = 'poisson',dfmax = neighborNum)
          cc<-coef(cv)
          temp.cc<-as(matrix(cc[-1,1], nrow = 1),'dgCMatrix')
          # temp.regression[i,]<-temp.cc
          # return(temp.regression)
          return(temp.cc)
        })
      regression.matrix<-do.call(rbind, regression.matrix)
      neighbor[['index']]<-neighbor_index
      # neighbor[['weight']]<-neighbor_weight
      # neighbor[['regression']]<-temp.regression
      neighbor[['regression']]<-regression.matrix
      return(neighbor)})

  return(neighbor.list)
}


#' @title initialize_variables_multi_gene
#'
#' @description Population initialization with multiple similarity matrices
#' @param data (genes) by column (cells) log-normalized expression matrix
#' @param counts.list A list geneSet contains gene's name in each gene subset
#' @param num.pop Number of populations
#' @param num.M Number of multiobjective
#' @param num.sample Number of cells
#' @param min_range lower range of each variable
#' @param max_range upper range of each variable
#' @param knn_id the k neighbors of each cell
#' @param paralle Whether implement parallel operation, 'FALSE'(default)
#' @param cores Number of cores, cores=1(default)
#' @param neighbor.methdod Method of construct similarity matrix, 'similarity'(default)
#' @return Return a list contains multiple similarity matrices.
#' @export
#'
#'

initialize_variables_multi_gene<-function (data,
                                           counts.list,
                                           num.pop,
                                           num.M,
                                           num.sample,
                                           min_range,
                                           max_range,
                                           knn.list,
                                           paralle=FALSE,
                                           cores=1,
                                           neighbor.method = 'similarity'){
  # num.M:number of multiobjective
  # num.sample: numbe of cells
  # min_range, max_range: lower and upper range of each variable
  # knn_id: the k neighboors of each cell
  # knn_weight: the weight of each edge

  # return: a list contains sparse popluation matrixs and objective fitness

  # data<-small_counts_scaled
  # num.pop<-40
  # num.genes<-geneNum
  # num.sample<-sampleNum
  # knn_id<-neighbor_index
  # knn_weight<-neighbor_weight
  # cores<-1
  # paralle<-TRUE
  #
  cat('estimate weight between neighbor cells using',neighbor.method,'.\n')
  num.counts<-length(counts.list)
  n<-1/num.counts
  variables<-list()
  ## initialize sparse pop matrix
  if (neighbor.method == 'similarity') {
    if (isFALSE(paralle)){
      sparse_matrix.list<-list()
      for (i in c(1:num.pop)){
        f<-matrix(0,nrow = num.sample, ncol = num.sample)
        de<-runif(1, min = 0, max = 1)
        m<-ceiling(de/n)
        knn_id=knn.list[[m]][['index']]
        knn_weight=knn.list[[m]][['weight']]

        for (j in c(1:num.sample)) {
          # f[knn_id[j,1:i],j]=knn_weight[j,1:i]
          f[knn_id[j,],j]=knn_weight[j,]
        }
        f<-SumoneNormalization(f)
        row.names(f)<-colnames(data)
        colnames(f)<-colnames(data)
        f<-as(f, 'dgCMatrix')
        sparse_matrix.list[[i]]<-f
      }

      ## calculate evaluate objective
      obj_value.list<-list()
      for (i in c(1:num.pop)) {
        value<-evaluatte_objective(data, sparse_matrix.list[[i]])
        colnames(value)<-c('obj1','obj2')
        obj_value.list[[i]]<-value
      }
      # cat('the initialization costs time:', (proc.time() -t1)[3],'\n')

    }else{## initialize sparse pop matrix with multi cores parallel
      sparse_matrix.list<-list()
      for (i in c(1:num.pop)){
        f<-matrix(0,nrow = num.sample, ncol = num.sample)
        de<-runif(1, min = 0, max = 1)
        m<-ceiling(de/n)
        knn_id=knn.list[[m]][['index']]
        knn_weight=knn.list[[m]][['weight']]

        for (j in c(1:num.sample)) {
          # f[knn_id[j,1:i],j]=knn_weight[j,1:i]
          f[knn_id[j,],j]=knn_weight[j,]
        }
        f<-SumoneNormalization(f)
        row.names(f)<-colnames(data)
        colnames(f)<-colnames(data)
        f<-as(f, 'dgCMatrix')
        sparse_matrix.list[[i]]<-f
      }

      # explain -----------------------------------------------------------------
      # parallel::mclapply function here is implemened slowly
      # sparse_matrix.list<-list()
      # sparse_matrix.list<-parallel::mclapply(
      #   X = 1:num.pop,
      #   mc.preschedule = FALSE,
      #   mc.cores = cores,
      #   FUN = function(round) {
      #
      #     f<-matrix(0,nrow = num.sample, ncol = num.sample)
      #     de<-runif(1, min = 0, max = 1)
      #     m<-ceiling(de/n)
      #     knn_id=knn.list[[m]][['index']]
      #     knn_weight=knn.list[[m]][['weight']]
      #     f[knn_id,]=knn_weight
      #     for (j in c(1:num.sample)) {
      #       f[knn_id[j,1:round],j]=knn_weight[j,1:round]
      #       # f[knn_id[j,],j]=knn_weight[j,]
      #     }
      #     # f<-SumoneNormalization(f)
      #     row.names(f)<-colnames(data)
      #     colnames(f)<-colnames(data)
      #     f<-as(f, 'dgCMatrix')
      #     return(f)
      #   })

      obj_value.list<-list()
      obj_value.list<-parallel::mclapply(
        X = 1:num.pop,
        mc.preschedule = FALSE,
        mc.cores = cores,
        FUN = function(round) {
          value<-evaluatte_objective(data, sparse_matrix.list[[round]])
          colnames(value)<-c('obj1','obj2')
          return(value)})
    }

  }else {
    if (isFALSE(paralle)){
      sparse_matrix.list<-list()
      for (i in c(1:num.pop)){
        f<-matrix(0,nrow = num.sample, ncol = num.sample)
        de<-runif(1, min = 0, max = 1)
        m<-ceiling(de/n)
        knn_id=knn.list[[m]][['index']]
        knn_weight=knn.list[[m]][['regression']]

        for (j in c(1:num.sample)) {
          # f[knn_id[j,1:i],j]=knn_weight[j,1:i]
          f[knn_id[j,],j]=knn_weight[j,]
          while(sum(f[,j]==0)==num.sample){
            q= ceiling(runif(1,min = 0, max = 1)*num.sample)
            f[knn_id[q,],j]=knn_weight[q,]
          }
        }
        f<-SumoneNormalization(f)
        row.names(f)<-colnames(data)
        colnames(f)<-colnames(data)
        f<-as(f, 'dgCMatrix')
        sparse_matrix.list[[i]]<-f
      }

      ## calculate evaluate objective
      obj_value.list<-list()
      for (i in c(1:num.pop)) {
        value<-evaluatte_objective(data, sparse_matrix.list[[i]])
        colnames(value)<-c('obj1','obj2')
        obj_value.list[[i]]<-value
      }

    }else{## initialize sparse pop matrix with multi cores parallel

      sparse_matrix.list<-parallel::mclapply(
        X = 1:num.pop,
        mc.preschedule = FALSE,
        mc.cores = cores,
        FUN = function(round) {

          f<-matrix(0,nrow = num.sample, ncol = num.sample)
          de<-runif(1, min = 0, max = 1)
          m<-ceiling(de/n)
          knn_id=knn.list[[m]][['index']]
          knn_weight=knn.list[[m]][['regression']]
          for (j in c(1:num.sample)) {
            # f[knn_id[j,1:round],j]=knn_weight[j,1:round]
            f[knn_id[j,],j]=knn_weight[j,]
            while(sum(f[,j]==0)==num.sample){
              q= ceiling(runif(1,min = 0, max = 1)*num.sample)
              f[knn_id[q,],j]=knn_weight[q,]
            }
          }
          f<-SumoneNormalization(f)
          row.names(f)<-colnames(data)
          colnames(f)<-colnames(data)
          f<-as(f, 'dgCMatrix')

          return(f)

        })

      obj_value.list<-parallel::mclapply(
        X = 1:num.pop,
        mc.preschedule = FALSE,
        mc.cores = cores,
        FUN = function(round) {

          value<-evaluatte_objective(data, sparse_matrix.list[[round]])
          colnames(value)<-c('obj1','obj2')

          return(value)

        })

    }

  }


  variables[['population']]<-sparse_matrix.list
  variables[['objective']]<-obj_value.list

  return(variables)
}


#' @title evaluatte_objective
#'
#' @description Evaluate objective function of each population
#' @param data (genes) by column (cells) log-normalized expression matrix
#' @param sparse.matrix sparse matrix
#' @return Return a matrix contains values of objective functions.
#' @export
#'
#'

evaluatte_objective<-function(data, sparse.matrix){
  predict.data<-data%*%sparse.matrix
  # library(Metrics)
  # obj1<-sum(colSums((data-predict.data)^2))
  if (is.numeric(sparse.matrix)&is.numeric(data)) {
    obj1<-Metrics::mse(data, predict.data)
  }else{
    obj1<-Metrics::mse(as.numeric(data), as.numeric(predict.data))
  }
  if(is.na(obj1)){obj1=10000}
  # obj1<-Metrics::mse(data, predict.data)
  obj2<-sum(colSums(as.matrix(sparse.matrix)!=0))
  if(is.na(obj2)){obj2=1e20}
  res<-cbind(obj1,obj2)
  return(res)
}

#' @title SumoneNormalization
#'
#' @description sum normalization
#' @param data data matrix
#' @return Return a sum normalized matrix
#' @export
#'
#'

SumoneNormalization<-function(data){
  data.norm<-apply(data, 2, function(x) {
    sum.x<-sum((x))
    norm.x<-x/sum.x
    return(norm.x)
  })
}

#' @title non_domination_sort_mod
#'
#' @param variables a list contains population matrices and objective values vector
#' @param num.M number of multi-objective
#' @param num.sample number of samples
#' @param num.pop number of population
#' @return Return a a list contains population matrices and objective values vector
#' @export
#'
#'

non_domination_sort_mod<-function(variables,
                                  num.M,
                                  num.sample,
                                  num.pop){
  # variables: a list contains two elements:
  # (1) pop_size sparse population matrices (num.sample * num.sample)
  # (2) pop_size objective values vector (1*num.M)
  # return: a list contains two elements:
  # (1) pop_size sparse population matrices (num.sample * num.sample)
  # (2) pop_size objective values vector:obj1, obj2, front and distance (1*(num.M+2))

  ##### algorithm starts
  front=1
  individal.p<-list()
  individal.n<-list()
  P.f<-list()
  P.f[[front]]<-c(0)
  obj_value.list<-variables[['objective']]
  pop<-matrix(as.matrix(as.data.frame(obj_value.list)),ncol = num.M)
  pop<-cbind(pop, matrix(0, nrow = num.pop, ncol = 2))
  colnames(pop)<-c('obj1','obj2','front','dist')

  for (i in c(1:num.pop)){
    individal.p[[i]]<-c(0)
    individal.n[[i]]=0
    # print(P.f[[front]])

    for (j in c(1:num.pop)) {
      dom_less=0
      dom_equal=0
      dom_more=0

      for (k in c(1:num.M)) {
        if (pop[i,k]<pop[j,k]) {dom_less=dom_less+1}
        else if (pop[i,k]==pop[j,k]) {dom_equal=dom_equal+1}
        else {dom_more=dom_more+1}
      }
      if (dom_less==0 & dom_equal != num.M)
      {
        # print(paste(i,'worser than',j))
        individal.n[[i]]=individal.n[[i]]+1}
      else if (dom_more==0 & dom_equal != num.M)
      {individal.p[[i]]<-append(individal.p[[i]],j)}

      if (individal.p[[i]][1]==0&length(individal.p[[i]])>1)
      {
        # print(paste('run_delete_0',i))
        individal.p[[i]]=individal.p[[i]][-1]}
    }
    if (individal.n[[i]]==0) {
      # print(paste('add',i))
      pop[i,(num.M+1)]=front
      # print(paste('before',P.f[[front]]))
      P.f[[front]]<-append(P.f[[front]],i)
      # print(paste('after',P.f[[front]]))
      # print(paste('run_add_front',i))
    }
    if (P.f[[front]][1]==0&length(P.f[[front]])>1){
      P.f[[front]]<-P.f[[front]][-1]
      # print(paste('after_delete_0',P.f[[front]]))
    }
  }

  while (P.f[[front]][1]!=0) {
    Q<-c(0)
    for (i in c(1:length(P.f[[front]]))) {

      if (individal.p[[P.f[[front]][i]]][1]!=0) {

        for (j in c(1:length(individal.p[[P.f[[front]][i]]]))) {

          individal.n[individal.p[[P.f[[front]][i]]][j]]=as.numeric(individal.n[individal.p[[P.f[[front]][i]]][j]])-1

          if (individal.n[individal.p[[P.f[[front]][i]]][j]]==0) {

            pop[individal.p[[P.f[[front]][i]]][j], num.M+1]=front+1
            Q<-append(Q,individal.p[[P.f[[front]][i]]][j])
          }
        }
      }
    }
    front=front+1
    P.f[[front]]=Q
    if (P.f[[front]][1]==0&length(P.f[[front]])>1)
    {P.f[[front]]<-P.f[[front]][-1]}
  }

  index_of_fronts<-order(pop[,num.M+1])
  sorted_based_on_front<-matrix(0,nrow = length(index_of_fronts),ncol = (num.M*2+2))

  for (i in c(1:length(index_of_fronts))){
    sorted_based_on_front[i,1:(num.M+2)]= pop[index_of_fronts[i],]}

  current_index=0
  next_pop<-matrix(0, nrow = num.pop, ncol = (num.M+2))

  # crowding distance
  y<-matrix(0,nrow = num.pop, ncol = (num.M*2+2))

  for (front in c(1:(length(P.f)-1))) {
    # distance=0
    # y<-matrix(0,nrow = num.pop, ncol = (num.M*2+2))
    previous_index = current_index + 1

    for (i in c(1:length(P.f[[front]]))) {
      y[P.f[[front]][i],]=sorted_based_on_front[current_index+i,]
      # print(paste('front:',front))
      # print(paste('i:',i))
      # print(paste('current+i:',current_index+i))
    }
    current_index=current_index+i


    if (length(P.f[[front]])==1) {
      y[P.f[[front]][1],(num.M+3):(num.M*2+2)]<-c(Inf,Inf)
    }
    else {
      for (i in c(1:num.M)) {
        index_of_objectives<-order(y[,i])
        sorted_based_on_objective<-matrix(0,nrow = length(index_of_objectives), ncol = (num.M*2+2))
        for (j in c(1:length(index_of_objectives))) {
          sorted_based_on_objective[j,]=y[index_of_objectives[j],]}

        f.max=sorted_based_on_objective[length(index_of_objectives),i]
        f.min=sorted_based_on_objective[1,i]
        y[index_of_objectives[length(index_of_objectives)],(num.M+2+i)]<-Inf
        y[index_of_objectives[1],(num.M+2+i)]<-Inf

        for (j in c(2:(length(index_of_objectives)-1))) {
          # print(paste('num 0f index_obj',j))
          # print(paste('num 0f obj',i))

          next_obj= sorted_based_on_objective[j+1,i]
          previous_obj=sorted_based_on_objective[j-1,i]

          if (f.max-f.min==0) {y[index_of_objectives[j],(num.M+2+i)]<-Inf}
          else {
            y[index_of_objectives[j],(num.M+2+i)]=(next_obj-previous_obj)/(f.max-f.min)}
        }
      }
    }

    distance<-matrix(0,nrow = num.pop,ncol = 1)
    if (length(P.f[[front]])==1) {
      distance[P.f[[front]],1]<-sum(y[P.f[[front]],(num.M+3):(2*num.M+2)])}
    else {
      distance[P.f[[front]],1]<-rowSums(y[P.f[[front]],(num.M+3):(2*num.M+2)])}

    y[P.f[[front]],num.M+2]=distance[P.f[[front]],1]
    next_pop[previous_index:current_index,]<- y[P.f[[front]],1:(num.M+2)]
  }
  temp.sparse_matrix<-list()
  temp.obj_value<-list()
  for (m in c(1:num.pop)) {
    temp.sparse_matrix[[m]]<-variables[['population']][[index_of_fronts[m]]]
    temp.names<-colnames(variables[['population']][[index_of_fronts[m]]])
    colnames(temp.sparse_matrix[[m]])<-temp.names
    row.names(temp.sparse_matrix[[m]])<-temp.names
    temp.obj_value[[m]]<-matrix(next_pop[m,],ncol = num.M+2)
  }
  sorted_variavles<-list()
  sorted_variavles[['population']]<-temp.sparse_matrix
  sorted_variavles[['objective+front+dist']]<-temp.obj_value

  return(sorted_variavles)
}


#' @title genetic_operator_with_multi_gene
#'
#' @param parent_variables a list contains parent population matrices and objective values vector
#' @param data expressed matrix
#' @param num.M Number of multi-objective
#' @param num.genes Number of genes
#' @param num.sample Number of samples
#' @param num.neighbor
#' @param delta Parameter delta
#' @param echo Parameter echo
#' @param l_limit Lower limitation of each genetic_operator
#' @param u_limit Upper limitation of each genetic_operator
#' @param crossover.p Probability perform crossover
#' @param knn_id the k neighbors of each cell
#' @param paralle Whether implement parallel operation, 'FALSE'(default)
#' @param cores Number of cores, cores=1(default)
#' @param neighbor.methdod Method of construct similarity matrix, 'similarity'(default)
#'
#' @return Return a pop_size by (genes+multiobjective) populaition matrix
#' @export
#'
#'

genetic_operator_with_multi_gene<-function(parent_variables,
                                           data,
                                           num.M,
                                           num.genes,
                                           num.sample,
                                           num.neighbor,
                                           delta, echo,
                                           l_limit,
                                           u_limit,
                                           crossover.p,
                                           knn.list,
                                           paralle=FALSE,
                                           cores = 1,
                                           neighbor.method = 'similarity'){

  # parent_variables:a list contains three elements:
  # (1) pop_size sparse populaition matries (num.sample * num.sample)
  # (2) pop_size objective values vector (1*num.M)
  # (3) pop_size operator crossover number and mutation number

  ##### algorithm starts
  num.pop<-length(parent_variables[[2]])
  offspring_variables<-list()

  ## dgCMatrix is slowly.
  #sparse_matrix.list<-parent_variables[[1]]

  temp.sparse_matrix.list<-parent_variables[[1]]
  if (!is.matrix(temp.sparse_matrix.list[[1]])) {
    sparse_matrix.list<-lapply(X = temp.sparse_matrix.list, FUN = function(X){
      f=as.matrix(X)
      return(f)
    })}

  obj_value.list<-parent_variables[[2]]
  num.counts<-length(knn.list)
  n<-1/num.counts

  if (neighbor.method == 'similarity') {

    if (isFALSE(paralle)){

      for (i in c(1:num.pop)) {
        # With crossover.p(default:90) % probability perform crossover
        # t1<-proc.time()
        sparse_matrix<-sparse_matrix.list[[i]]
        # sparse_matrix<-as.matrix(sparse_matrix.list[[i]])
        # cat('the GO costs time:', (proc.time() -t1)[3],'\n')

        if (runif(1,min = 0, max = 1)<crossover.p){
          # num.c1<-num.c1+1
          if (runif(1,min = 0, max = 1)<=0.5){
            # case1
            # t1<-proc.time()
            for (j in c(1:num.sample)){
              index.parent.1=round(num.pop*runif(1,min = 0,max = 1))
              if (index.parent.1<1) { index.parent.1=1}
              sparse_matrix[,j]=sparse_matrix.list[[index.parent.1]][,j]
              # sparse_matrix[,j]=as.matrix(sparse_matrix.list[[index.parent.1]])[,j]
            }
            # cat('the GO costs time:', (proc.time() -t1)[3],'\n')
          }
          else {
            # case 2
            index.parent.2=round(num.pop*runif(1,min = 0,max = 1))
            if (index.parent.2<1) { index.parent.2=1}
            sparse_matrix=delta*sparse_matrix+(1-delta)*sparse_matrix.list[[index.parent.2]]
            # 这里没有判断上下边界是否溢出，不知是否合理？
            sparse_matrix=SumoneNormalization(sparse_matrix)
          }

          obj_value.list[[i]]<-evaluatte_objective(data, sparse_matrix)
          colnames(obj_value.list[[i]])<-c('obj1','obj2')
          row.names(sparse_matrix)<-colnames(data)
          colnames(sparse_matrix)<-colnames(data)
          sparse_matrix.list[[i]]<-sparse_matrix

          was.crossover<-TRUE
          was.mutation<-FALSE
        }
        else {
          # With 10 % probability perform mutation. Mutation is based on
          # polynomial mutation.
          # num.c2<-num.c2+1

          de<-runif(1, min = 0, max = 1)
          m<-ceiling(de/n)
          knn_id=knn.list[[m]][['index']]
          knn_weight=knn.list[[m]][['weight']]

          for (j in c(1:num.sample)){
            for (k in c(1:(num.neighbor))){
              k_index<-knn_id[j,k]
              if (sparse_matrix[k_index,j]!=0) {
                if (runif(1,min = 0, max = 1)<=(k/num.neighbor)) {
                  sparse_matrix[k_index,j]=0
                }
                else {
                  sparse_matrix[k_index,j]=echo*sparse_matrix[k_index,j]
                }
              }
              else {
                if (runif(1,min = 0, max = 1)<=(1-(k/num.neighbor))) {
                  sparse_matrix[k_index,j]=knn_weight[j,k]
                }
                else {
                  sparse_matrix[k_index,j]=0
                }
              }
            }
          }
          row.names(sparse_matrix)<-colnames(data)
          colnames(sparse_matrix)<-colnames(data)
          obj_value.list[[i]]<-evaluatte_objective(data, sparse_matrix)
          sparse_matrix.list[[i]]<-sparse_matrix
          colnames(obj_value.list[[i]])<-c('obj1','obj2')
        }

        offspring_sparse_matrix.list<-sparse_matrix.list
        offspring_obj_value.list<-obj_value.list
      }

    }else{

      offspring_sparse_matrix.list<-parallel::mclapply(

        X = 1:num.pop,
        mc.preschedule = FALSE,
        mc.cores = cores,
        FUN = function(round){

          sparse_matrix<-sparse_matrix.list[[round]]

          if (runif(1,min = 0, max = 1)<crossover.p){
            # num.c1<-num.c1+1
            if (runif(1,min = 0, max = 1)<=0.5){
              # case1
              for (j in c(1:num.sample)){
                index.parent.1=round(num.pop*runif(1,min = 0,max = 1))
                if (index.parent.1<1) { index.parent.1=1}
                sparse_matrix[,j]=sparse_matrix.list[[index.parent.1]][,j]
              }
            }
            else {
              # case 2
              index.parent.2=round(num.pop*runif(1,min = 0,max = 1))
              if (index.parent.2<1) { index.parent.2=1}
              sparse_matrix=delta*sparse_matrix+(1-delta)*sparse_matrix.list[[index.parent.2]]
              # 这里没有判断上下边界是否溢出，不知是否合理？
              sparse_matrix=SumoneNormalization(sparse_matrix)
            }

            was.crossover<-TRUE
            was.mutation<-FALSE
          }
          else {
            # With 10 % probability perform mutation. Mutation is based on
            # polynomial mutation.
            # num.c2<-num.c2+1
            de<-runif(1, min = 0, max = 1)
            m<-ceiling(de/n)
            knn_id=knn.list[[m]][['index']]
            knn_weight=knn.list[[m]][['weight']]

            for (j in c(1:num.sample)){
              for (k in c(1:(num.neighbor))){
                k_index<-knn_id[j,k]
                if (sparse_matrix[k_index,j]!=0) {
                  if (runif(1,min = 0, max = 1)<=(k/num.neighbor)) {
                    sparse_matrix[k_index,j]=0
                  }
                  else {
                    sparse_matrix[k_index,j]=echo*sparse_matrix[k_index,j]
                  }
                }
                else {
                  if (runif(1,min = 0, max = 1)<=(1-(k/num.neighbor))) {
                    sparse_matrix[k_index,j]=knn_weight[j,k]
                  }
                  else {
                    sparse_matrix[k_index,j]=0
                  }
                }
              }
            }
            row.names(sparse_matrix)<-colnames(data)
            colnames(sparse_matrix)<-colnames(data)
          }
          return(sparse_matrix)
        })

      offspring_obj_value.list<-parallel::mclapply(
        X = 1:num.pop,
        mc.preschedule = FALSE,
        mc.cores = cores,
        FUN = function(round) {

          value<-evaluatte_objective(data, offspring_sparse_matrix.list[[round]])
          colnames(value)<-c('obj1','obj2')

          return(value)

        })

    }
  }else{
    if (isFALSE(paralle)){

      for (i in c(1:num.pop)) {
        # With crossover.p(default:90) % probability perform crossover
        sparse_matrix<-sparse_matrix.list[[i]]

        if (runif(1,min = 0, max = 1)<crossover.p){
          # num.c1<-num.c1+1
          if (runif(1,min = 0, max = 1)<=0.5){
            # case1
            for (j in c(1:num.sample)){
              index.parent.1=round(num.pop*runif(1,min = 0,max = 1))
              if (index.parent.1<1) { index.parent.1=1}
              sparse_matrix[,j]=sparse_matrix.list[[index.parent.1]][,j]
            }
          }
          else {
            # case 2
            index.parent.2=round(num.pop*runif(1,min = 0,max = 1))
            if (index.parent.2<1) { index.parent.2=1}
            sparse_matrix=delta*sparse_matrix+(1-delta)*sparse_matrix.list[[index.parent.2]]
            # 这里没有判断上下边界是否溢出，不知是否合理？
            sparse_matrix=SumoneNormalization(sparse_matrix)
          }
          obj_value.list[[i]]<-evaluatte_objective(data, sparse_matrix)
          colnames(obj_value.list[[i]])<-c('obj1','obj2')
          row.names(sparse_matrix)<-colnames(data)
          colnames(sparse_matrix)<-colnames(data)
          sparse_matrix.list[[i]]<-sparse_matrix

          was.crossover<-TRUE
          was.mutation<-FALSE
        }
        else {
          # With 10 % probability perform mutation. Mutation is based on
          # polynomial mutation.
          # num.c2<-num.c2+1

          de<-runif(1, min = 0, max = 1)
          m<-ceiling(de/n)
          knn_id=knn.list[[m]][['index']]
          knn_weight=knn.list[[m]][['regression']]

          for (j in c(1:num.sample)){
            for (k in c(1:(num.neighbor))){
              k_index<-knn_id[j,k]
              if (sparse_matrix[k_index,j]!=0) {
                if (runif(1,min = 0, max = 1)<=(k/num.neighbor)) {
                  sparse_matrix[k_index,j]=0
                }
                else {
                  sparse_matrix[k_index,j]=echo*sparse_matrix[k_index,j]
                }
              }
              else {
                if (runif(1,min = 0, max = 1)<=(1-(k/num.neighbor))) {
                  sparse_matrix[k_index,j]=knn_weight[j,k]
                }
                else {
                  sparse_matrix[k_index,j]=0
                }
              }
            }
          }
          row.names(sparse_matrix)<-colnames(data)
          colnames(sparse_matrix)<-colnames(data)
          sparse_matrix.list[[i]]<-sparse_matrix
          obj_value.list[[i]]<-evaluatte_objective(data, sparse_matrix)
          colnames(obj_value.list[[i]])<-c('obj1','obj2')
        }

        offspring_sparse_matrix.list<-sparse_matrix.list
        offspring_obj_value.list<-obj_value.list
      }

    }else{

      offspring_sparse_matrix.list<-parallel::mclapply(

        X = 1:num.pop,
        mc.preschedule = FALSE,
        mc.cores = cores,
        FUN = function(round){

          sparse_matrix<-sparse_matrix.list[[round]]

          if (runif(1,min = 0, max = 1)<crossover.p){
            # num.c1<-num.c1+1
            if (runif(1,min = 0, max = 1)<=0.5){
              # case1
              for (j in c(1:num.sample)){
                index.parent.1=round(num.pop*runif(1,min = 0,max = 1))
                if (index.parent.1<1) { index.parent.1=1}
                sparse_matrix[,j]=sparse_matrix.list[[index.parent.1]][,j]
              }
            }
            else {
              # case 2
              index.parent.2=round(num.pop*runif(1,min = 0,max = 1))
              if (index.parent.2<1) { index.parent.2=1}
              sparse_matrix=delta*sparse_matrix+(1-delta)*sparse_matrix.list[[index.parent.2]]
              # 这里没有判断上下边界是否溢出，不知是否合理？
              sparse_matrix=SumoneNormalization(sparse_matrix)
            }

            was.crossover<-TRUE
            was.mutation<-FALSE
          }
          else {
            # With 10 % probability perform mutation. Mutation is based on
            # polynomial mutation.
            # num.c2<-num.c2+1
            de<-runif(1, min = 0, max = 1)
            m<-ceiling(de/n)
            knn_id=knn.list[[m]][['index']]
            knn_weight=knn.list[[m]][['regression']]

            for (j in c(1:num.sample)){
              for (k in c(1:(num.neighbor))){
                k_index<-knn_id[j,k]
                if (sparse_matrix[k_index,j]!=0) {
                  if (runif(1,min = 0, max = 1)<=(k/num.neighbor)) {
                    sparse_matrix[k_index,j]=0
                  }
                  else {
                    sparse_matrix[k_index,j]=echo*sparse_matrix[k_index,j]
                  }
                }
                else {
                  if (runif(1,min = 0, max = 1)<=(1-(k/num.neighbor))) {
                    sparse_matrix[k_index,j]=knn_weight[j,k]
                  }
                  else {
                    sparse_matrix[k_index,j]=0
                  }
                }
              }
            }
            row.names(sparse_matrix)<-colnames(data)
            colnames(sparse_matrix)<-colnames(data)
          }
          return(sparse_matrix)
        })

      offspring_obj_value.list<-parallel::mclapply(
        X = 1:num.pop,
        mc.preschedule = FALSE,
        mc.cores = cores,
        FUN = function(round) {

          value<-evaluatte_objective(data, offspring_sparse_matrix.list[[round]])
          colnames(value)<-c('obj1','obj2')

          return(value)

        })

    }

  }
  offspring_variables[['population']]<-lapply(X = offspring_sparse_matrix.list, FUN = function(X){
    if(is.matrix(X)){
      f=as(X, 'dgCMatrix')
      return(f)}else{
        f = X
        return(f)}})
  offspring_variables[['objective']]<-offspring_obj_value.list
  # offspring_variables[['crossover_mutation']]<-c(num.c1,num.c2)

  return(offspring_variables)
}


#' @title PF_select
#'
#' @description select the optimal from Pareto front
#' @param variables Pareto front population
#' @param alpha Parameter alpha
#' @return Return the index of the optimal in Pareto front.
#' @export
#'
#'


PF_select<-function(variables,alpha,only.one = FALSE){
  if (isFALSE(only.one)) {
    sparse.matrix.list<-variables[['population']]
    obj_value.list<-variables[['objective+front+dist']]
    num.pop<-length(obj_value.list)
    obj_value.list<-do.call(rbind, obj_value.list)
    index.list<-which(obj_value.list[,3]==1)

    if (length(index.list)!=1) {
      obj_value.index.list<-obj_value.list[index.list,]
      obj_value.index.list.temp<-scale(obj_value.index.list[,c(1,2)], center = TRUE, scale = TRUE)
      obj_value.index.list[,c(1,2)]<-obj_value.index.list.temp
      obj_value.sum<-alpha*obj_value.index.list[,1]+(1-alpha)*obj_value.index.list[,2]
      index<-which.min(obj_value.sum)
      min.index<-index.list[index]
    }
    else {
      min.index<-index.list
    }
  }else{
    sparse.matrix.list<-variables[['population']]
    obj_value.list<-variables[['objective']]
    num.pop<-length(obj_value.list)
    obj_value.list<-do.call(rbind, obj_value.list)
    obj_value.list.temp<-scale(obj_value.list[,c(1,2)], center = TRUE, scale = TRUE)
    obj_value.list[,c(1,2)]<-obj_value.list.temp
    obj_value.sum<-alpha*obj_value.list[,1]+(1-alpha)*obj_value.list[,2]
    min.index<-which.min(obj_value.sum)
  }
  return(min.index)

}


#' @title tournament_selection
#'
#' @description tournament selection to obtain Pareto optimals
#' @param variables  a list contains sparse population matrices and objective values vector
#' @param num.M Number of multi-objective
#' @param num.pop Number of populations
#' @param pool_size Number of populations used for pooling
#' @param tour_size Number of populations used for tour
#' @return Return a list contains sparse population matrices and objective values vector
#' @export
#'
#'

tournament_selection<-function (variables,
                                num.M,
                                num.pop,
                                pool_size,
                                tour_size){

  # variables: a list contains two elements:
  # (1) pop_size sparse population matrices (num.sample * num.sample)
  # (2) pop_size objective values vector (1*num.M+front+dist)

  ##### build population matrix
  pool_size<-round(num.pop/2)
  tour_size<-2

  ##### algorithm starts
  pop<-matrix(0,nrow = num.pop, ncol = (num.M+2))
  obj_value.list<-variables[[2]]
  # obj_value.list<-variavles[[2]]
  for (p in c(1:num.pop)) {
    pop[p,]=matrix(obj_value.list[[p]],nrow=1)
  }

  front=num.M+1
  distance=num.M+2
  select_pop<-matrix(0, nrow = pool_size, ncol = (num.M+2))
  temp.sparse_matrix<-list()
  temp.obj_value<-list()

  candidate<-c()
  for (i in c(1:pool_size)) {
    for (j in c(1:tour_size)) {
      candidate[j]=round(num.pop*runif(1,min = 0, max = 1))
      if (candidate[j]==0) {
        candidate[j]=1
      }
      if (j>1&(sum(candidate[j]==candidate[1:(j-1)])>1)) {
        candidate[j]=round(num.pop*runif(1,min = 0, max = 1))
        if (candidate[j]==0) {
          candidate[j]=1
        }
      }
    }

    c_obj_rank=c()
    c_obj_distance=c()
    for (j in c(1:tour_size)) {
      c_obj_rank[j]=pop[candidate[j],front]
      c_obj_distance[j]=pop[candidate[j],distance]
    }
    min_candidate=order(c_obj_rank)[1:sum(c_obj_rank==min(c_obj_rank))]
    if (length(min_candidate)>1) {
      max_candidate = order(c_obj_distance[min_candidate],decreasing = TRUE)[1]
      select_pop[i,]=pop[candidate[min_candidate[max_candidate]],]
      temp.sparse_matrix[[i]]=variables[['population']][[candidate[min_candidate[max_candidate]]]]
    }
    else {
      select_pop[i,]=pop[candidate[min_candidate],]
      if (length(min_candidate)>1) {
        for (m in c(1:length(min_candidate))) {
          temp.sparse_matrix[[i+m-1]]=variables[['population']][[candidate[min_candidate[m]]]]}
      }
      else {
        temp.sparse_matrix[[i]]=variables[['population']][[candidate[min_candidate]]]
      }
    }
  }
  temp.names<-colnames(variables[['population']][[1]])
  for (m in c(1:dim(select_pop)[1])) {
    colnames(temp.sparse_matrix[[m]])<-temp.names
    row.names(temp.sparse_matrix[[m]])<-temp.names
    temp.obj_value[[m]]<-matrix(select_pop[m,],ncol = num.M+2)
  }
  select_variavles<-list()
  select_variavles[['population']]<-temp.sparse_matrix
  select_variavles[['objective+front+dist']]<-temp.obj_value

  return(select_variavles)
}


#' @title replace_chromosome
#'
#' @description replace populations to obtain Pareto optimals
#' @param intermediate_variables  a list contains sparse population matrices and objective values vector
#' @param num.M Number of multi-objective
#' @param num.sample Number of samples
#' @param num.pop Number of populations
#' @return Return a list contains sparse population matrices and objective values vector
#' @export
#'
#'

replace_chromosome<-function(intermediate_variables,
                             num.M,
                             num.sample,
                             num.pop){

  # intermediate_variables:  a list contains two elements:
  # (1) pop_size+pool_size sparse populaition matries (num.sample * num.sample)
  # (2) pop_size+pool_size objective values vector (1*num.M+front+dist)
  # return: a list contains two elements: (1) pop_size sparse populaition matries (num.sample * num.sample)
  # (2) pop_size objective values vector:obj1, obj2, front and distance (1*(num.M+2))

  ##### algorithm starts
  sorted_sparse.matrix.list<-list()
  intermediate_sparse.marix.list<-intermediate_variables[[1]]
  num.intermediate_pop<-length(intermediate_sparse.marix.list)
  intermediate_pop<-matrix(0,nrow = num.intermediate_pop, ncol = (num.M+2))
  intermediate_obj_value.list<-intermediate_variables[[2]]
  replace_variavles<-list()

  # obj_value.list<-variavles[[2]]
  for (p in c(1:num.intermediate_pop)) {
    intermediate_pop[p,]=matrix(intermediate_obj_value.list[[p]],nrow=1)
  }

  index_of_front<-order(intermediate_pop[,(num.M+1)])
  sorted_pop<-intermediate_pop[index_of_front[1],]
  sorted_sparse.matrix.list[[1]]<-intermediate_sparse.marix.list[[index_of_front[1]]]
  for (i in c(2:num.intermediate_pop)) {
    sorted_pop<-rbind(sorted_pop, intermediate_pop[index_of_front[i],])
    sorted_sparse.matrix.list[[i]]<-intermediate_sparse.marix.list[[index_of_front[i]]]
  }

  max.rank<-max(intermediate_pop[, num.M+1])
  previous.index=0
  current.index=0
  next.sparse.matrix.list<-list()
  next_pop<-matrix(0,nrow = num.pop, ncol = (num.M+2))

  for (i in c(1:max.rank)) {
    add.index=length(intermediate_pop[,(num.M+1)]==i)
    current.index=current.index+add.index
    if (current.index>num.pop) {
      temp.pop<-matrix(0,nrow = add.index, ncol = (num.M+2))
      temp.sparse_matrix.list<-list()
      remain.index=num.pop-previous.index
      temp.pop=sorted_pop[previous.index+1:current.index,]
      for (m in c(1:(current.index-previous.index))) {
        temp.sparse_matrix.list[[m]]<-sorted_sparse.matrix.list[[previous.index+m]]}
      index_of_distance<-order(temp.pop[,(num.M+2)], decreasing = TRUE)

      next.pop<-temp.pop[index_of_distance[1],]
      next.sparse.matrix.list[[1]]<-temp.sparse_matrix.list[[index_of_distance[1]]]
      for (j in c(2:remain.index)) {
        next.pop<-rbind(next.pop, temp.pop[index_of_distance[j],])
        next.sparse.matrix.list[[j]]<-temp.sparse_matrix.list[[index_of_distance[j]]]
      }
      break
    }
    else if (current.index<num.pop) {
      next.pop<-sorted_pop[previous.index+1,]
      next.pop<-rbind(next.pop, sorted_pop[previous.index+2:current.index,])
      for (m in c(1:(current.index-previous.index))) {
        next.sparse.matrix.list[[m]]<-sorted_sparse.matrix.list[[previous.index+m]]}
    }
    else {
      next.pop<-sorted_pop[previous.index+1,]
      next.pop<-rbind(next.pop, sorted_pop[previous.index+2:current.index,])
      for (m in c(1:(current.index-previous.index))) {
        next.sparse.matrix.list[[m]]<-sorted_sparse.matrix.list[[previous.index+m]]}
      break
    }
    previous.index=current.index

    # # temp.names<-colnames(intermediate_variables[['population']][[1]])
    # for (m in c(1:dim(next_pop)[1])) {
    #   # colnames(next.sparse.matrix.list[[m]])<-temp.names
    #   # row.names(next.sparse.matrix.list[[m]])<-temp.names
    #   temp.obj_value[[m]]<-next.pop[m,]
    # }
    # if (i==1) {
    #   replace_variavles[[1]]<-next.sparse.matrix.list
    #   replace_variavles[[2]]<-temp.obj_value
    # }
    # else {
    #   replace_variavles[[1]]<-append(replace_variavles[[1]],next.sparse.matrix.list)
    #   replace_variavles[[2]]<-append(replace_variavles[[2]],temp.obj_value)
    #
    # }

  }

  temp.names<-colnames(intermediate_variables[['population']][[1]])
  temp.obj_value<-list()
  # temp.names<-colnames(variables[['population']][[1]])
  for (m in c(1:dim(next_pop)[1])) {
    colnames(next.sparse.matrix.list[[m]])<-temp.names
    row.names(next.sparse.matrix.list[[m]])<-temp.names
    temp.obj_value[[m]]<-matrix(next.pop[m,],ncol = num.M+2)
  }

  replace_variavles<-list()
  replace_variavles[['population']]<-next.sparse.matrix.list
  replace_variavles[['objective+front+dist']]<-temp.obj_value

  return(replace_variavles)
}


#' @title predict_matrix
#'
#' @param variables  a list contains sparse population matrices and objective values vector
#' @param data Expressed matrix
#' @param index The index of the optimal
#' @return Return the imputed expressed matrix
#' @export
#'
#'

predict_matrix<-function(variables,
                         data,
                         index){
  sparse.matrix.list<-variables[['population']]
  sparse.matrix<-do.call(rbind, sparse.matrix.list[index])
  predict.data<-data%*%sparse.matrix
  return(predict.data)
}


#' @title dataNormalization
#'
#' @param small_counts  Expressed matrix
#' @param verbos Whether to print results
#' @param normMethod The normalization method, 'lib'(default)
#' @return Return the normalization expressed matrix
#' @export
#'
#'

dataNormalization<-function(small_counts,
                            verbose,
                            normMethod = 'lib'){
  # ## STANDARDIZE 1  ##
  if(normMethod == 'norm'){
    if(verbose) message('Norm normalization')
    GEOmean <- rep(NA,rowNum)
    for (i in 1:rowNum)
    {
      gene_NZ <- small_counts[i,small_counts[i,] > 0]
      GEOmean[i] <- exp(sum(log(gene_NZ), na.rm=TRUE)/length(gene_NZ))
    }
    S <- rep(NA, colNum)
    temp_counts <- small_counts
    for (j in 1:colNum)
    {
      sample_j <- small_counts[,j]/GEOmean
      S[j] <- median(sample_j[which(sample_j != 0)])
      temp_counts[,j] <- small_counts[,j]/S[j]
    }
    temp_counts <- ceiling(temp_counts)
    small_counts <- temp_counts
  }

  # ## STANDARDIZE 2 ##
  if(normMethod == 'scaled'){
    if(verbose) message('Scaled normalization')
    temp_counts <- Matrix::t(scale(Matrix::t(small_counts), center = T, scale = T))
    small_counts <- temp_counts
  }
  # ## STANDARDIZE 3 ##
  if(normMethod == 'lib'){
    if(verbose) message('library size normalization')
    temp_counts<-colSums(small_counts)
    temp_counts[temp_counts == 0] = 1
    small_counts = sweep(small_counts, MARGIN = 2, 10^6/temp_counts, FUN = "*")
  }
  return(small_counts)
}

#' @title logTransform
#'
#' @param small_counts  Expressed matrix
#' @return Return the log-normalization expressed matrix
#' @export
#'
#'

logTransform<-function(small_counts){
  small_counts = log10(small_counts + 1.01)
  return(small_counts)
}
