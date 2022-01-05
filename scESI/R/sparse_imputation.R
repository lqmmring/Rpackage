#' @title sparse_imputation_with_selected_genes
#'
#' @description sparse imputation with nearest neighbors
#' @param data.path the path of the data sets
#' @param data Row (genes) by column (cells) log-normalized expression matrix
#' @param verbose Whether to print results
#' @param processing Whether to process data
#' @param evalutation Whether to evaluate the results
#' @param expressed_genes_rate_low the filtering condition of low expressed genes
#' @param expressed_genes_rate_up the filtering condition of high expressed genes
#' @param expressed_cells_rate_low the filtering condition of low expressed cells
#' @param expressed_cells_rate_up the filtering condition of high expressed cells
#' @param subset_genes Number of genes
#' @param num.pop Number of populations
#' @param num.Iteration Number of iterations
#' @param num.M Number of multi-objective
#' @param mu Parameter mu
#' @param mum Parameter mum
#' @param crossover.p Probability perform crossover
#' @param set_num_list A list provided the number of selected top HVGs for each gene subset
#' @param K Number of nearest neighbors
#' @param delta Parameter delta
#' @param echo Parameter echo
#' @param paralle Whether implement parallel operation, 'FALSE'(default)
#' @param cores Number of cores
#' @param neighbor.method  Method of construct similarity matrix, 'similarity'(default)
#'
#' @return Return a list contains logNorm expression matrix and geneSet.
#' @export
#'
#'

sparse_imputation_with_selected_genes<-function(
  data.path=NULL,
  data,
  verbose=TRUE,
  processing = TRUE,
  evalutation = FALSE,
  expressed_genes_rate_low=NULL,
  expressed_genes_rate_up=NULL,
  expressed_cells_rate_low=NULL,
  expressed_cells_rate_up=NULL,
  subset_genes=NULL,
  num.pop=30,
  num.Iteration=30,
  num.M=2,
  mu=5,
  mum=20,
  crossover.p=0.5,
  set_num_list=c(50,250,500,800,1000),
  K = 50,
  delta=0.5,
  echo=0.7,
  paralle = FALSE,
  cores =1,
  neighbor.method = 'similarity'
){
  if (is.null(data.path)) {
    counts<- as.matrix(data)
    data.sc<-c()
  }else{
    data.path<-c('/home/qmliu/projects/impuation/sparse_imputation/real_data/rds/Data_Tasic.rds')
    data.sc<-readRDS(data.path)
    counts<- data.sc@assays@data@listData[["counts"]]
  }

  if (is.null(row.names(counts))) {row.names(counts)<-paste0('genes-',c(1:dim(counts)[1]))}
  if (is.null(colnames(counts))) {colnames(counts)<-paste0('cells-',c(1:dim(counts)[2]))}

  ## processing dat a ##
  if (isTRUE(processing)) {
    message("Removing ", sum(rowSums(counts) == 0), " rows of genes with all zero counts")
    counts <- counts[rowSums(counts) != 0,]

    if(!is.null(expressed_genes_rate_low)){
      message("Removing ", sum(rowSums(counts>0)/dim(counts)[1]<expressed_genes_rate_low), " rows of genes expressed rate are too low")
      counts <- counts[rowSums(counts>0)/dim(counts)[1]>=expressed_genes_rate_low,]
    }
    if(!is.null(expressed_genes_rate_up)){
      message("Removing ", sum(rowSums(counts>0)/dim(counts)[1]>expressed_genes_rate_up), " rows of genes expressed rate are too high")
      counts <- counts[rowSums(counts>0)/dim(counts)[1]<=expressed_genes_rate_up,]
    }

    if(!is.null(expressed_cells_rate_low)){
      message("Removing ", sum(colSums(counts>0)/dim(counts)[2]<expressed_cells_rate_low), " cols of cell expressed rate are too lower")
      counts <- counts[,colSums(counts>0)/dim(counts)[2]>=expressed_cells_rate_low]
      if (isTRUE(evalutation)) {
        cell_types<- data.sc$cell_type1[which(colSums(counts>0)/dim(counts)[2]>=expressed_cells_rate_low)]
      }
    }

    if(!is.null(expressed_cells_rate_up)){
      message("Removing ", sum(colSums(counts>0)/dim(counts)[2]>expressed_cells_rate_up), " cols of cell expressed rate are too high")
      counts <- counts[,colSums(counts>0)/dim(counts)[2]<=expressed_cells_rate_up]
      if (isTRUE(evalutation)) {
        cell_types<- data.sc$cell_type1[which(colSums(counts>0)/dim(counts)[2]<=expressed_cells_rate_up)]
      }
    }
  }

  if(is.null(subset_genes)) {
    subset_genes=row.names(counts)}

  if(verbose) cat('subsetting genes \n \n')
  small_counts <- counts[subset_genes, ]

  rowNum<-dim(small_counts)[1]
  colNum<-dim(small_counts)[2]

  # gene selection and data with different gene sets ------------------------
  num.cores = length(set_num_list)
  counts.list<-gene_set_selection(sc.counts = small_counts, # matrix: gene*cell
                                  set_num_list)

  ## Memory management ##
  remove(counts, small_counts,data.sc)
  gc()

  # num.pop<-30
  # num.Iteration<- 20
  num.genes<-dim(counts.list$logNorm)[1]
  num.samples<-dim(counts.list$logNorm)[2]
  # num.M<-2
  min_range<-rep(c(0),num.samples)
  max_range<-rep(c(1),num.samples)
  # mu<-5
  # mum<-20
  # crossover.p<-0.5
  # delta<-runif(1,min = 0,max = 1)
  # echo<-runif(1,min = 0,max = 1)


  # similarity with different gene sets -------------------------------------
  if (neighbor.method =='similarity') {
    t1<-proc.time()
    neighbor.list<-similarity_multi_gene(counts.list = counts.list$geneSets, #each element of list is a matrix: gene*cell
                                         neighborNum = K,
                                         cores = num.cores # equal to the number of elemets in list
    )
    cat('the similarity construction costs time:', (proc.time() -t1)[3],'\n')
  }else{
    t1<-proc.time()
    neighbor.list<-regression_multi_gene(counts.list = counts.list$geneSets,
                                         neighborNum = K,
                                         num.sample = colNum,
                                         cores = num.cores)
    cat('the regression costs time:', (proc.time() -t1)[3],'\n')
  }

  t<-proc.time()
  cat('the optimization processing starts.\n')

  t1<-proc.time()
  chromosome = initialize_variables_multi_gene(data = counts.list$logNorm,
                                               counts.list = counts.list$geneSets,
                                               num.pop = num.pop,
                                               num.M = num.M,
                                               num.sample = num.samples,
                                               min_range = min_range,
                                               max_range = max_range,
                                               knn.list = neighbor.list,
                                               # knn.list = regression.list,
                                               paralle = paralle, cores=num.cores,
                                               # neighbor.method = 'similarity')
                                               neighbor.method = neighbor.method)

  cat('the initialization costs time:', (proc.time() -t1)[3],'\n')
  # return: a list contains   sparse popluation matrixs and objective fitness

  t1<-proc.time()
  chromosome = non_domination_sort_mod(chromosome, num.M, num.samples, num.pop)
  # non_domination_sort_mod(variables, num.M, num.sample, num.pop)
  # return: a list contains two elements: (1) pop_size sparse populaition matries (num.sample * num.sample)
  # (2) pop_size objective values vector:obj1, obj2, front and distance (1*(num.M+2))

  cat('the non_domination coste time:', (proc.time() -t1)[3],'\n')

  Iteration.chromosome.list<-list()
  Iteration.chromosome.list[['population']]<-list()
  Iteration.chromosome.list[['objective']]<-list()
  t3<-proc.time()

  for (i in c(1:num.Iteration)) {

    t2<-proc.time()
    pool = round(num.pop/2)
    tour = 2
    parent_chromosome = tournament_selection(chromosome, num.M, num.pop, pool, tour)
    # tournament_selection(variables, num.M, num.pop, pool_size, tour_size)
    # return: a list contains two elements: (1) pop_size sparse populaition matries (num.sample * num.sample)
    # (2) pop_size objective values vector:obj1, obj2, front and distance (1*(num.M+2))
    cat('the tournament selection costs time:', (proc.time() -t2)[3],'\n')

    t1<-proc.time()
    offspring_chromosome = genetic_operator_with_multi_gene(parent_variables = parent_chromosome,
                                                            data = counts.list$logNorm,
                                                            num.M = num.M,
                                                            num.genes = num.genes,
                                                            num.sample = num.samples,
                                                            num.neighbor = K,
                                                            delta = delta,
                                                            echo = echo,
                                                            l_limit = min_range,
                                                            u_limit = max_range,
                                                            crossover.p = crossover.p,
                                                            knn.list = neighbor.list,
                                                            # knn.list = regression.list,
                                                            paralle = paralle,
                                                            cores = cores,
                                                            neighbor.method = neighbor.method)
    cat('the GO costs time:', (proc.time() -t1)[3],'\n')
    # genetic_operator(parent_variables, data, num.M, num.genes,num.sample, mu, mum, delta, echo, l_limit, u_limit, crossover.p, max_iter, current_iter, knn_id, knn_weight)
    # return: parent_variables:a list contains three elements: (1) pop_size sparse populaition matries (num.sample * num.sample)
    # (2) pop_size objective values vector (1*num.M) (3) pop_size operator crossover number and mutation number

    intermediate_chromosome<-list()

    intermediate_chromosome[['population']]<-append(chromosome[['population']], offspring_chromosome[['population']])
    num.intermediate.pop<-length(intermediate_chromosome[[1]])

    temp_chromosome<-list()
    for (q in c(1:num.pop)) {
      temp_chromosome[['objective']][[q]]<-matrix(chromosome[['objective+front+dist']][[q]][1:2],nrow = 1)}

    intermediate_chromosome[['objective']]<-append(temp_chromosome[['objective']], offspring_chromosome[['objective']])


    intermediate_chromosome = non_domination_sort_mod(intermediate_chromosome, num.M, num.samples, num.intermediate.pop)
    # non_domination_sort_mod(variables, num.M, num.sample, num.pop)
    # return: a list contains two elements: (1) pop_size sparse populaition matries (num.sample * num.sample)
    # (2) pop_size objective values vector:obj1, obj2, front and distance (1*(num.M+2))

    chromosome = replace_chromosome(intermediate_chromosome, num.M, num.samples, num.pop)
    # replace_chromosome(intermediate_variables, num.M, num.sample, num.pop)
    # return: a list contains two elements: (1) pop_size sparse populaition matries (num.sample * num.sample)
    # (2) pop_size objective values vector:obj1, obj2, front and distance (1*(num.M+2))

    if ((i%%5)==0) {
      cat('the iteration is', i,'\n')
      cat('the costs time:', (proc.time() -t2)[3],'\n')
    }
    temp.select.index<-PF_select(chromosome, 0.7)

    Iteration.chromosome.list[['population']][[i]]<-chromosome[['population']][[temp.select.index]]
    Iteration.chromosome.list[['objective']][[i]]<-chromosome[['objective+front+dist']][[temp.select.index]][,c(1,2)]
  }

  cat('the  optimization costs:', (proc.time() -t3)[3],'\n')
  pf.select.index<-PF_select(Iteration.chromosome.list,0.7, only.one = TRUE)

  imputation.matrix<-predict_matrix(chromosome, counts.list$logNorm, pf.select.index)

  times<-(proc.time() -t)[3]
  cat('optimization process costs:', (proc.time() -t)[3],'\n')
  return(list(predictCount=imputation.matrix,
              IterProcess=Iteration.chromosome.list,
              Index=pf.select.index,
              Times=times))
}
