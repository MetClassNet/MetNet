#' @name combine
#'
#' @aliases combine
#'
#' @title Combine structural and statistical adjacency matrix
#'
#' @description
#' The function `combine` takes as
#' input the structural and statistical adjacency matrix, created in former
#' steps, adds them together and will report a connection between metabolites
#' in the returned when the sum exceeds the `threshold`.
#' \code{combine} returns this consensus matrix supported
#' by the structural and statistical adjacency matrices.
#' 
#'
#' @param structural list containing `numeric` structural adjacency matrix in
#' the first entry and `character` structural adjanceny matrix in the second
#' entry
#'
#' @param statistical matrix containing `numeric` statistical adjacency matrix
#'
#' @param threshold numeric, threshold value to be applied to define a
#' connection as present
#' 
#' @param model character, defines model that is used for building the consensus 
#' matrix. The model is part of `statistical`. The default is "Consensus" which creates
#' the consensus matrix by using a combination of models present in `statistical`. 
#' Besides, "pearson" and "spearman" models were implemented.
#'
#' @details The matrices will be added and a unweighted connection will
#' be reported when the value exceeds a certain value.
#'
#' @return `list`, in the first entry `matrix` of type `numeric`containing the
#' consensus adjacency matrix as described
#' above harbouring connections reported by the structual and
#' statistcal adjacency matrices. In the second entry a `matrix` of type
#' `character` the corresonding type/putative link at this position.
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' Liesa Salzer, \email{liesa.salzer@@helmholtz-muenchen.de}
#'
#' @examples
#' data("x_test", package = "MetNet")
#' x_test <- as.matrix(x_test)
#' functional_groups <- rbind(
#'     c("Monosaccharide (-H2O)", "C6H10O5", "162.0528234315"),
#'     c("Disaccharide (-H2O)", "C12H20O11", "340.1005614851"),
#'     c("Trisaccharide (-H2O)", "C18H30O15", "486.1584702945"))
#' functional_groups <- data.frame(group = functional_groups[, 1],
#'      formula = functional_groups[, 2],
#'      mass = as.numeric(functional_groups[, 3]))
#' struct_adj <- structural(x_test, functional_groups, ppm = 5)
#' stat_adj_l <- statistical(x_test,
#'     model = c("pearson", "spearman"),
#'     correlation_adjust = "bonferroni")
#' stat_adj <- threshold(stat_adj_l, type = "top2", args = list(n = 10))
#' combine(struct_adj, stat_adj)
#'
#' @export
combine <- function(structural, statistical, threshold = 1, model = "Consensus") {
  
  
  if (!is.matrix(structural[[1]]) | !is.numeric(structural[[1]]))
    stop("structural[[1]] is not a numeric matrix")
  
  
  if (!is.matrix(statistical[["Consensus"]]) | !is.numeric(statistical[["Consensus"]]))
    stop("statistical is not a numeric matrix")
  
  if (!all(rownames(structural[[1]]) == rownames(structural[[2]])))
    stop(c("rownames of structural[[1]] are not identical to rownames of ",
           "structural[[2]]"))
  
  if (!all(colnames(structural[[1]]) == colnames(structural[[2]])))
    stop(c("colnames of structural[[1]] are not identical to colnames of ",
           "structural[[2]]"))
  
  if (!all(rownames(structural[[1]]) == rownames(statistical)))
    stop("rownames are not identical")
  
  if (!all(colnames(structural[[1]]) == colnames(statistical)))
    stop("colnames are not identical")
  
  if (!is.numeric(threshold)) stop("threshold is not numeric")
  
  ## create list to store results
  res <- list()
  
  

  ## create the first entry of the list
  ## sum the matrices structural and statistical, if the value is above
  ## threshold then assign 1, otherwise 0
  cons_num <- structural[[1]] + statistical[[model]]
  cons_num <- ifelse(cons_num > threshold, 1, 0)
  
  ## create the second entry of the list
  ## if element in cons_num is equal to 1, take the element in structural[[2]]
  ## (the type of link), otherwise ""
  cons_char <- ifelse(cons_num == 1, structural[[2]], "")
  
  ## assign to list
  res[[1]] <- cons_num
  res[[2]] <- cons_char
  
  return(res)
}


#' @name adjacency_list
#'
#' @aliases adjacency_list
#'
#' @title Create a list from adjacency matrices
#'
#' @description
#' The function `adjacency_list` creates a dataframe from the adjacency matrix.
#' Unique values are displayed by the dataframe. Missing values and the upper triangle of 
#' the adjacency matrix are left out.
#' The function returns a dataframe containing different column outputs depending 
#' on the source of the adjacency matrix, which is defined by `from`.
#' 
#' @param x
#' `data.frame`, adjacency matrix, that has previously been generated 
#' by the functions `structural`, `statistical`, or `combine`
#' 
#' @param from
#' `character` defines the function from which the adjacency matrix originated.
#' This can be either "structural", "statistical", or "combine".
#'
#' @details
#' The function `adjacency_list` includes functionality to generate a dataframe from
#' adjacency matrices. The output depends on the original function from which the adjacency
#' matrix was generated. 
#'
#' @return `data.frame` containing the listed values each pair of features from 
#' adjacency matrix `x` specified by `from`.
#'
#' @author Liesa Salzer, \email{liesa.salzer@@helmholtz-muenchen.de}
#'
#' @examples
#' data("x_test", package = "MetNet")
#' x <- x_test[1:10, 3:ncol(x_test)]
#' x <- as.matrix(x)
#' stat_adj <- statistical(x = x, model = c("pearson", "spearman"))
#' adjacency_list(x = stat_adj, from = "statistical")
#'
#' @export
adjacency_list <- function(x, from){
    
    if (!(all(from %in% c("structural", "statistical", "combine"))))
        stop("'from' not implemented in adjacency_list")
    
    if ("structural" %in% from) {
        
        
        x[[2]][upper.tri(x[[2]])] <- ''
        x[[3]][upper.tri(x[[3]])] <- ''
        
        list_type <- reshape2::melt(x[[2]]) %>% filter(Var1 != Var2) %>% filter(value != '')
        list_mass <- reshape2::melt(x[[3]]) %>% filter(Var1 != Var2) %>% filter(value != '')
        combine <- merge(list_type, list_mass, by = c("Var1", "Var2"))
        colnames(combine) <- c("Var1", "Var2", "value", "mass-difference")
        return(combine)
    }
    else if ("statistical" %in% from) {
        
        for (i in seq_along(x)) {
            if (i == 1) {
                x[[i]][upper.tri(x[[i]])] <- ''
                list_corr <- reshape2::melt(x[[i]]) %>% filter(Var1 != Var2) %>% filter(value != '') %>% 
                    select(Var1, Var2, value) 
                colnames(list_corr) <- c("Var1", "Var2", names(x[i]))
            }
            if (i != 1){
                model = names(x[i])
                x[[i]][upper.tri(x[[i]])] <- ''
                list_corr2 <- reshape2::melt(x[[i]]) %>% filter(Var1 != Var2) %>% filter(value != '')
                list_comb <- merge(list_corr, list_corr2, by = c("Var1", "Var2") )
                colnames(list_comb)[i+2] <- c(names(x[i]))
                list_corr <- list_comb
            } 
        } 
        return(list_corr)
    }
    else if ("combine" %in% from){
        x[[2]][upper.tri(x[[2]])] <- ''
        
        list_mass <- reshape2::melt(x[[2]]) %>% filter(Var1 != Var2) %>% filter(value != '')
        colnames(list_mass) <- c("Var1", "Var2", "value")

        return(list_mass)
        
    }
    
}


#' @name summary_mz
#'
#' @aliases summary_mz
#'
#' @title Create a summary from adjacency list containing mass-differences with possible filter 
#'
#' @description
#' The function `summary_mz` creates a summary from the adjacency list.
#' Individual mass differences are count over all features. The input may be 
#' adjacency lists originating from the function `structural`,
#' or `combine`. The parameter `filtered` will define if data will be 
#' filtered above a certain threshold or not. 
#' 
#' @param adjacency_list
#' `data.frame`, a list of the mass-difference adjacency mateix, that has previously 
#' been generated by the function `adjacency_list(x, from = "statistical")`or 
#' `adjacency_list(x, from = "combine")`
#' 
#' @param filter
#' `number`/`FALSE`, leave empty or set to `FALSE` if unfiltered data are 
#' required. Select a `number` as a threshold where mz differences were filtered.
#' May be usefull at visualization (plotting) for big data. 
#' 
#'
#' @details
#' Summarises the adjacency list containing mass difference values, 
#' i.e. either adjacency list from structural or combine may be used.
#' Also the plots will be displayed. 
#' The default is filter = F, so the unfiltered summary will be returned. 
#' If filter is set to a `number`, e.g. 1000 only mz differences above 
#' this threshold will be displayed. 
#' The function can be applied for adjacency lists from `structural` and `combine`
#' 
#' @return 
#' `data.frame` containing the numbers of present mz differences and
#' corresponding name. Also a plot will be displayed.
#'
#' @author Liesa Salzer, \email{liesa.salzer@@helmholtz-muenchen.de}
#'
#' @examples
#' data("x_test", package = "MetNet")
#' x <- x_test[1:10, 3:ncol(x_test)]
#' x <- as.matrix(x)
#' stat_adj <- statistical(x = x, model = c("pearson", "spearman"))
#' adj_l <- adjacency_list(x = stat_adj, from = "statistical")
#' summary_mz(adjacency_list = adj_l) 
#' summary_mz(adjacency_list = adj_l, filter = 100) # filtered
#'
#'
#'@export
summary_mz <- function(adjacency_list, filter = F, ...){
    if("mass difference" %in% names(adjacency_list)){
        sum_mass <- adjacency_list %>% group_by(`mass difference`) %>% summarise(count=n()) %>%
            as.data.frame()
        sum_comb <- adjacency_list %>% group_by(`value`) %>% summarise(count=n()) %>%
            as.data.frame()  %>% add_column(sum_mass$`mass difference`)
        colnames(sum_comb) <- c("Type", "Counts", "Mass Difference")
        sum_comb <- sum_comb %>% select(Type, `Mass Difference`, Counts)}
    else{
        sum_comb <- adjacency_list %>% group_by(`value`) %>% summarise(count=n()) %>%
            as.data.frame()
        colnames(sum_comb) <- c("Type", "Counts")
        
    }
    
    
    if (filter == F) {
        plot_list <- ggplot(sum_comb, aes(x=Type, y=Counts)) + 
            geom_bar(stat = "identity") + 
            theme_minimal() + 
            coord_flip() + 
            labs(title = "Numbers of determined mass differences")
        
        plot(plot_list)
        
        return(sum_comb)
    }
    
    else if (filter != F) {
        sum_f <- filter(sum_comb,sum_comb$Counts >= filter)
        plot_list_f <- ggplot(sum_f, aes(x=Type, y=Counts)) + 
            geom_bar(stat = "identity") + 
            theme_minimal() + 
            coord_flip() + 
            labs(title = "Numbers of determined mass differences", 
                 subtitle = "(filtered)")
        
        plot(plot_list_f)
        
        return(sum_f)
        
    }
    
   
}






#' @name annotaionNames 
#'
#' @aliases annotaionNames 
#'
#' @title Add annotations to adjacency list
#'
#' @description
#' The function `annotationNames` adds available annotation values for all features
#' to the adjacency list. 
#' 
#' @param
#' `list` is the adjacency list, and `names` is a dataframe containing the feature names as rows,
#' mz values, RT values and annotations in columns
#' 
#'
#' @details
#' annotaionNames adds annotation to an adjacency list as additional column. 
#'
#' @return `data.frame` containing the the `list`-values and additional column with annotations
#'
#' @author Liesa Salzer, \email{liesa.salzer@@helmholtz-muenchen.de}
#'
#' @examples
#' data("x_test", package = "MetNet")
#' x <- x_test[1:10, 3:ncol(x_test)]
#' x <- as.matrix(x)
#' stat_adj <- statistical(x = x, model = c("pearson", "spearman"))
#' adj_l <- adjacency_list(x = stat_adj, from = "statistical")
#' annotationNames(list = adj_l, names = annotationFile)
#' 
#' @export
annotaionNames <- function (list, names) {
    Var1_names <- c()
    
    for (i in 1:length(list$Var1)) {
        
        
        value <- names[row.names(names) == list$Var1[i], "manualAnnotation"]
        value <- as.character(value)
        Var1_names <- c(Var1_names, value)
        
    }  
    Var2_names <- c()
    
    for (i in 1:length(list$Var2)) {
        
        
        value <- names[row.names(names) == list$Var2[i], "manualAnnotation"]
        value <- as.character(value)
        Var2_names <- c(Var2_names, value)
        
    }
    
    
    list$Var1_annotation <- Var1_names
    list$Var2_annotation <- Var2_names
    
    return(list)
}