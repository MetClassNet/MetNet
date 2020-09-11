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
#' Changes to MetNet: New attributes added, if model = "combined" (default) the result will be the same 
#' as in MetNet, except the output list items were named to "combined" and "Character" 
#' if model = "pearson" or "spearman" than also the corresponding weighted statistical 
#' adjacency matrix is required as attribute (weighted_statistical = XY)
#' The output in this case will be a list containing 4 listitems, where combination relies on the
#' unweighted adjacency matrix of either Pearson or Spearman. 
#' Moreover corresponding Correlation and p-values will be displayes as listitems
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
combine <- function(structural, statistical, threshold = 1, model = "combined", weighted_statistical) {
    
    ## Is changed since structural list is now lenght 3
    if (!is.list(structural) | length(structural) != 3)
        stop("structural is not a list of length 3")
    
    if (!is.matrix(structural[[1]]) | !is.numeric(structural[[1]]))
        stop("structural[[1]] is not a numeric matrix")
    
    if (!is.matrix(structural[[2]]) | !is.character(structural[[2]]))
        stop("structural[[2]] is not a character matrix")
    
    ## Additional to MetNet:
    if (!is.matrix(statistical[[3]]) | !is.numeric(statistical[[3]]))
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
    
    
    ## Changes to MetNet: Distinguish between default (model == "combined") and model = "pearson" or "spearman"
    if(model == "pearson"){
        
        ## create the first entry of the list
        ## sum the matrices structural and statistical, if the value is above
        ## threshold then assign 1, otherwise 0
        cons_num <- structural[[1]] + statistical[["pearson"]]
        cons_num <- ifelse(cons_num > threshold, 1, 0)
        
        ## if p-values have previously been calculated
        if("Correlation Value" %in% names(weighted_statistical[["pearson"]][1])){
            cons_corr <- ifelse(cons_num == 1, weighted_statistical[["pearson"]][["Correlation Value"]], "")
            cons_p <- ifelse(cons_num == 1, weighted_statistical[["pearson"]][["p-Value"]], "")}
        else {
            cons_corr <- ifelse(cons_num == 1, weighted_statistical[["pearson"]], "")
            cons_p <- NaN
        }}
    
    if(model == "spearman"){
        ## create the first entry of the list
        ## sum the matrices structural and statistical, if the value is above
        ## threshold then assign 1, otherwise 0
        cons_num <- structural[[1]] + statistical[["spearman"]]
        cons_num <- ifelse(cons_num > threshold, 1, 0)
        
        ## if p-values have previously been calculated
        if("Correlation Value" %in% names(weighted_statistical[["spearman"]][1])){
            cons_corr <- ifelse(cons_num == 1, weighted_statistical[["spearman"]][["Correlation Value"]], "")
            cons_p <- ifelse(cons_num == 1, weighted_statistical[["spearman"]][["p-Value"]], "")}
        else {
            cons_corr <- ifelse(cons_num == 1, weighted_statistical[["spearman"]], "")
            cons_p <- NaN
        }}
    
    
    if(model == "combined"){
        ## create the first entry of the list
        ## sum the matrices structural and statistical, if the value is above
        ## threshold then assign 1, otherwise 0
        cons_num <- structural[[1]] + statistical[[3]]
        cons_num <- ifelse(cons_num > threshold, 1, 0)
        cons_corr <- NaN
        cons_p <- NaN
    }
    
    ## create the second entry of the list
    ## if element in cons_num is equal to 1, take the element in structural[[2]]
    ## (the type of link), otherwise ""
    cons_char <- ifelse(cons_num == 1, structural[[2]], "")
    
    ## assign to list
    ## Compared to MetNet names were assigned
    res[[model]] <- cons_num
    res[["Character"]] <- cons_char
    
    ## assign Correlation and p-values to list if model is "pearson" or "spearman"
    if(model == "pearson" || model == "spearman"){
        res[["Correlation Value"]] <- cons_corr
        res[["p-Value"]] <- cons_p
    }
    
    
    return(res)
}


#' exportNet2gml is a function that exports adjacency matrices to gml using igraph
#' Needs following attributes:
#' x: adjacency matrix that needs to be exported
#' from: originated from wich function, possible values are 
#' - "structural" Produces a gml file with edge attributes containing mass difference values, data saved as "structural_type.gml"
#' - "statistical+p" produces a gml file with edge attributes containing correlation values and p-values, saved as "statistical.'model'.gml"
#' (TO DO: ADD "statistical")
#' - "combine" produces a gml file with edge attributes containing correlation and p-values for pearson/ spearman correlations, 
#' saved as "combined.gml"
exportNet2gml <- function (x, from, adjacency_list, ...) {
    #Check arguments
    if (!from %in% c("structural", "statistical+p", "combine"))
        stop("type not in 'structural', 'statistical+p', 'combine'")
    
    if ("structural" %in% from) {
        
        mat <- x[[1]]
        net       <- 
            graph_from_adjacency_matrix(mat, mode = "undirected")
        E(net)$sourceID <-  as.character(adjacency_list$`Var1`)
        E(net)$targetID <-  as.character(adjacency_list$`Var2`)
        E(net)$`mass difference` <-  adjacency_list$`mass difference`
        
        write_graph(net, "structural.gml", format = c("gml"))
    }
    else if ("statistical+p" %in% from) {
        
        for (i in 1:length(x)) {
            cor_list <- x[[i]]
            ##Plot structural adjacency matrix and export to gml
            net_cor <-
                igraph::graph_from_adjacency_matrix(cor_list[[1]], mode = "undirected", weighted = T)
            net_p   <-
                igraph::graph_from_adjacency_matrix(cor_list[[2]], mode = "undirected", weighted = T)
            net_comb <- union(net_cor, net_p)
            names(edge_attr(net_comb))[1] <- "correlation"
            names(edge_attr(net_comb))[2] <- "p"
            # #net_plot <- plot(net_type, edge.width = 5, vertex.label.cex = 0.5, edge.color = "grey")
            q <- names(x[i])
            #E(net_cor)$`correlation_tst` <-  adjacency_list[,"pearson"]
            
            # rownames(cor_list[[1]])
            write_graph(net_comb, file = sprintf('statistical.%s.gml', q), format = c("gml"))
            
        }
        
        
        net <- graph_from_data_frame(adjacency_list, directed = TRUE, vertices = NULL)
        E(net)$sourceID <-  as.character(adjacency_list$`Feature1`)
        E(net)$targetID <-  as.character(adjacency_list$`Feature2`)
        write_graph(net, file = "statistical.gml", format = c("gml"))
        
        
    }
    
    else if ("combine" %in% from) {
        if ("pearson" %in% names(x[1]) | "spearman" %in% names(x[1]) ){
            net <- graph_from_adjacency_matrix(x[[1]], mode = "undirected")
            E(net)$sourceID <-  as.character(adjacency_list$`Var1`)
            E(net)$targetID <-  as.character(adjacency_list$`Var2`)
            E(net)$massdifference <-  adjacency_list$`value`
            E(net)$correlation <-  adjacency_list$`Correlation Value`
        }
        else { #if "combined" or other model
            net       <- 
                graph_from_adjacency_matrix(x[[1]], mode = "undirected", weighted = T)
            
        }
        write_graph(net, "combined.gml", format = c("gml"))
        
        
    }
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
#' on the source `from` of the adjacency matrix.
#' 
#'
#' @param
#' x is an adjacency matrix, that has previously been generated by the functions 
#' `structural`, `statistical`, or `combine`
#' 
#
#' @param
#' `from` defines the function from which the adjacency matrix originated.
#' This can be either `structural`, `statistical`, or `combine`
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
        
        list_type <- melt(x[[2]]) %>% filter(Var1 != Var2) %>% filter(value != '')
        list_mass <- melt(x[[3]]) %>% filter(Var1 != Var2) %>% filter(value != '')
        combine <- add_column(list_type,  `mass difference`= list_mass$value) %>% as.data.frame()
        return(combine)
    }
    else if ("statistical" %in% from) {
        
        for (i in seq_along(x)) {
            if (i == 1) {
                x[[i]][upper.tri(x[[i]])] <- ''
                list_corr <- melt(x[[i]]) %>% filter(Var1 != Var2) %>% filter(value != '') %>% 
                    select(Var1, Var2, value) 
                colnames(list_corr) <- c("Var1", "Var2", names(x[i]))
            }
            if (i != 1){
                model = names(x[i])
                x[[i]][upper.tri(x[[i]])] <- ''
                list_corr2 <- melt(x[[i]]) %>% filter(Var1 != Var2) %>% filter(value != '')
                list_comb <- add_column(list_corr, list_corr2$value)
                list_comb <- as.data.frame(list_comb) 
                colnames(list_comb)[i+2] <- c(names(x[i]))
                list_corr <- list_comb
            } 
        } 
        return(list_corr)
    }
    else if ("combine" %in% from){
        x[[2]][upper.tri(x[[2]])] <- ''
        x[[3]][upper.tri(x[[3]])] <- ''
        x[[4]][upper.tri(x[[4]])] <- ''
        
        list_mass <- melt(x[[2]]) %>% filter(Var1 != Var2) %>% filter(value != '')
        list_corr <- melt(x[[3]]) %>% filter(Var1 != Var2) %>% filter(value != '')
        list_p    <- melt(x[[4]]) %>% filter(Var1 != Var2) %>% filter(value != '')
        listed <- add_column(list_mass, `Correlation Value` = list_corr$value)
        listed <- add_column(listed, `p-Value` = list_p$value)
        return(listed)
        
    }
    
}

#' sum_mass summarises the adjacency list containing mass difference values, 
#' i.e. either adjacency list from structural or combine may be used
sum_mass <- function(adjacency_list){
    
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
    
    plot_list <- ggplot(sum_comb, aes(x=Type, y=Counts, fill=Type)) + geom_bar(stat = "identity") + theme_minimal() + 
        labs(title = "Numbers of destinct type of a biochemical reaction")+ scale_fill_brewer(palette = "Blues") + theme(legend.position = "right")
    plot(plot_list)
    
    return(sum_comb)
}




# Function to add annotation names
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