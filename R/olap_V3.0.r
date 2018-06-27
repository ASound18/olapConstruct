##-------------------------------------------------------
##    Discovery-driven Exploration of OLAP Data Cubes
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


olapConstruct <- function(data, measure, dimensions, user_agg_function="mean", exception_tau=2.5, one_side_trim=0, output_dir, output_csv=FALSE){

  ##-------------------------
  ##        OLAP
  ##~~~~~~~~~~~~~~~~~~~~~~~~~

  ## user-aggregate functions
  if(user_agg_function=="mean"){
    user_agg_function = function(x, ...) UseMethod("mean")
  }else{
    user_agg_function = function(x, ...) .Primitive("sum")
  }

  ## generate unique combinations from all possible values in dimensions
  for(i in 1:length(dimensions)){
    assign(paste0("dim",i), c(unique(data[,dimensions[i]]), "ALL"))
  }
  levels_Grid <- expand.grid(lapply(1:length(dimensions), function(x) get(paste0("dim",x))), stringsAsFactors = FALSE)
  names(levels_Grid) <- dimensions

  ## add layer and group-by
  levels_Grid$layer <- as.character(apply(levels_Grid, 1, function(x) length(dimensions)+1-sum(x[1:length(dimensions)]=="ALL")))
  levels_Grid$groupby <- apply(levels_Grid, 1, function(x) paste(dimensions[x[1:length(dimensions)]!="ALL"], collapse = ", "))

  ## initialize computed, aggregated value, coefficient, anticipated value, exception
  levels_Grid$computed <- "N"
  levels_Grid$user_agg <- NA
  levels_Grid$g <- NA
  levels_Grid$residual <- NA
  levels_Grid$standardized_residual <- NA
  # levels_Grid$exception <- NA

  ## initialize base-level user-agg from data
  temp <- merge(levels_Grid, data, by = dimensions, all.x = T)
  temp$user_agg <- temp[,get("measure")]
  levels_Grid <- temp[,names(temp)!=get("measure")]
  levels_Grid$g <- levels_Grid$user_agg

  ## sort by layers (necessary)
  levels_Grid <- levels_Grid[order(levels_Grid$layer, levels_Grid$groupby, decreasing = TRUE),]

  ## OLAP function
  olapCompute <- function(groupby, levels_Grid){
    ## prepare
    curr_dim = strsplit(groupby, ", ")[[1]]
    levels_Grid$residual[levels_Grid$groupby==groupby] = levels_Grid$g[levels_Grid$groupby==groupby]
    ## mark current groupby as computed
    levels_Grid$computed[levels_Grid$groupby==groupby] = "Y"
    ## loop through all immediate children (not marked computed)
    for(i in length(curr_dim):1) {
      ## (inverse) prepare
      curr_dim_child = curr_dim[-i]
      groupby_child = paste(curr_dim_child, collapse=", ")
      if(unique(levels_Grid$computed[levels_Grid$groupby==groupby_child])=="Y") next
      ## compute user-agg
      levels_Grid$user_agg[levels_Grid$groupby==groupby_child] = apply(levels_Grid[levels_Grid$groupby==groupby_child,], 1, function(x) user_agg_function(merge(data, t(x[1:length(dimensions)]), by=curr_dim_child)[,get("measure")], na.rm=TRUE) )
      ## compute avg-g
      avg_g = aggregate(as.formula(paste("residual ~", ifelse(length(curr_dim_child)!=0, paste(curr_dim_child, collapse="+"), "1"))), data = levels_Grid[levels_Grid$groupby==groupby,], FUN = mean, trim = one_side_trim)
      names(avg_g) = c(names(avg_g)[-ncol(avg_g)],"avg")
      ## subtract the avg-g from g
      #### update parent-groupby
      temp1 = merge(levels_Grid[levels_Grid$groupby==groupby,], avg_g, by = curr_dim_child, all.x = TRUE)
      temp1$residual = temp1$residual - temp1$avg
      #### update child-groupby
      temp2 = merge(levels_Grid[levels_Grid$groupby==groupby_child,], avg_g, by = curr_dim_child, all.x = TRUE)
      temp2$g = temp2$avg
      #### resember levels_grid
      levels_Grid = rbind(levels_Grid[!levels_Grid$groupby %in% c(groupby,groupby_child), ], temp1[,names(temp1)!="avg"], temp2[,names(temp2)!="avg"])
      levels_Grid = levels_Grid[order(levels_Grid$layer, levels_Grid$groupby, decreasing = TRUE),]
    }
    ## derive standardized residual and exception at the current groupby
    residual_std = sqrt(var(levels_Grid$residual[levels_Grid$groupby==groupby], na.rm = TRUE))
    levels_Grid$standardized_residual[levels_Grid$groupby==groupby] = abs(levels_Grid$residual[levels_Grid$groupby==groupby]) / residual_std
    # levels_Grid$exception[levels_Grid$groupby==groupby] = levels_Grid$standardized_residual[levels_Grid$groupby==groupby] > exception_tau

    ## loop through the immediate children above to call olapCompute recursively
    for(j in length(curr_dim):1) {
      ## (inverse) prepare
      curr_dim_child = curr_dim[-j]
      groupby_child = paste(curr_dim_child, collapse=", ")
      ## recursive call
      if(groupby_child != "") levels_Grid = olapCompute(groupby=groupby_child, levels_Grid=levels_Grid)
    }
    return(levels_Grid)
  }

  ## main algorithm to call the core function
  unique_groupby <- unique(levels_Grid[,"groupby"])
  base_groupby <- unique_groupby[1]
  ans = olapCompute(groupby=base_groupby, levels_Grid=levels_Grid)


  ## define SelfExp, InExp, and PathExp
  #### SelfExp: the same as the standardized residual
  #### InExp: maximum of all SelfExp underneath the cell (should be calculated from bottom to top)
  #### PathExp: maximum of all SelfExp for each groupby
  ans$SelfExp <- ans$standardized_residual
  ans$InExp <- NA
  ans$PathExp <- NA
  max_layer <- max(ans$layer, na.rm = TRUE)
  ## calculate InExp
  for(i in 1:nrow(ans)) {
    curr_layer = ans$layer[i]
    if(curr_layer == max_layer) next
    ## merge with the one-layer-lower layer to get their SelfExp
    temp = merge(ans[ans$layer==as.character(as.numeric(curr_layer)+1),], ans[i,dimensions], by=dimensions[ans[i,dimensions]!="ALL"])
    ans$InExp[i] = max(temp$SelfExp, na.rm=TRUE)
  }
  ## calculate PathExp
  for(i in 1:length(unique_groupby)) {
    curr_groupby = unique_groupby[i]
    ans$PathExp[ans$groupby==curr_groupby] = max(c(ans$SelfExp[ans$groupby==curr_groupby], -9999), na.rm=TRUE)
  }

  ## remove rows whose user-agg is NA
  ans = ans[!is.na(ans$user_agg),]

  ## output to csv
  if(output_csv) write.csv(ans, paste0(output_dir, "olap_output.csv"))
  saveRDS(ans, paste0(output_dir, "olap_output.rds"))

  ## return ans
  return(ans)
}




