var.rank.xgb.all <- function(spltree, var.name = c("z1", "z2"),p=50) {
  var.imp.xgb.shap = spltree$imporant_variables_xgb_all
  var.imp.xgb.spltree = spltree$imporant_variables_sptree$Variable
  var.imp.xgb.spltree.gain = spltree$imporant_variables_sptree_gain$Variable
  
  ranks.var.imp.xgb.shap = ranks.var.imp.sptrees = ranks.var.imp.sptrees.gain = rep(NaN, length(var.name))
  for (i in 1:length(var.name)) {
    ranks.var.imp.xgb.shap[i] = ifelse(length(which(var.imp.xgb.shap == var.name[i])) == 0,
                                       (p+length(var.imp.xgb.shap))/2,
                                       which(var.imp.xgb.shap == var.name[i]))
    ranks.var.imp.sptrees[i] = ifelse(length(which(var.imp.xgb.spltree == var.name[i])) == 0,
                                       (p+length(var.imp.xgb.spltree))/2,
                                       which(var.imp.xgb.spltree == var.name[i]))
    ranks.var.imp.sptrees.gain[i] = ifelse(length(which(var.imp.xgb.spltree.gain == var.name[i])) == 0,
                                      (p+length(var.imp.xgb.spltree.gain))/2,
                                      which(var.imp.xgb.spltree.gain == var.name[i]))
  }
  
  return(
    list(
      sptrees.xgb = ranks.var.imp.sptrees,
      sptrees.gain = ranks.var.imp.sptrees.gain,
      xgb.shap = ranks.var.imp.xgb.shap
      # xgb.gain = ranks.var.imp.xgb.gain
    )
  )
  
}