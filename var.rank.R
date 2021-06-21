var.rank <- function(spltree, var.name="s1"){
  
  var.imp.xgb.gain <- as.vector(xgboost::xgb.importance(spltree$xgbModel$feature_names, spltree$xgbModel)[,"Feature"])
  
  var.imp.xgb.shap = spltree$imporant_variables_xgb_all
  var.imp.sptrees = names(spltree$imporant_variables_sptree)
  var.imp.sptrees.gain = names(spltree$imporant_variables_sptree_gain)
  
  ranks.var.imp.xgb.gain = ifelse(length(which(var.imp.xgb.gain==var.name))==0, 0, which(var.imp.xgb.gain==var.name))
  ranks.var.imp.xgb.shap = ifelse(length(which(var.imp.xgb.shap==var.name))==0, 0, which(var.imp.xgb.shap==var.name))
  ranks.var.imp.sptrees = ifelse(length(which(var.imp.sptrees==var.name))==0, 0, which(var.imp.sptrees==var.name))
  ranks.var.imp.sptrees.gain = ifelse(length(which(var.imp.sptrees.gain==var.name))==0, 0, which(var.imp.sptrees.gain==var.name))
  
  return(list(sptrees.xgb = ranks.var.imp.sptrees, 
              sptrees.gain = ranks.var.imp.sptrees.gain, 
              xgb.shap= ranks.var.imp.xgb.shap,
              xgb.gain= ranks.var.imp.xgb.gain))
  
}