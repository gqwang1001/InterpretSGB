var.rank.xgb.all <- function(spltree, var.name = c("z1", "z2"),p=50) {
  
  ranks.var.imp.xgb.shap =
    ranks.var.imp.xgb.gain =
    ranks.var.imp.sptrees = 
    ranks.var.imp.lm = Variable
    ranks.var.imp.sptrees.gain =
    ranks.var.imp.lm.gain=
    rep(NaN, length(var.name))
  
  var.imp.xgb.gain = spltree$imporant_variables_gain$Feature
  var.imp.xgb.shap = spltree$imporant_variables_xgb_all$Variable
  var.imp.xgb.spltree = spltree$simplified_tree$vi$Variable
  var.imp.xgb.spltree.gain = spltree$simplified_tree_gain$vi$Variable
  var.imp.xgb.spltree.lm = spltree$lmModel$vi$data$Variable
  var.imp.xgb.spltree.lm.gain = spltree$lmModel_gain$vi$data$
  
  for (i in 1:length(var.name)) {
    ranks.var.imp.xgb.shap[i] = ifelse(length(which(var.imp.xgb.shap == var.name[i])) == 0,
                                       (p+length(var.imp.xgb.shap))/2,
                                       which(var.imp.xgb.shap == var.name[i]))
    ranks.var.imp.xgb.gain[i] = ifelse(length(which(var.imp.xgb.gain == var.name[i])) == 0,
                                       (p+length(var.imp.xgb.gain))/2,
                                       which(var.imp.xgb.gain == var.name[i]))
    
    ranks.var.imp.sptrees[i] = ifelse(length(which(var.imp.xgb.spltree == var.name[i])) == 0,
                                       (p+length(var.imp.xgb.spltree))/2,
                                       which(var.imp.xgb.spltree == var.name[i]))
    ranks.var.imp.sptrees.gain[i] = ifelse(length(which(var.imp.xgb.spltree.gain == var.name[i])) == 0,
                                      (p+length(var.imp.xgb.spltree.gain))/2,
                                      which(var.imp.xgb.spltree.gain == var.name[i]))
    ranks.var.imp.lm[i] = ifelse(length(which(var.imp.xgb.spltree.lm == var.name[i])) == 0,
                                 (p+length(var.imp.xgb.spltree.lm))/2,
                                 which(var.imp.xgb.spltree.lm == var.name[i]))
    ranks.var.imp.lm.gain[i] = ifelse(length(which(var.imp.xgb.spltree.lm.gain == var.name[i])) == 0,
                                 (p+length(var.imp.xgb.spltree.lm.gain))/2,
                                 which(var.imp.xgb.spltree.lm.gain == var.name[i]))
  }
  
  return(
    list(
      sptrees.xgb = ranks.var.imp.sptrees,
      sptrees.gain = ranks.var.imp.sptrees.gain,
      xgb.shap = ranks.var.imp.xgb.shap,
      xgb.gain = ranks.var.imp.xgb.gain,
      lm = ranks.var.imp.lm,
      lm.gain = ranks.var.imp.lm.gain)
  )
  
}