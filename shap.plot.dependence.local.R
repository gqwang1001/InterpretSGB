shap.plot.dependence.local=
  function (data_long, x, y = NULL, color_feature = NULL, data_int = NULL, 
          dilute = FALSE, smooth = TRUE, size0 = NULL, add_hist = FALSE, 
          add_stat_cor = FALSE, alpha = NULL, jitter_height = 0, jitter_width = 0, 
          ...) 
{
  if (is.null(y)) 
    y <- x
  data0 <- data_long[variable == y, .(variable, value)]
  data0$x_feature <- data_long[variable == x, rfvalue]
  if (!is.null(color_feature) && color_feature == "auto") {
    color_feature <- strongest_interaction(X0 = data0, Xlong = data_long)
  }
  if (!is.null(color_feature)) {
    data0$color_value <- data_long[variable == color_feature, 
                                   rfvalue]
  }
  if (!is.null(data_int)) 
    data0$int_value <- data_int[, x, y]
  nrow_X <- nrow(data0)
  if (is.null(dilute)) 
    dilute = FALSE
  if (dilute != 0) {
    dilute <- ceiling(min(nrow(data0)/10, abs(as.numeric(dilute))))
    set.seed(1234)
    data0 <- data0[sample(nrow(data0), min(nrow(data0)/dilute, 
                                           nrow(data0)/2))]
  }
  if (x == "dayint") {
    data0[, `:=`(x_feature, as.Date(data0[, x_feature], 
                                    format = "%Y-%m-%d", origin = "1970-01-01"))]
  }
  if (is.null(size0)) {
    size0 <- if (nrow(data0) < 1000L) 
      1
    else 0.4
  }
  if (is.null(alpha)) {
    alpha <- if (nrow(data0) < 1000L) 
      1
    else 0.6
  }
  plot1 <- ggplot(data = data0, aes(x = x_feature, y = if (is.null(data_int)) 
    value
    else int_value, color = if (!is.null(color_feature)) 
      color_value
    else NULL)) + 
    geom_jitter(size = size0, width = jitter_width, 
                              height = jitter_height, alpha = alpha, ...) +
    labs(y = if (is.null(data_int)) 
                                paste0("SHAP value for ", label.feature(y))
                                else paste0("SHAP interaction values for\n", label.feature(x), 
                                            " and ", label.feature(y)), x = label.feature(x), color = if (!is.null(color_feature)) 
                                              paste0(label.feature(color_feature), "\n", "(Feature value)")
                                else NULL) + 
    # scale_color_gradient(low = "blue", high = "red") +
    scale_color_gradientn(colours = rainbow(5))+
    theme_bw() + 
    theme(legend.position = "right", legend.title = element_text(size = 10), 
                       legend.text = element_text(size = 8))
  if (smooth) {
    plot1 <- plot1 + geom_smooth(method = "loess", color = "red", 
                                 size = 0.4, se = FALSE)
  }
  plot1 <- plot.label(plot1, show_feature = x)
  if (add_stat_cor) {
    plot1 <- plot1 + ggpubr::stat_cor(method = "pearson")
  }
  if (add_hist) {
    plot1 <- ggExtra::ggMarginal(plot1, type = "histogram", 
                                 bins = 50, size = 10, color = "white")
  }
  plot1
  }

label.feature=function (x) 
{
  labs <- SHAPforxgboost::labels_within_package
  if (!is.null(new_labels)) {
    if (!is.list(new_labels)) {
      message("new_labels should be a list, for example,`list(var0 = 'VariableA')`.\n")
    }
    else {
      message("Plot will use your user-defined labels.\n")
      labs = new_labels
    }
  }
  out <- rep(NA, length(x))
  for (i in 1:length(x)) {
    if (is.null(labs[[x[i]]])) {
      out[i] <- x[i]
    }
    else {
      out[i] <- labs[[x[i]]]
    }
  }
  return(out)
}

plot.label=function (plot1, show_feature) 
{
  if (show_feature == "dayint") {
    plot1 <- plot1 + scale_x_date(date_breaks = "3 years", 
                                  date_labels = "%Y")
  }
  else if (show_feature == "AOT_Uncertainty" | show_feature == 
           "DevM_P1km") {
    plot1 <- plot1 + scale_x_continuous(labels = function(x) paste0(x * 
                                                                      100, "%"))
  }
  else if (show_feature == "RelAZ") {
    plot1 <- plot1 + scale_x_continuous(breaks = c((0:4) * 
                                                     45), limits = c(0, 180))
  }
  plot1
}

