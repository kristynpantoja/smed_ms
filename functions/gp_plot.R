

plotGP = function(
  x.seq, function.values, x.data, y.data, kernel, length.scale, period,
  xmin = 0, xmax = 1, signal.var, measurement.var = NULL
){
  gpfit = getGPPredictive(
    x = x.seq, x.input = x.data, y.input = y.data, type = kernel, 
    l = length.scale, p = period, signal.var = signal.var, 
    measurement.var = measurement.var)
  
    ggstar = data.table(
      x = x.seq, 
      p_mean = gpfit$pred_mean,
      p_sd = sqrt(diag(gpfit$pred_var))
    )
    
    ggdata = data.table(
      x = x.data,
      y = y.data
    )
    
    plt = ggplot() + 
      geom_line(
        data = data.frame(x = x.seq, y = function.values), aes(x = x, y = y)) +
      geom_line(data = ggstar, aes(x = x, y = p_mean), color = 3) +
      geom_ribbon(
        data = ggstar, aes(
          x = x, ymin = p_mean - 2 * p_sd,ymax = p_mean + 2 * p_sd),
        fill = 3, alpha = 0.2) +
      geom_point(data = ggdata, aes(x = x, y = y))
  
  return(plt)
}
