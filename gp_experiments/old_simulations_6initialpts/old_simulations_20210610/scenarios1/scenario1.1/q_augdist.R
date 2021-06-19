f = function(x){
  1 / (x^(1/4))
}

f_augdist = function(x){
  b = 1e-10
  1 / ((x^(1/2) + b)^(1/2))
}

curve(f, from = 0, to = 5)
curve(f_augdist, add = TRUE, col = "red", lty = 2)
legend("topright", legend = c("q", "q with augdist"), 
       col = c(1, 2), lty = c(1, 2))   

f_diff = function(x){
  b = 1e-10
  (1 / (x^(1/4)) - 1 / ((x^(1/2) + b)^(1/2)))^2
}
curve(f_diff, from = 0, to = 0.5)
