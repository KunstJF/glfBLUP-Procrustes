procustes_smoother <- function(rotation_matrix){
  m <- dim(rotation_matrix)[1]
  perm_mat = matrix(0, nrow = m, ncol = m)
  for (k in 1:m){
    k_plus <- which.max(rotation_matrix[, k])
    k_minus <- which.min(rotation_matrix[, k])
    r_plus <- max(rotation_matrix[k, ])
    r_minus <- min(rotation_matrix[k, ])
    if (r_minus < r_plus & r_minus > 0 & r_plus > 0){
      perm_mat[k_plus, k] = 1
    } else if (abs(r_minus) < r_plus & r_minus < 0 & r_plus > 0){
      perm_mat[k_plus, k]  = 1
    } else if (abs(r_minus) > r_plus & r_minus < 0 & r_plus > 0){
      perm_mat[k_minus, k] = -1
    } else{
      perm_mat[k_minus, k] = -1
    }
  }
  return(perm_mat)
}
