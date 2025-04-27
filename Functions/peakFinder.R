# Simple peak finding function
# taken from https://github.com/stas-g/findPeaks

find_peaks <- function (x, m){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  
  # Remove closely called peaks (for some reason is not ignored by the m variable above)
  space <- m / 2
  filtered_peaks <- c(pks[1])
  for (i in 2:length(pks)) {
    if ((pks[i] - pks[i - 1]) > space) {
      filtered_peaks <- append(filtered_peaks, pks[i])
    }
  }
  
  filtered_peaks
}
