num_to_roman <- function(num) {
  split <- substr(num, 4, nchar(num))
  roman <- as.roman(split)
  string <- paste0("chr", roman)
  return(string)
}

roman_to_num <- function(roman) {
  split <- substr(roman, 4, nchar(roman))
  num <- as.numeric(as.roman(split))
  string <- paste0("chr", num)
  return(string)
}
