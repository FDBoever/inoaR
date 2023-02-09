#functions

#' Title
#'
#' @param data
#' @param filename
#'
#' @return
#' @export
#'
#' @examples
writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,1], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,2]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}
