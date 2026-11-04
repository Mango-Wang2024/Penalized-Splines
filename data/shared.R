# Input loading only: R constructs and fits its own spline model.
load_case <- function(inputs,name) {
  path <- file.path(inputs,name)
  dat <- read.csv(file.path(path,"data.csv"),check.names=FALSE)
  terms <- read.delim(file.path(path,"terms.tsv"),stringsAsFactors=FALSE)
  stopifnot(nrow(terms)>0L,!anyDuplicated(terms$predictor),
            all(grepl("^[A-Za-z][A-Za-z0-9_]*$",terms$predictor)),
            all(terms$predictor%in%names(dat)))
  knots <- setNames(lapply(terms$predictor,function(key)scan(file.path(path,paste0("knots_",key,".csv")),quiet=TRUE)),terms$predictor)
  supplied <- knots[terms$knot_mode!="default"]
  labels <- vapply(seq_len(nrow(terms)),function(j)
    sprintf("s(%s,bs='bs',m=c(3,2),k=%d)",terms$predictor[j],terms$k[j]),"")
  list(name=name,path=path,data=dat,terms=terms,knots=knots,supplied=supplied,
       formula=as.formula(paste("y~",paste(labels,collapse="+"))))
}
save_matrix <- function(M,name,path) {
  M <- as.matrix(M)
  strings <- matrix(sprintf("%.17g",as.numeric(M)),nrow(M),ncol(M))
  write.table(strings,file.path(path,paste0(name,".csv")),sep=",",quote=FALSE,row.names=FALSE,col.names=FALSE)
}
grid_inputs <- function(case) {
  keys <- case$terms$predictor
  read.csv(file.path(case$path,"grid.csv"),check.names=FALSE)[,keys,drop=FALSE]
}
