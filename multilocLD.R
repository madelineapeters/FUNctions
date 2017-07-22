multiloc.LD<-function(df, var=c(D), loci=c(1,2)){
  if (length(loci) == 2) {
    offspring_sub<-select(df, loci)
    p.a<-sum(offspring_sub[,1])/(dim(df)[1]) #probability of having A allele
    p.b<-sum(offspring_sub[,2])/(dim(df)[1])  #probability of having B allele
    c.ab<-0
    for (r in 1:dim(df)[1]){
      if (sum(offspring_sub[r,]) == 2) {c.ab<-c.ab+1}
    }
    p.ab<-c.ab/dim(df)[1]

    D<-(p.ab-(p.a*p.b)) #disequilibrium coefficient
    r<-D/sqrt(p.a*(1-p.a)*p.b*(1-p.b)) #correlation coefficient

  } else if (length(loci) == 3) {

    offspring_sub<-select(df, loci)
    p.a<-sum(offspring_sub[,1])/(dim(df)[1]) #probability of having A allele
    p.b<-sum(offspring_sub[,2])/(dim(df)[1])  #probability of having B allele
    p.c<-sum(offspring_sub[,3])/(dim(df)[1])  #probability of having C allele
    c.abc<-0
    for (r in 1:dim(df)[1]){
      if (sum(offspring_sub[r,]) == 3) {c.abc<-c.abc+1}
    }
    p.abc<-c.abc/dim(df)[1]

    D<-(p.abc-(p.a*p.b*p.c)) #disequilibrium coefficient
    r<-D/sqrt(p.a*(1-p.a)*p.b*(1-p.b)*p.c*(1-p.c)) #correlation coefficient
  } else if (length(loci) == 4) {
    offspring_sub<-select(df, loci)
    p.a<-sum(offspring_sub[,1])/(dim(df)[1]) #probability of having A allele
    p.b<-sum(offspring_sub[,2])/(dim(df)[1]) #probability of having B allele
    p.c<-sum(offspring_sub[,3])/(dim(df)[1]) #probability of having C allele
    p.d<-sum(offspring_sub[,4])/(dim(df)[1]) #probability of having D allele
    c.abcd<-0
    for (r in 1:dim(df)[1]){
      if (sum(offspring_sub[r,]) == 4) {c.abcd<-c.abcd+1}
    }
    p.abcd<-c.abcd/dim(df)[1]

    D<-(p.abcd-(p.a*p.b*p.c*p.d)) #disequilibrium coefficient
    r<-D/sqrt(p.a*(1-p.a)*p.b*(1-p.b)*p.c*(1-p.c)*p.d*(1-p.d)) #correlation coefficient
  } else if (length(loci) == 5) {
    offspring_sub<-select(df, A, B, C, D, E)
    p.a<-sum(offspring_sub[,1])/(dim(df)[1]) #probability of having A allele
    p.b<-sum(offspring_sub[,2])/(dim(df)[1]) #probability of having B allele
    p.c<-sum(offspring_sub[,3])/(dim(df)[1]) #probability of having C allele
    p.d<-sum(offspring_sub[,4])/(dim(df)[1]) #probability of having D allele
    p.e<-sum(offspring_sub[,5])/(dim(df)[1]) #probability of having E allele
    c.abcde<-0
    for (r in 1:dim(df)[1]){
      if (sum(offspring_sub[r,]) == 5) {c.abcde<-c.abcde+1}
    }
    p.abcde<-c.abcde/dim(df)[1]

    D<-(p.abcde-(p.a*p.b*p.c*p.d*p.e)) #disequilibrium coefficient
    r<-D/sqrt(p.a*(1-p.a)*p.b*(1-p.b)*p.c*(1-p.c)*p.d*(1-p.d)*p.e*(1-p.e)) #correlation coefficient
  } else {var<-"Linkage can only be calculated for between 0 and 5 loci. Please adjust 'loci'. Default includes loci 1 and 2"}

  out<-var
  return(out)
} #end function
