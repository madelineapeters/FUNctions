corr.multiloc.LD<-function(df, var=c(D), loci=c(1,2)){
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
    } #end subfunction that outputs D

    offspring_sub<-select(df, loci)

    if (length(loci) ==  2) {return(D)} else if (length(loci) == 5) {

      p.A<-sum(offspring_sub[,1])/(pop_size)
      p.B<-sum(offspring_sub[,2])/(pop_size)
      p.C<-sum(offspring_sub[,3])/(pop_size)
      p.D<-sum(offspring_sub[,4])/(pop_size)
      p.D<-sum(offspring_sub[,5])/(pop_size)

      #calculate frequency of full haplotype
      c.abcde<-0
      for (r in 1:pop_size){
        if (sum(offspring_sub[r,]) == 5) {c.abcde<-c.abcde+1}
      }
      p.ABCDE<-c.abcde/pop_size


      #calculate frequency for 4/5 loci
      c.abcd<-0
      for (r in 1:pop_size){
        if (sum(offspring_sub[r,1:4]) == 4) {c.abcd<-c.abcd+1}
      }
      p.ABCD<-c.abcd/pop_size

      #calculate frequency for 3/5 loci
      c.abc<-0
      for (r in 1:pop_size){
        if (sum(offspring_sub[r,1:3]) == 3) {c.abc<-c.abc+1}
      }
      p.ABC<-c.abc/pop_size

      D.ABCD<-multiloc.LD(offspring_sub, loci=c(1, 2, 3, 4))
      D.ABCE<-multiloc.LD(offspring_sub, loci=c(1, 2, 3, 5))
      D.ACDE<-multiloc.LD(offspring_sub, loci=c(1, 3, 4, 5))
      D.ABCE<-multiloc.LD(offspring_sub, loci=c(1, 2, 3, 5))
      D.BCDE<-multiloc.LD(offspring_sub, loci=c(2, 3, 4, 5))
      D.ABDE<-multiloc.LD(offspring_sub, loci=c(1, 2, 4, 5))

      D.ABC<-multiloc.LD(offspring_sub, loci=c(1, 2, 3))
      D.ABD<-multiloc.LD(offspring_sub, loci=c(1, 2, 4))
      D.ABE<-multiloc.LD(offspring_sub, loci=c(1, 2, 5))
      D.ACE<-multiloc.LD(offspring_sub, loci=c(1, 3, 5))
      D.ACD<-multiloc.LD(offspring_sub, loci=c(1, 3, 4))
      D.ADE<-multiloc.LD(offspring_sub, loci=c(1, 4, 5))
      D.BCD<-multiloc.LD(offspring_sub, loci=c(2, 3, 4))
      D.BCE<-multiloc.LD(offspring_sub, loci=c(2, 3, 5))
      D.BDE<-multiloc.LD(offspring_sub, loci=c(2, 4, 5))
      D.CDE<-multiloc.LD(offspring_sub, loci=c(3, 4, 5))

      D.AB<-multiloc.LD(offspring_sub, loci=c(1, 2))
      D.AC<-multiloc.LD(offspring_sub, loci=c(1, 3))
      D.AD<-multiloc.LD(offspring_sub, loci=c(1, 4))
      D.AE<-multiloc.LD(offspring_sub, loci=c(1, 5))
      D.BC<-multiloc.LD(offspring_sub, loci=c(2, 3))
      D.BD<-multiloc.LD(offspring_sub, loci=c(2, 4))
      D.BE<-multiloc.LD(offspring_sub, loci=c(2, 5))
      D.CD<-multiloc.LD(offspring_sub, loci=c(3, 4))
      D.CE<-multiloc.LD(offspring_sub, loci=c(3, 5))
      D.DE<-multiloc.LD(offspring_sub, loci=c(4, 5))

      D.5<-
        p.ABCDE - p.A*D.BCDE - p.B*D.ACDE - p.C*D.ABDE - p.D*D.ABCE - p.E*D.ABCD
      p.A*p.B*D.CDE - p.A*p.C*D.BDE - p.A*p.D*D.BCE - p.A*p.E*D.BCD -
        p.B*p.C*D.ADE - p.B*p.D*D.ACE - p.B*p.E*D.ACD -
        p.C*p.D*D.ABE - p.C*p.E*D.ABD -
        p.D*p.E*D.ABC -
        p.A*p.B*p.C*D.DE - p.A*p.B*p.D*D.CE - p.A*p.B*p.E*D.CD - p.A*p.C*p.D*D.BE - p.A*p.C*p.E*D.BD - p.A*p.D*p.E*D.BC -
        p.B*p.C*p.D*D.AE - p.B*p.D*p.E*D.AC - p.B*p.C*p.E*D.AD -
        p.C*p.D*p.E*D.AB -
        D.AB*(D.CD+D.CE+D.DE)-D.AC*(D.BD+D.BE+D.DE)-D.AD*(D.BC+D.BE+D.CE)-D.AE*(D.BC+D.BD+D.CD)-D.BC*(D.DE)-D.BD*(D.CE)-D.BE*(D.CD)
    }

    else if (length(loci) == 4) {

      p.A<-sum(offspring_sub[,1])/(pop_size)
      p.B<-sum(offspring_sub[,2])/(pop_size)
      p.C<-sum(offspring_sub[,3])/(pop_size)
      p.D<-sum(offspring_sub[,4])/(pop_size)

      #calculate frequency for 4/5 loci
      c.abcd<-0
      for (r in 1:pop_size){
        if (sum(offspring_sub[r,1:4]) == 4) {c.abcd<-c.abcd+1}
      }
      p.ABCD<-c.abcd/pop_size

      #calculate frequency for 3/5 loci
      c.abc<-0
      for (r in 1:pop_size){
        if (sum(offspring_sub[r,1:3]) == 3) {c.abc<-c.abc+1}
      }
      p.ABC<-c.abc/pop_size

      D.ABC<-multiloc.LD(offspring_sub, loci=c(1, 2, 3))
      D.ABD<-multiloc.LD(offspring_sub, loci=c(1, 2, 4))
      D.ACD<-multiloc.LD(offspring_sub, loci=c(1, 3, 4))
      D.BCD<-multiloc.LD(offspring_sub, loci=c(2, 3, 4))

      D.AB<-multiloc.LD(offspring_sub, loci=c(1, 2))
      D.AC<-multiloc.LD(offspring_sub, loci=c(1, 3))
      D.AD<-multiloc.LD(offspring_sub, loci=c(1, 4))
      D.BC<-multiloc.LD(offspring_sub, loci=c(2, 3))
      D.BD<-multiloc.LD(offspring_sub, loci=c(2, 4))
      D.CD<-multiloc.LD(offspring_sub, loci=c(3, 4))

      D.corr<-
        p.ABCD - p.A*D.BCD - p.B*D.ACD - p.C*D.ABD - p.D*D.ABC -
        p.A*p.B*D.CD - p.A*p.C*D.BD - p.A*p.D*D.BC -
        p.B*p.C*D.AD - p.B*p.D*D.AC -
        p.C*p.D*D.AB -
        D.AB*(D.CD)-D.AC*(D.BD)-D.AD*(D.BC)

    }

    else if (length(loci) == 3) {

      #calculate frequency for 3/5 loci
      c.abc<-0
      for (r in 1:pop_size){
        if (sum(offspring_sub[r,1:3]) == 3) {c.abc<-c.abc+1}
      }
      p.ABC<-c.abc/pop_size

      p.A<-sum(offspring_sub[,1])/(pop_size)
      p.B<-sum(offspring_sub[,2])/(pop_size)
      p.C<-sum(offspring_sub[,3])/(pop_size)

      D.AB<-multiloc.LD(offspring_sub, loci=c(1, 2))
      D.AC<-multiloc.LD(offspring_sub, loci=c(1, 3))
      D.BC<-multiloc.LD(offspring_sub, loci=c(2, 3))

      D.corr<-p.ABC - p.A*D.BC - p.B*D.AC - p.C*D.AB

    }

    return(D.corr)
} #end corr function
