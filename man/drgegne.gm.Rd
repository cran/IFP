\name{drgegne}
\alias{drgegne}
\title{causal models with three possible causal factors: G, G*E and E}
\description{
 provides concordance probabilities of relative pairs for a causal model with G, G*E and E components
}
\usage{
 drgegne(fdg,frg,fdge,frge,eg,e)
}
\arguments{
 \item{fdg}{an array (size=number of dominant genes+recessive genes) of dominant gene frequencies including 0 values of recessive genes of G component}
 \item{frg}{an array (size=number of dominant genes+recessive genes) of recessive gene frequencies including 0 values of dominant genes of G component}
 \item{fdge}{an array (size=number of dominant genes+recessive genes) of dominant gene frequencies including 0 values of recessive genes of G*E component}
 \item{frge}{an array (size=number of dominant genes+recessive genes) of recessive gene frequencies including 0 values of dominant genes of G*E component}
 \item{eg}{a proportion of population who are exposed to environmental cause of G*E interactiong the genetic cause of G*E during their entire life}
 \item{e}{a proportion of population who are exposed to environmental cause during their entire life}
}
\value{
 matrix of NN, ND, and DD probabilities of 9 relative pairs: 1:mzt,2:parent-offspring,3:dzt,4:sibling,5:2-direct(grandparent-grandchild),6:3rd(uncle-niece),7:3-direct(great-grandparent-great-grandchild),8:4th (causin),9:4d(great-great-grandparent-great-great-grandchild)
}
\seealso{
 drgn
 drgene
}
\examples{
### PLI=0.01.
ppt<-0.01

pg<-0.002  # the proportion of G component in total populations
pge<-0.005  # the proportion of G*E component in total populations
e<-1-(1-ppt)/(1-pg)/(1-pge)   
  # the proportion of E component in total populations

fd<-0.001  # one dominant gene
tt<-2      # the number of recessive genes

temp<-sqrt(1-((1-pg)/(1-fd)^2)^(1/tt))
fr<-c(array(0,length(fd)),array(temp,tt))
fd<-c(fd,array(0,tt))

ppe<-0.5
ppg<-pge/ppe

fdge<-0.002
ttge<-2      # the number of recessive genes

temp<-sqrt(1-((1-ppg)/(1-fdge)^2)^(1/ttge))
frge<-c(array(0,length(fdge)),array(temp,ttge))
fdge<-c(fdge,array(0,ttge))


drgegne(fd,fr,fdge,frge,ppe,e)


}
