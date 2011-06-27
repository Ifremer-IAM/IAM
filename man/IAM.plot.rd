\name{IAM.plot}
\alias{IAM.plot}
\alias{IAM.barplot}
\alias{IAM.barIC}
\alias{IAM.bwplot}
\alias{IAM.zone}

\title{Plotting functions for IAM package based on lattice methods}
\description{ToDo}

\usage{
IAM.plot(formula,data,\dots)
IAM.barplot(formula,data,\dots)
IAM.barIC(formula,data,\dots)
IAM.bwplot(formula,data,\dots)
IAM.zone(formula,data,\dots)
}

\arguments{
  \item{formula}{Typical 'Lattice' formula. No 'groups' parameter, so multiple terms separated by a 'plus' sign is required for grouping (cf \emph{Extended formula interface} part in 'xyplot' documentation)}
  \item{data}{A data frame containing values for any variables in the formula (for example, an \emph{IAM.unite} output)}
  \item{\dots}{Further graphical parameters (see below for details)}
    
}

\section{Graphical parameters}{
\tabular{ll}{
\bold{type} \tab Character \cr
\bold{pch} \tab Numeric \cr
\bold{cex.pch} \tab Numeric \cr  
\bold{lty} \tab Numeric \cr
\bold{lwd} \tab Numeric \cr
\bold{relation} \tab Character \cr
\bold{key} \tab Logical \cr
\bold{pch.leg} \tab Numeric \cr
\bold{pchSize.leg} \tab Numeric \cr 
\bold{lty.leg} \tab Numeric \cr
\bold{lwd.leg} \tab Numeric \cr
\bold{txt.leg} \tab Character \cr
\bold{font.leg} \tab Numeric \cr
\bold{cex.leg} \tab Numeric \cr
\bold{space} \tab Character \cr
\bold{alpha} \tab Numeric \cr
\bold{add.v} \tab Numeric \cr
\bold{add.h} \tab Numeric \cr
\bold{add.lwd} \tab Numeric \cr
\bold{add.lty} \tab Numeric \cr
\bold{add.col} \tab Character \cr
\bold{cex.lab.x} \tab Numeric \cr
\bold{cex.lab.y} \tab Numeric \cr
\bold{font.lab.x} \tab Numeric \cr
\bold{font.lab.y} \tab Numeric \cr
\bold{cex.axis.x} \tab Numeric \cr
\bold{cex.axis.y} \tab Numeric \cr
\bold{font.axis.x} \tab Numeric \cr
\bold{font.axis.y} \tab Numeric \cr
\bold{rot.x} \tab Numeric \cr
\bold{rot.y} \tab Numeric \cr
\bold{cex.strip} \tab Numeric \cr
\bold{font.strip} \tab Numeric \cr
\bold{col.strip} \tab Character \cr
\bold{col} \tab Character \cr
\bold{fill} \tab Character \cr
\bold{xlab} \tab Character \cr
\bold{ylab} \tab Character \cr
\bold{origin} \tab Numeric \cr          
\bold{as.table} \tab Logical \cr
}
}

\details{
}

\examples{

data(iamData)
IAM.bwplot(un+deux+trois+quatre+cinq+gp~t|scen,data=DFiam_iter,txt.leg=c("1yr","2yrs","3yrs","4yrs","5yrs","6yrs"),fill="whitesmoke")

IAM.zone(un+deux+trois+quatre+cinq+gp~t|scen,data=DFiam_iter,txt.leg=c("1yr","2yrs","3yrs","4yrs","5yrs","6yrs"),alpha=0.2)

IAM.barIC(un+deux+trois+quatre+cinq+gp~t|scen,data=DFiam_iter,txt.leg=c("1yr","2yrs","3yrs","4yrs","5yrs","6yrs"))

IAM.barplot(un+deux+trois+quatre+cinq+gp~t,data=DFiam,txt.leg=c("1yr","2yrs","3yrs","4yrs","5yrs","6yrs"))

IAM.plot(un+deux+trois+quatre+cinq+gp~t,data=DFiam,txt.leg=c("1yr","2yrs","3yrs","4yrs","5yrs","6yrs"))

}

\keyword{methods}


