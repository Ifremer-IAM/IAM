setMethod(
  f = "summary", signature("iamInput"),
  function(object, ...) {
    spe <- object@specific
    cat(object@desc, "(IAM input)\n")
    cat(sprintf("Simulation of %d dynamic species, %d static species and %d fleet\n",
                length(spe$Species), length(spe$StaticSpp), length(spe$Fleet) ))
    cat(rep("-",24), "\nDynamic Species | Model \n", sep = "")
    cat(sprintf("%15s | %5s\n", spe$Species, ifelse(spe$Q == 0, "XSA", "SS3")), sep ="")
    cat(rep("-",24),"\n", sep = "")
    cat(sprintf("Simulation start in %d and end in %d (%d steps)\n",
                spe$t_init, tail(spe$times,1), spe$NbSteps ))
  }
)

summary(input)


