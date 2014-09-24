setClass("fd+", slots = c(coefs = "matrix", basis = "basis+"),
         validity = function(object) {
             if(ncol(object$coefs) != object@basis@ncoef)
                 return("ncol(coefs) must be equal to the number of basis functions")
             return(TRUE)
         }
)

