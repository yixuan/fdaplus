setClassUnion("vecmat", c("numeric", "matrix"))

setClass("fd+", slots = c(coefs = "vecmat",
                          basis = "basis+")
)

