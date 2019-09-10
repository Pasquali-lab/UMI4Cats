#' @export
setMethod("dgram", "UMI4C",
          function(x, ...) {
            out <- x@dgram

            out
          })

#' @exportMethod "dgram<-"
setReplaceMethod("dgram",
                 signature(object="UMI4C", value="list"),
                 function(object, value) {

                   object@dgram <- value
                   return(object)
                 })

#' @export
setMethod("trend", "UMI4C",
          function(x, ...) {
            out <- x@trend

            out
          })

#' @exportMethod "trend<-"
setReplaceMethod("trend",
                 signature(object="UMI4C", value="data.frame"),
                 function(object, value) {

                   object@trend <- value
                   return(object)
                 })
