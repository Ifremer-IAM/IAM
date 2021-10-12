#' If x is `NULL`, return y, otherwise return x
#'
#' from ThinkR style
#'
#' @param x,y Two elements to test, one potentially `NULL`
#'
#' @noRd
#'
#'
#' @examples
#' NULL %||% 1
"%||%" <- function(x, y){
  if (is.null(x)) {
    y
  } else {
    x
  }
}

#' return `TRUE` if in `devloppment mode`
#'
#' from ThinkR style
#'
#' @export
app_dev <- function(){
  getOption( "IAM.dev" ) %||% FALSE
}
