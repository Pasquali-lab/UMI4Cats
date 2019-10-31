#' umi4Cats virtualenv
#'
#'@description
#' Install automatically virtualenv from requirements.txt
#'
#'@usage
#'umi4catsVenvInst(pathVenv)
#'@param pathVenv Path where to install the virtualenv
#'
#'\dontrun{
#'@examples
# umi4catsVenvInst('/imppc/labs/lplab/share/env/umi4CatsVenv')
#'}
#'
#'@export


umi4catsVenvInst <- function(pathVenv){

  # define paths
  requirements <-  system.file("python/requirements.txt", package = "UMI4Cats")

  # create virtualenv
  system(paste('python3 -m venv',
               file.path(pathVenv)))

  # upgrade pip
  system(paste(file.path(pathVenv, 'bin', 'pip'),
               'install --upgrade pip'))

  # install cython module first
  system(paste(file.path(pathVenv, 'bin', 'pip3'),
               'install Cython'))

  # install required modules
  system(paste(file.path(pathVenv, 'bin', 'pip3'),
               'install -r', requirements))

}
