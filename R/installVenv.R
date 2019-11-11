#' Install UMI4Cats python virtual environment
#'
#' Automatically install python virtual environment and all necessary packages from \code{python/requirements.txt}.
#' @param path_venv Path where to install the virtual environment.
#' @details Automatilly creates a python virtual environment inside the defined \code{path_env}. User needs to have
#' \link[=https://www.python.org/download/releases/3.0/]{Python 3} installed and added to your path as \code{python 3}.
#' @examples
#' \dontrun{
#'  installVenv('~/venvs/UMI4Cats_venv/')
#' }
#' @export
installVenv <- function(path_venv){
  # TODO: Create & install python packages using {reticulate}

  # define paths
  requirements <- system.file("python/requirements.txt", package = "UMI4Cats")

  # TODO: Check if venv already exists

  # create virtualenv
  system(paste('python3 -m venv',
               file.path(path_venv)))

  # upgrade pip
  system(paste(file.path(path_venv, 'bin', 'pip'),
               'install --upgrade pip'))

  # install cython module first
  system(paste(file.path(path_venv, 'bin', 'pip3'),
               'install Cython'))

  # install required modules
  system(paste(file.path(path_venv, 'bin', 'pip3'),
               'install -r', requirements))

  # TODO: Return TRUE of FALSE if the venv creation was successful
}
