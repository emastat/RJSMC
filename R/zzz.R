.onUnload <- function (libpath) {
  # Whenever you use C++ code in your package, you need to clean up after yourself when your package is unloaded
  # see http://r-pkgs.had.co.nz/src.html#cpp
  library.dynam.unload("RJSMC", libpath)
}
