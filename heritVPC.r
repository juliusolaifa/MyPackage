herit.vpc <- function(x, ...) {
  UseMethod("my_function", x)
}
herit.vpc.default <- function(x, ...) {
  stop("Unsupported input type")
}
herit.vpc.fitted_model <- function(x, ...) {
  # Extract parameters from the fitted model object
  params <- extract_params(x)
  # Call the function with the parameters
  my_function.params(params, ...)
}

herit.vpc.fitted_model <- function(x, ...) {
  # Extract parameters from the fitted model object
  params <- extract_params(x)
  # Call the function with the parameters
  my_function.params(params, ...)
}

fitted_model_obj <- structure(list(), class = "fitted_model")
params_obj <- structure(list(), class = "params")