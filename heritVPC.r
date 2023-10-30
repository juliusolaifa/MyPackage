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
  herit.vpc.params(params, ...)
}

herit.vpc.params <- function(x, ...) {
  # Implement the functionality using the parameters
  # x is expected to be a list or similar structure containing the parameters
  ...
}

extract_params <- function(x) {
}

fitted_model_obj <- structure(list(), class = "fitted_model")
params_obj <- structure(list(), class = "params")
