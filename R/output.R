#' Format vaccine model output
#'
#' Takes raw odin vaccine model output from \code{run()} and returns a long
#' data frame. Selected compartments are extracted (and cumulative outputs are
#' converted to daily increments), with the option to collapse over age groups
#' and to add calendar dates.
#'
#' @param x nimue_simulation object returned by \code{run()}.
#' @param compartments Vector of compartment names present in the model output
#'   (e.g. \code{c("S", "R")}) or sub-compartment names (e.g.
#'   \code{c("S", "E1", "E2")}).
#' @param reduce_age Logical; if \code{TRUE} (default) collapse the age
#'   dimension, otherwise retain age-specific counts.
#' @param date_0 Optional Date for time 0 (e.g. "2020-03-01"); if provided a
#'   \code{date} column is added.
#' @param replicate Replicate to format. Defaults to 1.
#'
#' @return A data.frame in long format containing the selected compartments,
#'   time, replicate, and (optionally) age group labels.
#' @export
format <- function(x,
                   compartments = c("S", "E",
                                    "IMild", "ICase", "IHosp",
                                    "R", "D"),
                   reduce_age = TRUE,
                   date_0 = NULL,
                   replicate = 1){

  # Arg checks
  assert_custom_class(x, "nimue_simulation")
  assert_logical(reduce_age)

  # Standardise output dimensions
  if(length(dim(x$output)) == 4){
    x$output <- abind::adrop(x$output, drop = c(FALSE, FALSE, FALSE, TRUE))
  }

  # Get columns indices of variables
  index <- odin_index(x$model)
  if(!all(compartments %in% names(index))){
    stop("Some compartments specified not output by model")
  }

  # Extract time
  time <- x$output[,index$time, replicate]

  output <- format_internal(x = x, compartments = compartments,
                            reduce_age = reduce_age, index = index,
                            time = time, replicate = replicate)

  # Set levels (order) of output variables
  output$compartment <- factor(output$compartment, levels = c(compartments))

  # Add date
  if(!is.null(date_0)){
    assert_date(date_0)
    output$date <- as.Date(output$t + as.Date(date_0),
                           format = "%Y-%m-%d")
  }

  # Add age-groups if present
  if("age_index" %in% names(output)){
    ag <- c(paste0(seq(0, 75, 5), "-", seq(5, 80, 5)), "80+")
    output$age_group = factor(ag[output$age_index], levels = ag)
    output <- output  %>%
      dplyr::select(-.data$age_index)
  }

  return(output)
}

#' Internal helper for formatting single-location output
#' @inheritParams format
#' @param index odin output index
#' @param time Time vector from the odin model output
#' @param replicate Output replicate number
format_internal <- function(x, compartments, reduce_age, index, time,
                            replicate){

  # Convert cumulative outputs
  i_convert <- unlist(index[grepl("_cumu", names(index))])
  x$output[, i_convert, replicate] <- apply(x$output[,i_convert, replicate], 2, function(x){
    x - dplyr::lag(x)
  })
  names(index)[grepl("_cumu", names(index))] <- sapply(strsplit(names(index)[grepl("_cumu", names(index))], "_"), `[`, 1)

  # Select variables and summary outputs
  get <- c(compartments)
  get <- get[get %in% names(index)]
  i_select <- index[get]
  # Select outputs, collapsing over vaccine dimension where required
  o <- lapply(i_select, function(x, y){
    if(is.matrix(x)){
      apply(x, 1, function(a, b){
        rowSums(b[,a,replicate])
      }, b = y)
    } else {
      y[,x,replicate]
    }
  }, y = x$output)

  # Collapse age
  if(reduce_age){
    o <- lapply(o, collapse_age)
  } else {
    o <- lapply(o, add_age)
  }

  # Add names of variables
  for(i in 1:length(o)){
    o[[i]] <- data.frame(o[[i]]) %>%
      dplyr::mutate(compartment = names(o)[i])
  }

  # Add time and replicate columns
  o <- dplyr::bind_rows(o) %>%
    dplyr::mutate(t = rep(time, dplyr::n() / length(time)),
                  replicate = replicate)

  return(o)
}

#' Keep age groups
#'
#' @param x age-disaggregated odin output matrix
#'
#' @return age-disaggregated output matrix
add_age <- function(x){
  m <- matrix(c(rep(1:ncol(x), each = (nrow(x))), as.vector(x)), ncol = 2)
  colnames(m) <- c("age_index", "value")
  return(m)
}

#' Collapse age groups
#'
#' @param x age-disaggregated odin output matrix
#'
#' @return age-aggregated output matrix
collapse_age <- function(x){
  m <- matrix(rowSums(x), ncol = 1)
  colnames(m) <- "value"
  return(m)
}

## Index locations of outputs in odin model
#' @noRd
odin_index <- function(model) {
  len <- length(model$.__enclos_env__$private$ynames)
  model$transform_variables(seq_len(len))
}

#' @noRd
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

