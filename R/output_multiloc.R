#' Format multi-location vaccine model output
#'
#' Takes raw odin vaccine model output from \code{run()} (with multiple
#' locations) and returns a long data frame. Selected compartments are extracted
#' (with cumulative outputs converted to daily increments) while optionally
#' collapsing the age dimension. Location indices are preserved in the output.
#'
#' @param x nimue_simulation object returned by \code{run()}.
#' @param compartments Vector of compartment names present in the model output,
#'   e.g. \code{c("S", "R")}, or sub-compartment names (e.g.
#'   \code{c("S", "E1", "E2")}).
#' @param reduce_age Logical; if \code{TRUE} (default) collapse the age
#'   dimension while retaining location indices, otherwise keep age-specific
#'   counts.
#' @param date_0 Optional Date for time 0 (e.g. "2020-03-01"); if provided a
#'   \code{date} column is added.
#' @param replicate Replicate to format. Defaults to 1.
#'
#' @return A data.frame in long format containing the selected compartments,
#'   time, replicate, and location (and optionally age group) labels.
#' @export
format_multiloc <- function(x,
                   compartments = c("S", "E",
                                    "IMild", "ICase", "IHosp",
                                    "R", "D"),
                   reduce_age = TRUE,
                   date_0 = NULL,
                   replicate = 1){

  # Arg checks
  assert_custom_class(x, "nimue_simulation")
  assert_logical(reduce_age)

  # Standardise output dimensions (drop unused dim if present)
  if (length(dim(x$output)) == 4) {
    x$output <- abind::adrop(x$output, drop = c(FALSE, FALSE, FALSE, TRUE))
  }

  # Get columns indices of variables
  index <- odin_index(x$model)
  if (!all(compartments %in% names(index))) {
    stop("Some compartments specified not output by model")
  }

  # Extract time
  time <- x$output[, index$time, replicate]

  output <- format_internal_multiloc(x = x,
                            compartments = compartments,
                            reduce_age = reduce_age,
                            index = index,
                            time = time,
                            replicate = replicate)

  # Set levels (order) of output variables
  output$compartment <- factor(output$compartment, levels = c(compartments))

  # Add date
  if (!is.null(date_0)) {
    assert_date(date_0)
    output$date <- as.Date(output$t + as.Date(date_0),
                           format = "%Y-%m-%d")
  }

  # Add age-groups if present
  if ("age_index" %in% names(output)) {
    ag <- c(paste0(seq(0, 75, 5), "-", seq(5, 80, 5)), "80+")
    output$age_group <- factor(ag[output$age_index], levels = ag)
    output <- output %>%
      dplyr::select(-.data$age_index)
  }

  # You now also have `location_index` in the output if multiple locations
  # (you can map this to country/region names outside this function).

  return(output)
}

# -------------------------------------------------------------------------
# Internals
# -------------------------------------------------------------------------

#' Internal helper for formatting multi-location output
#' @inheritParams format
#' @param index odin output index
#' @param time Time vector from the odin model output
#' @param replicate Output replicate number
format_internal_multiloc <- function(x, compartments, reduce_age, index, time,
                            replicate){

  # Convert cumulative outputs to increments
  i_convert <- unlist(index[grepl("_cumu", names(index))])
  if (length(i_convert) > 0) {
    x$output[, i_convert, replicate] <- apply(
      x$output[, i_convert, replicate, drop = FALSE],
      2,
      function(z) {
        z - dplyr::lag(z)
      }
    )

    names(index)[grepl("_cumu", names(index))] <- sapply(
      strsplit(names(index)[grepl("_cumu", names(index))], "_"),
      `[`,
      1
    )
  }

  # Select variables and summary outputs
  get <- compartments[compartments %in% names(index)]
  if (length(get) == 0L) {
    stop("No requested compartments are present in model output")
  }
  i_select <- index[get]

  # Extract outputs, collapsing over vaccine dimension where required,
  # but preserving age and location dimensions if present.
  o <- lapply(i_select, function(idx, y) {
    extract_compartment_multi_loc(idx = idx,
                                  output = y,
                                  replicate = replicate)
  }, y = x$output)

  # Collapse age or keep age dimension
  if (reduce_age) {
    o <- lapply(o, collapse_age)
  } else {
    o <- lapply(o, add_age)
  }

  # Add names of variables
  for (i in seq_along(o)) {
    o[[i]] <- data.frame(o[[i]]) %>%
      dplyr::mutate(compartment = names(o)[i])
  }

  # Add time and replicate columns
  o <- dplyr::bind_rows(o) %>%
    dplyr::mutate(
      t = rep(time, dplyr::n() / length(time)),
      replicate = replicate
    )

  return(o)
}

# -------------------------------------------------------------------------
# Helper to extract arrays from odin index with multiple locations
# -------------------------------------------------------------------------

#' @noRd
extract_compartment_multi_loc <- function(idx, output, replicate) {

  # No dim -> simple vector of indices
  if (is.null(dim(idx))) {
    return(output[, idx, replicate])
  }

  # 2D index: age x "strata" (e.g. vaccine strata OR already-aggregated)
  if (length(dim(idx)) == 2L) {
    # As in original nimue code: sum over columns of the index matrix
    # for each age row.
    return(
      apply(
        idx,
        1,
        function(a, b) {
          rowSums(b[, a, replicate])
        },
        b = output
      )
    )
  }

  # 3D index: age x strata x location
  if (length(dim(idx)) == 3L) {
    n_age   <- dim(idx)[1]
    n_loc   <- dim(idx)[3]
    n_time  <- dim(output)[1]

    res <- array(NA_real_, dim = c(n_time, n_age, n_loc))

    for (loc in seq_len(n_loc)) {

      # idx[ , , loc] is age x strata; sum across strata as before
      res[, , loc] <- apply(
        idx[, , loc],
        1,
        function(a, b) {
          rowSums(b[, a, replicate])
        },
        b = output
      )
    }

    return(res)
  }

  stop("extract_compartment_multi_loc: unsupported index dimensions")
}

# -------------------------------------------------------------------------
# Age / age+location helpers
# -------------------------------------------------------------------------

#' Keep age groups (and locations, if present)
#'
#' @param x age-disaggregated odin output matrix or array
#'   - 2D: time x age
#'   - 3D: time x age x location
#'
#' @return long matrix with columns:
#'   - age_index, value              (single location)
#'   - age_index, location_index, value (multiple locations)
add_age <- function(x){

  # Single-location case: time x age
  if (length(dim(x)) == 2L) {

    m <- matrix(
      c(
        rep(seq_len(ncol(x)), each = nrow(x)),
        as.vector(x)
      ),
      ncol = 2
    )
    colnames(m) <- c("age_index", "value")
    return(m)
  }

  # Multi-location case: time x age x location
  if (length(dim(x)) == 3L) {

    n_time <- dim(x)[1]
    n_age  <- dim(x)[2]
    n_loc  <- dim(x)[3]

    age_index <- rep(rep(seq_len(n_age), each = n_time), times = n_loc)
    location_index <- rep(seq_len(n_loc), each = n_time * n_age)

    m <- matrix(
      c(age_index, location_index, as.vector(x)),
      ncol = 3
    )
    colnames(m) <- c("age_index", "location_index", "value")
    return(m)
  }

  stop("add_age: unsupported input dimensions")
}

#' Collapse age groups (preserving locations if present)
#'
#' @param x age-disaggregated odin output matrix or array
#'   - 2D: time x age
#'   - 3D: time x age x location
#'
#' @return age-aggregated output:
#'   - 2D input: matrix with columns "value"
#'   - 3D input: matrix with columns "location_index", "value"
collapse_age <- function(x){

  # Single-location: time x age
  if (length(dim(x)) == 2L) {
    m <- matrix(rowSums(x), ncol = 1)
    colnames(m) <- "value"
    return(m)
  }

  # Multi-location: time x age x location
  if (length(dim(x)) == 3L) {
    n_time <- dim(x)[1]
    n_loc  <- dim(x)[3]

    # Sum over age dimension -> time x location
    x_sum <- apply(x, c(1, 3), sum)  # time x location

    m <- matrix(
      c(
        rep(seq_len(n_loc), each = n_time),
        as.vector(x_sum)
      ),
      ncol = 2
    )
    colnames(m) <- c("location_index", "value")
    return(m)
  }

  stop("collapse_age: unsupported input dimensions")
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
