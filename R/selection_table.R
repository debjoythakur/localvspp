## Category keywords in display order; within a category, rows follow the
## keyword order and then the radius (50, 100, 200, 500, 1000 m).
default_categories <- function() {
  list(
    "Road Types"           = c("intersection", "motorway", "trunk", "primary", "residential",
                               "secondary", "tertiary", "service", "road|street|highway|junction"),
    "Buildings / Land Use" = c("building", "commercial", "industrial", "retail", "\\bbars?\\b",
                               "\\bpubs?\\b", "nightclub", "restaurant", "fast food", "cafe",
                               "cinema", "convenience", "supermarket", "fuel"),
    "Transportation"       = c("\\bbus", "parking", "rail|tram|subway|transit"),
    "Education"            = c("school", "universit|college"),
    "Healthcare"           = c("hospital", "pharmac|clinic|doctor|dentist"),
    "Recreation"           = c("\\bparks?\\b", "playground"))
}

radius_of <- function(v) {
  suppressWarnings(as.integer(ifelse(grepl("_[0-9]+m$", v), sub(".*_([0-9]+)m$", "\\1", v), NA)))
}
stem_of <- function(v) sub("_[0-9]+m$", "", v)

#' Variable-selection table across methods
#'
#' One row per covariate selected by at least one method, grouped into
#' categories and sorted by radius, with a logical column per method (the
#' layout of Table S1). A covariate counts as selected by a local method if
#' any of its Haar atoms has a non-zero coefficient, as in
#' [selected_covariates()].
#'
#' @param object A `localvs` fit.
#' @param categories Named list of regular expressions (matched against the
#'   covariate stem with underscores replaced by spaces) defining the
#'   categories and the row order within each.
#' @param labels Optional named character vector mapping covariate stems to
#'   display names, e.g. `c(fuel = "Fuel station", bus_stops = "Bus stops")`.
#' @return A data frame with columns `Category`, `Predictor`, `variable` and
#'   one logical column per fitted method. Covariates that match no category
#'   are put in `"Other"`.
#' @export
selection_table <- function(object, categories = default_categories(), labels = character(0)) {
  sel <- selected_covariates(object)
  all_sel <- sort(unique(unlist(sel)))
  classify <- function(stem) {
    s <- tolower(gsub("_", " ", stem))
    for (ci in seq_along(categories)) for (ki in seq_along(categories[[ci]]))
      if (grepl(categories[[ci]][ki], s, perl = TRUE)) return(c(ci, ki))
    c(length(categories) + 1, 0)
  }
  make_label <- function(stem, radius) {
    s <- ifelse(stem %in% names(labels), labels[stem], gsub("_", " ", stem))
    s <- paste0(toupper(substr(s, 1, 1)), substring(s, 2))
    ifelse(is.na(radius), s, paste0(s, " within ", radius, "m"))
  }
  stem <- stem_of(all_sel); radius <- radius_of(all_sel)
  cls  <- matrix(vapply(stem, classify, numeric(2)), ncol = 2, byrow = TRUE)
  out  <- data.frame(Category = c(names(categories), "Other")[cls[, 1]],
                     Predictor = unname(make_label(stem, radius)), variable = all_sel,
                     stringsAsFactors = FALSE)
  for (m in names(sel)) out[[m]] <- out$variable %in% sel[[m]]
  out <- out[order(cls[, 1], cls[, 2], stem, radius), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' LaTeX version of the variable-selection table
#'
#' Writes the output of [selection_table()] as a `table*` environment with
#' check marks and crosses, in the layout of Table S1. Needs the LaTeX
#' packages `booktabs`, `multirow` and `amssymb`.
#'
#' @param tab Output of [selection_table()].
#' @param file Optional file to write to.
#' @param caption,label Table caption and LaTeX label.
#' @return Invisibly, the LaTeX lines.
#' @export
selection_latex <- function(tab, file = NULL,
                            caption = paste0("Comparison of selected variables. A check mark ($\\checkmark$) ",
                                             "indicates that the corresponding predictor was selected by the method, ",
                                             "whereas ($\\times$) indicates that it was not selected."),
                            label = "tab:variable_selection") {
  meth <- setdiff(names(tab), c("Category", "Predictor", "variable"))
  ck <- function(x) ifelse(x, "$\\checkmark$", "$\\times$")
  blocks <- lapply(split(tab, factor(tab$Category, levels = unique(tab$Category))), function(g) {
    marks <- apply(as.matrix(as.data.frame(lapply(g[meth], ck))), 1, paste, collapse = " & ")
    c(sprintf("\\multirow{%d}{*}{%s}", nrow(g), g$Category[1]), "",
      paste0("& ", g$Predictor, "\n& ", marks, " \\\\\n"), "\\midrule", "")
  })
  body <- unlist(blocks, use.names = FALSE)
  body <- body[seq_len(length(body) - 2)]            # drop the trailing \midrule
  out <- c("\\begin{table*}[!htbp]", "\\centering", paste0("\\caption{", caption, "}"),
           paste0("\\label{", label, "}"), "",
           "\\renewcommand{\\arraystretch}{1.00}", "\\setlength{\\tabcolsep}{3pt}", "",
           "\\resizebox{0.6\\textwidth}{!}{%", paste0("\\begin{tabular}{ll", strrep("c", length(meth)), "}"),
           "\\toprule",
           paste0("\\textbf{Category} & \\textbf{Predictor} & ",
                  paste0("\\textbf{", meth, "}", collapse = " & "), " \\\\"),
           "\\midrule", "", body, "\\bottomrule", "\\end{tabular}%", "}", "\\end{table*}")
  if (!is.null(file)) writeLines(out, file)
  invisible(out)
}

#' Read a variable-selection table back from LaTeX
#'
#' Parses a table in the layout written by [selection_latex()] (for example,
#' the table pasted from the manuscript) so it can be checked with
#' [compare_selection()].
#'
#' @param path Path to a `.tex` file.
#' @param methods Method codes of the last columns, in order.
#' @return A data frame with `Predictor` and one logical column per method.
#' @export
read_selection_latex <- function(path, methods = names(method_labels)) {
  txt  <- paste(readLines(path, warn = FALSE), collapse = "\n")
  txt  <- gsub("(?<!\\\\)%[^\n]*", "", txt, perl = TRUE)          # strip LaTeX comments
  rows <- strsplit(txt, "\\\\", fixed = TRUE)[[1]]                  # split on row ends
  rows <- gsub("\\\\multirow\\{[^}]*\\}\\{[^}]*\\}\\{[^}]*\\}", "", rows)
  k <- length(methods)
  out <- lapply(rows, function(r) {
    f <- trimws(strsplit(r, "&", fixed = TRUE)[[1]])
    if (length(f) < k + 1) return(NULL)
    mk <- utils::tail(f, k)
    if (!all(grepl("checkmark|times", mk))) return(NULL)
    d <- data.frame(Predictor = f[length(f) - k], stringsAsFactors = FALSE)
    d[methods] <- as.list(grepl("checkmark", mk))
    d
  })
  out <- do.call(rbind, out)
  if (is.null(out)) stop("No table rows found in ", path, call. = FALSE)
  out
}

#' Compare two variable-selection tables
#'
#' Matches rows by predictor label (case and spacing ignored) and reports,
#' per method, the number of selected covariates in each table and the
#' cells that disagree.
#'
#' @param code_tab Output of [selection_table()].
#' @param paper_tab Output of [read_selection_latex()].
#' @param verbose Print the comparison?
#' @return Invisibly, a list with `summary`, `diffs` (cell-level mismatches),
#'   `only_code` and `only_paper` (rows present in only one table).
#' @export
compare_selection <- function(code_tab, paper_tab, verbose = TRUE) {
  meth <- intersect(setdiff(names(code_tab), c("Category", "Predictor", "variable")), names(paper_tab))
  key  <- function(x) tolower(gsub("\\s+", " ", trimws(x)))
  a <- data.frame(key = key(code_tab$Predictor),  code_tab[meth],  check.names = FALSE)
  b <- data.frame(key = key(paper_tab$Predictor), paper_tab[meth], check.names = FALSE)
  keys <- union(a$key, b$key)
  get  <- function(d, m) { v <- d[[m]][match(keys, d$key)]; v[is.na(v)] <- FALSE; v }
  diffs <- do.call(rbind, lapply(meth, function(m) {
    cc <- get(a, m); pp <- get(b, m); i <- which(cc != pp)
    data.frame(Predictor = keys[i], Method = rep(m, length(i)), Code = cc[i], Paper = pp[i],
               stringsAsFactors = FALSE)
  }))
  summ <- data.frame(Method = meth, Code_selected = colSums(a[meth]), Paper_selected = colSums(b[meth]),
                     Mismatches = vapply(meth, function(m) sum(diffs$Method == m), numeric(1)),
                     row.names = NULL)
  only_code <- setdiff(a$key, b$key); only_paper <- setdiff(b$key, a$key)
  if (verbose) {
    cat("==== Code vs. manuscript table ====\n")
    print(summ, row.names = FALSE)
    cat("\nRows only in code table :", if (length(only_code))  paste(only_code,  collapse = "; ") else "none", "\n")
    cat("Rows only in paper table:", if (length(only_paper)) paste(only_paper, collapse = "; ") else "none", "\n")
    if (nrow(diffs)) { cat("\nCell-level mismatches:\n"); print(diffs, row.names = FALSE) }
    else cat("\nAll check marks agree.\n")
  }
  invisible(list(summary = summ, diffs = diffs, only_code = only_code, only_paper = only_paper))
}
