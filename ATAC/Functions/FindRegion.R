# convert region argument to genomic coordinates
# region can be a string, name of a gene, or GRanges object
FindRegion <- function(
  object,
  region,
  sep = c("-", "-"),
  assay = NULL,
  extend.upstream = 0,
  extend.downstream = 0
) {
  if (!is(object = region, class2 = "GRanges")) {
    # if separators are present in the string and we can convert the
    # start to a number, assume we're using genomic coordinates
    if (all(sapply(X = sep, FUN = grepl, x = region))) {
      region <- StringToGRanges(regions = region, sep = sep)
    } else {
      region <- LookupGeneCoords(object = object, assay = assay, gene = region)
      if (is.null(x = region)) {
        stop("Gene not found")
      }
    }
  }
  region <- suppressWarnings(expr = Extend(
    x = region,
    upstream = extend.upstream,
    downstream = extend.downstream
  )
  )
  return(region)
}