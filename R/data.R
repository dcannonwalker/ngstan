#' Translational Efficiency Ribo-Seq and RNA-Seq data
#'
#' Ribo-Seq and RNA-Seq counts per gene from
#' Report ...
#'
#' @format ## `dte`
#' A data frame with 27,628 rows and 37 columns:
#' \describe{
#'   \item{geneid}{Gene identifier}
#'   \item{ribo_..., rna_...}{Ribo-Seq and RNA-Seq counts for treatment
#'   and control samples at 5 and 8 days post infection}
#' }
#' @source <https://www.ncbi.nlm.nih.gov/bioproject/PRJNA950066>
"dte"

#' Gene Expression for MOV10 data
#'
#' Counts per gene
#'
#' @format ## `mov10`
#' A data frame with 14,675 rows and 9 columns:
#' \describe{
#'  \item{geneid}{Gene identifier}
#'  \item{low_..., high_..., control_...}{Gene expression counts
#'  for each sample}
#' }
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50499>
"mov10"

#' Tuch tumor data
#'
#' Counts per gene
#'
#' @format ## `tuch_tumor`
#' A data frame with 15,668 rows and 9 columns:
#' \describe{
#'  \item{id_ref_seq}{Gene identifier}
#'  \item{name_of_gene}{Gene name}
#'  \item{number_of_exons}{Exons per gene}
#'  \item{8n, 8t, 33n, 33t, 51n, 51t}{Gene expression counts
#'  for each sample}
#' }
#' @source <https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009317#s5>
"tuch_tumor"
