#' 26 double-blind, randomized, active, or placebo-controlled clinical trials on patients with primary hypercholesterolemia sponsored by Merck & Co., Inc., Kenilworth, NJ, USA.
#' 
#' A data set containing clinical trial on hypercholesterolemia including 26 trials and 2 treatment arms each, and other attributes of the participants
#' 
#' @format A data frame with 52 rows and 19 variables
#' \describe{
#'		\item{Study}{study identifier}
#'		\item{Trial}{trial identifier}
#'		\item{Npt}{the number of participants per trial}
#'		\item{pldlc}{mean percentage difference in LDL-C}
#'		\item{phdlc}{mean percentage difference in HDL-C}
#'		\item{ptg}{mean percentage difference in triglycerides (TG)}
#'		\item{sdldl}{sample standard deviation of percentage difference in LDL-C}
#'		\item{sdhdl}{sample standard deviation of percentage difference in HDL-C}
#'		\item{sdtg}{sample standard deviation of percentage difference in triglycerides (TG)}
#'		\item{onstat}{whether the participants were on Statin prior to the trial}
#'		\item{trt}{treatment indicator for Statin or Statin+Ezetimibe}
#'		\item{bldlc}{baseline LDL-C}
#'		\item{bhdlc}{baseline HDL-C}
#'		\item{btg}{baseline triglycerides (TG)}
#'		\item{age}{age in years}
#'		\item{white}{the proportion of white participants}
#'		\item{male}{the proportion of male participants}
#'		\item{dm}{the proportion of participants with diabetes mellitus}
#'		\item{durat}{duration in weeks}
#' }
#' @usage data(cholesterol)
#' @keywords datasets
#' @examples
#' data(cholesterol)
"cholesterol"