#### sacurine ####

#' Analysis of the human adult urinary metabolome variations with age, body
#' mass index and gender
#'
#' Urine samples from 183 human adults were analyzed by liquid chromatography
#' coupled to high-resolution mass spectrometry (LTQ Orbitrap) in the negative
#' ionization mode. A total of 109 metabolites were identified or annotated at
#' the MSI level 1 or 2. After retention time alignment with XCMS, peaks were
#' integrated with Quan Browser. After signal drift and batch effect correction
#' of intensities, each urine profile was normalized to the osmolality of the
#' sample. Finally, the data were log10 transformed.
#'
#' @name sacurine
#' @docType data
#' @format A list with the following elements:
#' \itemize{
#' \item dataMatrix a 183 samples x
#' 109 variables matrix of numeric type corresponding to the intensity profiles
#' (values have been log10-transformed)
#' \item sampleMetadata a 183 x 3 data frame, with the volunteers' age
#' ('age', numeric), body mass index ('bmi',
#' numeric), and gender ('gender', factor)
#' \item variableMetadata a 109 x 3 data frame, with the metabolites'
#' MSI identification level ('msiLevel':
#' either 1 or 2), HMDB ID when available ('hmdb', character), chemical class
#' according to the 'super class' taxonomy of HMDB ('chemicalClass', character)
#' }
#' @return List containing the 'dataMatrix' matrix (numeric) of data (samples
#' as rows, variables as columns), the 'sampleMetadata' data frame of sample
#' metadata, and the variableMetadata data frame of variable metadata. Row
#' names of 'dataMatrix' and 'sampleMetadata' are identical. Column names of
#' 'dataMatrix' are identical to row names of 'variableMetadata'. For details
#' see the 'Format' section above.
#' @references Thevenot E.A., Roux A., Xu Y., Ezan E. and Junot C. (2015).
#' Analysis of the human adult urinary metabolome variations with age, body
#' mass index and gender by implementing a comprehensive workflow for
#' univariate and OPLS statistical analyses. Journal of Proteome Research, DOI:
#' 10.1021/acs.jproteome.5b00354
#' @keywords datasets
NULL


#### foods ####

#' Food consumption patterns accross European countries (FOODS)
#'
#' The relative consumption of 20 food items was compiled for 16 countries. The
#' values range between 0 and 100 percent and a high value corresponds to a
#' high consumption. The dataset contains 3 missing data.
#'
#'
#' @name foods
#' @docType data
#' @format A data frame with the following parameters:
#' \itemize{
#' \item Country Name of the country
#' \item Gr_CoffeGround Coffee
#' \item Inst_Coffe Instant Coffee
#' \item Tea Tea \item Sweetner Sweetner
#' \item Biscuits Biscuits \item Pa_Soup Powder Soup
#' \item Ti_Soup Tin Soup \item In_Potat Instant Potatoes
#' \item Fro_Fish Frozen Fish
#' \item Fro_Veg Frozen Vegetables
#' \item Apples Apples
#' \item Oranges Oranges
#' \item Ti_Fruit Tin Fruit
#' \item Jam Jam
#' \item Garlic Garlic
#' \item Butter Butter
#' \item Margarine Margarine
#' \item Olive_Oil Olive Oil
#' \item Yoghurt Yoghurt
#' \item Crisp_Brea Crisp Bread
#' }
#' @return Data frame (numeric type except the first column, which can be
#' transformed into row names) with 16 rows and 21 columns, corresponding to
#' the 'Country' (column 1), followed by the consumption of each of the 20 food
#' items (columns 2-21). For details see the 'Format' section above.
#' @references Eriksson et al. (2006). Multi- and Megarvariate Data Analysis.
#' Umetrics Academy. pp.10, 33, 48.
#' @keywords datasets
NULL


#### lowarp ####

#' A multi response optimization data set (LOWARP)
#'
#' This example concerns the development of a polymer similar to that used in
#' the plastic covering of mobile phones. The desired profile of the polymer
#' was low warp and high strength. Four constituents (glas, crtp, mica, and
#' amtp) were varied in the polymer formulation by means of a 17 run mixture
#' design. For each new polymer, i.e., each new experiment in the mixture
#' design, 14 responses relating to both warp and strength were measured on the
#' product. The objective of the data analysis was to uncover which combination
#' of factors (the four ingredients) gave polymers with low warp and high
#' strength. The data set contains 10 missing values (NA).
#'
#'
#' @name lowarp
#' @docType data
#' @format A data frame with the following parameters:
#' \itemize{
#' \item num mixture number
#' \item glas glas constituent
#' \item crtp crtp constituent
#' \item mica mica constituent
#' \item amtp amtp constituent
#' \item wrp1 warp response 1
#' \item wrp2 warp response 2
#' \item wrp3 warp response 3
#' \item wrp4 warp response 4
#' \item wrp5 warp response 5
#' \item wrp6 warp response 6
#' \item wrp7 warp response 7
#' \item wrp8 warp response 8
#' \item st1 strength response 1
#' \item st2 strength response 2
#' \item st3 strength response 3
#' \item st4 strength response 4
#' \item st5 strength response 5
#' \item st6 strength response 6
#' }
#' @return Data frame (numeric type only; the first column can be transformed
#' into row names) with 17 rows and 19 columns corresponding to the subject's
#' 'num'ber (column 1), the 4 constituent variables (columns 2-5), the 8 warp
#' responses (columns 6-13), and the 6 strength responses (columns 14-19). For
#' details see the 'Format' section above.
#' @references Eriksson et al. (2006). Multi- and Megarvariate Data Analysis.
#' Umetrics Academy. pp.16, 77, 209.
#' @keywords datasets
NULL


#### cornell ####

#' Octane of various blends of gasoline
#'
#' Twelve mixture component proportions of the blend are analysed
#'
#'
#' @name cornell
#' @docType data
#' @format A data frame with the following parameters:
#' \itemize{
#' \item num mixture number
#' \item x1 proportion of
#' component 1
#' \item x2 proportion of component 2
#' \item x3 proportion of component 3
#' \item x4 proportion of component 4
#' \item x5 proportion of component 5
#' \item x6 proportion of component 6
#' \item x7 proportion of component 7 Note: the 7 variables are correlated since they sum up to 1
#' \item y octane (quantitative variable)
#' }
#' @return Data frame (numeric type only; the first column can be transformed
#' into row names) with 12 rows and 9 columns corresponding to the 'num'ber of
#' the mixture (column 1), the proportion of each of the 7 'x' components
#' within the mixture (columns 2-8), and the octane indice 'y' (column 9). For
#' details see the 'Format' section above.
#' @references Tenenhaus (1998). La regression PLS: theorie et pratique. Paris:
#' Editions Technip.
#' @source Tenenhaus (1998), Table 6, page 78.
#' @keywords datasets
NULL