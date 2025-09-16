#' Example HES diagnosis ages
#'
#' A tiny illustrative subset of Hospital Episode Statistics (HES) data
#' with participant IDs and ages at diagnosis, used in examples and tests.
#'
#' @format A data frame/tibble with example rows. Typical columns include:
#' \describe{
#'   \item{eid}{Participant identifier (integer or character).}
#'   \item{age_diag}{Age at diagnosis (numeric).}
#'   \item{diag_icd10}{ICD-10 diagnosis code (character).}
#' }
#' @examples
#' head(HES_age_example)
"HES_age_example"

#' Example HES ICD-10 diagnoses
#'
#' A tiny illustrative subset of HES diagnoses with participant IDs and ICD-10 codes.
#'
#' @format A data frame/tibble with example rows. Typical columns include:
#' \describe{
#'   \item{eid}{Participant identifier (integer or character).}
#'   \item{diag_icd10}{ICD-10 diagnosis code (character).}
#'   \item{age_diag}{ICD-10 diagnosis age point (double).}
#' }
#' @examples
#' head(HES_icd10_example)
"HES_icd10_example"

#' SNOMED ↔ ICD-10(-CM) mapping (excerpt)
#'
#' A small mapping table used by functions such as \code{\link{icd2phecode}}
#'
#' @format A data frame/tibble. Common columns include:
#' \describe{
#'   \item{SNOMED}{SNOMED CT concept identifier (character).}
#'   \item{ICD10}{ICD-10 code (character), and/or}
#'   \item{ICD10_name}{ICD-10-CM code (character).}
#'   \item{SNOMED_description}{SNOMED readable explanation}
#'   \item{occ}{ICD10 occurence in UKB}
#' }
#' @examples
#' head(SNOMED_ICD10CM)
"SNOMED_ICD10CM"

#' List of 349 UK Biobank diseases (example)
#'
#' A character vector or table listing the set of disease phenotypes used
#' in examples/vignettes.
#'
#' @format A data frame/tibble containing disease identifiers/names. Columns include:
#' \describe{
#'   \item{diag_icd10}{Phecode (character).}
#'   \item{occ}{number of distinct patient in UKB}
#' } @examples
#' head(UKB_349_disease)
"UKB_349_disease"

#' Example topic model output (10 topics, UKB HES)
#'
#' An illustrative result object/table from a 10-topic model fit to UKB HES-like data;
#' used for examples, plotting, and tests.
#'
#' @format An array for UKB topic loadings.
#' Dimention is age, disease, topics. the ordering of disease is the same as UKB_349_disease.
#' @examples
#' head(UKB_HES_10topics)
"UKB_HES_10topics"

#' Disease information linking PheCodes and ICD-10
#'
#' A helper table with disease metadata to support mapping between PheCodes
#' and ICD-10.
#'
#' @format A data frame/tibble. Common columns include:
#' \describe{
#'   \item{phecode}{PheCode as character.}
#'   \item{ICD10}{Alternative PheCode column name (if present).}
#'   \item{exclude_range}{ancestor PheCode range (character).}
#'   \item{phenotype}{Human-readable phenotype/label (character), if available.}
#'  \item{exclude_name}{ancestor PheCode name (character).}
#' }
#' @examples
#' head(disease_info_phecode_icd10)
"disease_info_phecode_icd10"

#' ICD-10 ↔ PheCode mapping
#'
#' Mapping table between ICD-10 codes and PheCodes.
#'
#' @format A data frame/tibble. Common columns include:
#' \describe{
#'   \item{ICD10}{ICD-10 code (character).}
#'   \item{PheCode}{PheCode (character).}
#'   \item{Excl..Phecodes}{ancestor PheCode range (character).}
#'   \item{Excl..Phenotypes}{ancestor PheCode name (character).}
#' }
#' @examples
#' head(phecode_icd10)
"phecode_icd10"

#' ICD-10-CM ↔ PheCode mapping
#'
#' Mapping table between ICD-10-CM codes and PheCodes.
#'
#' @format A data frame/tibble. Common columns include:
#' \describe{
#'   \item{ICD10}{ICD-10-CM code (character).}
#'   \item{phecode}{PheCode (character).}
#'   \item{exclude_range}{ancestor PheCode range (character).}
#'   \item{exclude_name}{ancestor PheCode name (character).}
#' }
#' @examples
#' head(phecode_icd10cm)
"phecode_icd10cm"

#' Short labels (at most first for letters/digits) for ICD-10 codes
#'
#' A lookup table mapping ICD-10 codes to concise human-readable labels.
#'
#' @format A data frame/tibble. Common columns include:
#' \describe{
#'   \item{ICD10}{ICD-10 code (character).}
#'   \item{parent_phecode}{phecode of parent node (character).}
#'   \item{Excl..Phecodes}{ancestor PheCode range (character).}
#'   \item{Excl..Phenotypes}{ancestor PheCode name (character).}
#'   \item{occ}{number of distinct patient in UKB}
#' }
#' @examples
#' head(short_icd10)
"short_icd10"

#' Short labels (at most first for letters/digits) for ICD-10-CM codes
#'
#' A lookup table mapping ICD-10-CM codes to concise human-readable labels.
#'
#' @format A data frame/tibble. Common columns include:
#' \describe{
#'   \item{ICD10}{ICD-10 code (character).}
#'   \item{parent_phecode}{phecode of parent node (character).}
#'   \item{exclude_range}{ancestor PheCode range (character).}
#'   \item{exclude_name}{ancestor PheCode name (character).}
#'   \item{occ}{number of distinct patient in UKB}
#' }
#' @examples
#' head(short_icd10cm)
"short_icd10cm"
