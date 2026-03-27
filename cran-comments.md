## Resubmission

This resubmission addresses the incoming checks:

* A CRAN package should not be larger than 5 MB. Please reduce the size. ->

*If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year, ISBN:...)

Added: Jiang (2023) <doi:https://doi.org/10.1038/s41588-023-01522-8>

*Please write TRUE and FALSE instead of T and F. Please don't use "T" or
"F" as vector names.

Updated wrapper_ATM and wrapper_LFA. 


## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* The package uses a mapping table which is ~3MB to map between differnt disease coding systems. We decided it is too complicated to ask the users to 
download this table from the NHS website, which require them to register a new account. We decided it is much more
user friendly to include the table as an internal data object.
