## Resubmission

This resubmission addresses the incoming checks:

* DESCRIPTION spelling: added full expression of acronym (e.g., “EHR”).
* HTML manual warnings: restructured roxygen details in
  - simulate_genetic_disease_from_topic
  - simulate_topics
  to avoid lists immediately following headings; the HTML validator warnings are gone.
* Examples runtime: kept small, always-run examples and gated longer demonstrations with NOT_CRAN for the following functions:
  - loading2weights
  - wrapper_ATM
  - prediction_OR
  Examples now complete well under the recommended runtime on CRAN machines.

No API changes; only documentation and example adjustments.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* The package uses a mapping table which is ~3MB to map beteewn differnt disease coding systems. We decided it is too complicated to ask the users to 
download this table from the NHS website, which require them to register a new account. We decided it is much more
user friendly to include the table as an internal data object.
