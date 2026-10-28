## Resubmission

This is a resubmission. In this version I have:

* Removed "Tools for" from the start of the Description field.
* Added the reference in the Description field in the requested form,
  Allen and Allen (2013, ISBN:978-0-470-67377-5).
* Changed plot_event_chart() so that it saves all graphical parameters with
  `oldpar <- par(no.readonly = TRUE)` and restores them with an immediate
  `on.exit(par(oldpar))`. plot_basin_model() already did this. No function
  changes the user's options or working directory.

Following the CRAN checklist I also shortened the Title to 58 characters
("Burial and Thermal History Modelling of Sedimentary Basins"), added the
copyright holder role [cph] to Authors@R, corrected the full title of
Allen and Allen (2013) in the documentation and added the CRAN installation
line to the README.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.
* Possibly misspelled words in DESCRIPTION: geotherm, radiogenic,
  stratigraphers, stratigraphic. These are correct geoscience terms.
