# IAM 2.0.0

* Added a NEWS.md file to track changes to the package.

* Add new functions *summary()* for in iamInput and iamArgs class.

* Add new functions to modify iamArgs class, with format *IAM.editArgs_X()*, X being "Eco", "Scenar", "Gest" or "Rep".

* Add reproductible example dataset for documentation purpose.

## Breaking change

* IAM.format is totally modified inf argument and computation. 

* IAM.args change from tcltk to shiny for the GUI

## Code modifications

* Price flexibility is extracted as a function and is activated only if the input excell file contain the price flexibility sheet.

* Replaced Gestion module after this module was deleted in earlier versions.

* Removed deprecated arguments in iamInput and iamArgs.

## License

* GPL-2 to CeCILL-2.1

## Organisation

* Documentation is now set with `Roxygen2` package format.

* Added a git repository to track down every change

* Added a `NEWS.md` file to track changes to the package.

* Split single C++ file in differents files. Each file is now a function or more related to the module it's used in.

* C++ function not related with BioEcoPar class are in a specific file.

* Notice is written in LaTeX and include vignettes compiled as pdf.


# IAM 1.7.0.9999

## **Warning** Prior versions.

Prior versions of IAM where specially developed for project usage and no NEWS.md file or equivalent was used.
Information below was compiled from code reading and adding previous versions on git.
Prior versions are archived **for reading purpose only** and using them could result in unpredictable results.

Prior versions were not using the same package name and version. Maxime edited the descriptions to build back the git historic. A summary table with original source location, version name and number is accessible for the dev team in gitlab repository **IAM/good_practice**

# IAM 1.7.0

Source communicated by Mathieu

* R :
  * Minor modification of the code, mostly argument for C++ and export of TACbyF and TACtot variables

* C++ :
  * Add new modules : EstimationTACfromF and SRmod2


# IAM 1.6.3

Source found in IAM/5 versions : test3

* Definitive choice of the ``openxlsx`` package for reading input file
* Add possibility to modify Stock Recrutement parameters during a simulation


# IAM 1.6.2

Source found in IAM/5 version : test2

Mainly modification of the GestionF2 function and related variables in R and C++

Personnal modification : added missing R/1_input_integration.r file to complete the package


# IAM 1.6.1

Source found in IAM/5 versions : test

Possibility to choose between ``XLConnetct`` and ``openxlsx`` package for reading the input file.

Intern modification of the C++ file


# IAM 1.6.0

Source found in IAM/5 versions

Dependency modification rom ``XLConnect`` to ``openxlsx``

Add over quotat discard computation


# IAM 1.5.1

Source found in IAM/5 versions

Add over quota discard


# IAM 1.5.0

Source found in IAM/5 versions

Add capacity to input **Spict species**
Mainly C++ modification


# IAM 1.4.3

Source found in IAM/5 versions

Minor modification of C++ mainly in GestionF2 function


# IAM 1.4.2

Source found in IAM/5 versions

Minor fix to prevent cases if 0 dynamic species in input and NA names.


# IAM 1.4.1

Source found in IAM/5 versions

Add updateE, which specify if the update of gestion is referenced from starting point of the simulation or the previous point.


# IAM 1.4.0

Source found in IAM/5 versions

Modification of TACbyF and set of this parameter in IAM.args() function

Add Management for static species

Landing obligation is now taken into account in model


# IAM 1.3.1

Source found in IAM versions SaisineRejets_mars2015

Add mortality for SS3 species


# IAM 1.3.0

Source found in IAM versions SaisinRejets_mars2015

Major update with SS3 and Static species added in R and C++
add effSup in Gestion


# IAM 1.2.0

Source found in IAM version ICES.
Special request.

Add smooth-HS SR model and modification of SR C++
New version of scenario function
Add TypeGest output
Default value for TACcontrol

Personnal comment : I'm not sure the parent version is 1.1.2, maybe a more ancient one (because of stochprice return).


# IAM 1.1.3

Source found in SOCIOEC directory

Fix in C++ for GestionF
Add lambda for TACbyF


# IAM 1.1.2

Source found in IAM/5 versions

* R :
  * Minor fix in input and IAM.model.

Delete TAC control options and P output.


# IAM 1.1.1

Source found in personnal computer (Mathieu)

Added stochasticity in prices.

* R :
  * Minor fix in R


# IAM 1.1.0

Source found in personnal computer (Mathieu)

* R :
  * Add Folder Fleet possibility
  * Replace xlsReadWrite for XLConnect
  * Add TACbyF argument into IAM.model function

* C++ :
  * Add functions: GestionF2, QuotaExchV2 and more.
  * Big refactoring of the C++ file.


# IAM 1.0.0

Source found in GT_BioEco

* R :
  * This version xlsReadWrite and use FLR format for the SIAD website.


# IAM 0.4.0

Source found in personnal computer (Mathieu)

* R :
  * Modify input functions.
  * Replace dependency from xlsReadWrite to XLConnect


# IAM v0.3.0

Source found in personnal computer (Mathieu)

* R :
  * Add ventilation process of input in R and not in C++ anymore.
  * Add IAMtoFLR methods

* C++ :
  * Add Normal noise in SR
  * Add mOTH and M_fm for mortality. Management possibility at different levels from global to fleet.

* Fix :
  * Now mortality and mu multiplicator in Management is minored to 0.


# IAM 0.2.1

Source found in SIAD project

* R :
  * Small fix in git plot default values

* C++ :
  * Modification of debug file output in C++ and a single formula.

Data was ignored cause not secret.


# IAM 0.2.0

Source found in SIAD project

Added option in argum to select outputRep variables in a specified list
Added plot methods.


# IAM 0.1.0

First source found in SIAD project.
