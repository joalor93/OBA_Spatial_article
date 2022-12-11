OBJECTIVE BAYESIAN ANALYSIS FOR GEOSTATISTICAL STUDENT-T PROCESSES

MAIN AUTHOR: JOSE ALEJANDRO ORDONEZ

email: ordonezjosealejandro@gmail.com

University: State University of Campinas.

These are the instructions for running the codes attached to the folder.

The scripts are organized in two internal folders:

*model_fitting *summarized_results

model_fitting: Contains scripts used for running the models. The sample obtained for each model were saved as .rds archives, and their names are listed below, for each kappa and prior

*** kappa = 0.3 ***

** 'newcadeia11kappa0.3.rds' contains the sample obtained under the reference prior** 'newcadeiakappa0.3vague.rds' contains the sample obtained under the vague prior ** 'newcadeiakappa0.3jefrul.rds' contains the sample obtained under the Jeffreys' prior ** 'newcadeiakappa0.3jefrul.rds' contains the sample obtained under the Jeffreys' prior ** 'newcadeiakappa0.3jefind.rds' contains the sample obtained under the  indepedent Jeffreys' prior


*** kappa = 0.5 ***

** 'newcadeia11kappa0.5.rds' contains the sample obtained under the reference prior** 'newcadeiakappa0.5vague.rds' contains the sample obtained under the vague prior ** 'newcadeiakappa0.5jefrul.rds' contains the sample obtained under the Jeffreys' prior ** 'newcadeiakappa0.5jefrul.rds' contains the sample obtained under the Jeffreys' prior ** 'newcadeiakappa0.5jefind.rds' contains the sample obtained under the  indepedent Jeffreys' prior


*** kappa = 0.7 ***

** 'newcadeia11kappa0.7.rds' contains the sample obtained under the reference prior** 'newcadeiakappa0.7vague.rds' contains the sample obtained under the vague prior ** 'newcadeiakappa0.7jefrul.rds' contains the sample obtained under the Jeffreys' prior ** 'newcadeiakappa0.7jefrul.rds' contains the sample obtained under the Jeffreys' prior ** 'newcadeiakappa0.7jefind.rds' contains the sample obtained under the  indepedent Jeffreys' prior


summarized results: contains the .rds archives as well as the scripts used for processing the saved samples (.rds archives), its content are listed below.


** 'summarized_results' was the script used for processing the results for all the models under different priors.

** .rds archives generated during the model fitting stage (described above).

** auxiliar functions_deviance.R: Auxiliar functions used for the computation of different model deviance selection criteria.