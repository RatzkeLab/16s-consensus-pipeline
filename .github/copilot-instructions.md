
this is a snakemake 9 pipeline
use snakemake 9 best practices
use snakemake 9 syntax
avoid lambdas in rule inputs, rather put functions in common.smk
i am currently totally overhauling the codebase, as it does not work in its current form. don't worry about backward compatibility.
rules are found in workflow/rules
scripts (including python) are found in workflow/scripts
workflow/Snakefile is the main snakemake file
before running anything with snakemake run this command separately conda activate snakemake_env