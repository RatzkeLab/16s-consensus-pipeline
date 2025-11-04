
this is a snakemake 9 pipeline that runs per sample given by the user, aggrigating only at the very end (to help with parallelization)
use snakemake 9 best practices
use snakemake 9 syntax
avoid lambdas in rule inputs, rather put functions in common.smk
i am currently working on this project by myself and iterating quickly, so don't worry about backward compatibility.
rules are found in workflow/rules
scripts (including python) are found in workflow/scripts
workflow/Snakefile is the main snakemake file
before running anything with snakemake run this command separately conda activate snakemake_env
all pipeline variables (including paths or folder names etc) should be declared in common.smk so there is a single source of truth
most variables and parameters should be configurable from the default configfile: "demo/default_config.yaml" or from a user configfile passed to snakemake with --configfile
use conda for all software dependencies, with environments declared in workflow/envs
it is preferable to use external tools instead of writing custom scripts, such as mafft or nanofilt
if there are two ways of doing something, prefer the simpler way, but also feel free to show me both options and ask what I prefer. 
if there is unclarity about what is meant when I give you instructions, ask for clarification
for now, try to be brief and compact as possible in each rule and each function you write - use good practice, but leave out bells, whistles, and optimisations, let's just get the basic minimum viable product working first
running snakemake -c 4 --use-conda should run the pipeline on the demo data in demo/ with 4 cores, and is a good way to check that everything is working on real data
when i ask a question (in "Ask" mode rather than in "Agent" mode), i want you to give a brief answer to my question, and only provide examples with code when I specifically ask for them.
if there are more than a few changes to be made in more than a few files, please break the task into multiple steps, and check in with me at each step before proceeding to the next one