# Hello Variant Hunter!

This is VariantHunter version 0.5.0

VariantHunter is written by Alex Thorn.

The scripts folder contain code from other people, see the README in the scripts folder.

If you want to run VariantHunter on computerome then you can find  the correct modules to load in modules_load.txt.

Fill out the config file and test that everything works by typing "snakemake -np" in the commandline. 

To run the pipeline for real type "snakemake --cores ?", where ? is how many cores that are available.

More commanline instructions can be found here: https://snakemake.readthedocs.io/en/stable/executing/cli.html

Notes on running the pipeline:
-Mapping with kma is the ratelimiting step.
-The output folders are ordered 1-5 according to when they are created in the pipeline.
-To rerun part of the pipeline without running the mapping again simply delete the relavant output folders.
-For example to run the consensus sequence filtering again with new parameters simply delete folders 3, 4, and 5 and update the config file.
-Then run snakemake again. This will take shorter time than the entire run since the mapping has already been done.

Log files can be found in the hidden snakemake folder. Try "cd .snakemake/log"

Outout folder 5_results contain the actual results that you may be interested in.
