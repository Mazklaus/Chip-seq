# Readme

This script was made to analyse result of expression microarray experiment from 3' affymetrix chip.
It has been design to work on cluster since it is ram costly.

## Summary

- [Environment installation](https://github.com/Mazklaus/Chip-seq/blob/dev/microArray/readme.md#installing-the-environment)
- [Script use](https://github.com/Mazklaus/Chip-seq/blob/dev/microArray/readme.md#using-the-script)

## Installing the environment

First you will need to install the environment provide by the yml file.
This type of file is used to provide a virtual environment with all the dependency needed to run a script.
Those environment are managed with [miniconda](https://conda.io/miniconda.html) (or anaconda wich is the full package).
It maybe already install in the cluster where you are working but you may have to install it yourself.
If it is the case just dowload the correct file from the link above and follow explication.

Once it has been install you just have to run this command line to create the new environment :

```
conda env -n myEnvName -f pathTo/??.yml
```

Once the environment is created you can use :

```
conda activate myEnvName
```

to place yourself into the environment.

You can quit the environment by using :

```
conda deactivate
```

Just be aware that even if your environment is different all modification on file are persistant, conda just manage packages used. For instance, a file deleted during the use of the virtual environment will be permentaly remove from your pc or the cluster.
For a more complete explanation of how conda work checkout this [link](https://conda.io/docs/commands.html#conda-general-commands).

## Using the script
