# Installation instructions to use the R and python libraries to be used during the ACCORD 2025 visit at GeoSphere
Instructions below are focused on using the scripts
in the ECMWF atos machine.
Since the installation makes heavy use of github, it is
a good idea to create a personal token in github. 

## Creating a personal token in github

Occasionally the harp installation can fail with a message like:
``` r
Downloading GitHub repo andrew-MET/harp@master
Error: HTTP error 403.
 API rate limit exceeded for 130.226.71.190. (But here's the good news: Authenticated requests get a hi
gher rate limit. Check out the documentation for more details.)
```
or:
``` r
Error: Failed to install 'unknown package' from GitHub:
```

If this happens, if you already have a personal token, this must be added to the file '.Renviron' in your home:

``` r
touch ~/.Renviron and edit to add GITHUB_PAT=thetokenabove
```

If you do not have a token, follow the instructions [here](https://happygitwithr.com/https-pat.html#get-a-pat) 

Note that to add the token to .Renviron, you can also follow these instructions:

``` r
gitcreds::gitcreds_set()
<Enter token here>
```
Then use:
``` r
usethis::edit_r_environ()
```
This will open an editor to edit the file `.Renviron`, where you can write GITHUB_PAT=thetokenabove.

The system will ask you to restart R to take effect. Then try to install harp again.

Keep the personal token in your `.Renviron` for later use. The system will detect it an use it every time you use `install_github`


## Steps followed for Installation

In order to use the harp libraries it is recommended
to use a local environment created with the `renv` library.
The local environment can be created from scratch using the instructions below.


**On ATOS**
(tested)

Module was used (load before running R)  
```
module load R
module load ecmwf-toolbox  # this one is necessary to install Rgrib2 dependencies
module load proj
module load hdf5
module load netcdf4

This gives the following default moduels (July 2025)
R/4.4.3
hdf5/1.14.6
ecmwf-toolbox/2025.04.0.0
proj/9.5.1
netcdf4/4.9.3
```

**Using a conda environment** 
(...to be tested and updated by Polly...)

Use `conda_environment_R_4.3.3.txt` file to get conda environment.

Or do it step by step:

```
conda create --name <env_name>
conda activate <env_name>
conda install conda-forge::r-base proj r-hdf5r r-mass metview r-matrix r-ragg

Start R session

install.packages("renv")  
Exit R (ctrl-D or exit)
```

## install HARP

Choose the path where you want to have your local installation:

cd <harp_local_installation>

**Create renv from scratch**
```
Start R session

library(renv, bare = TRUE)
renv::init()
Exit R (ctrl-D or exit)
```

Since the renv environment has just been initiated, after starting R in your <harp_local_installation> you should now see a message like this, pointing
to a local R installation and not the standard `$HOME/x86_64-pc-linux-gnu-library`:

```
- Project '/etc/ecmwf/nfs/dh1_perm_b/miag/ACCORD_VS/testing/installHarp' loaded. [renv 1.0.7]
```

Enter R and install the following libraries
Note that we will be using some of the development versions here.
`harpCore` development version is important to
be able to read the INCA data.


```
Start R session

install.packages("remotes")
library(remotes)
install_github("harphub/harp")
install_github("harphub/harpIO", "develop")
install_github("harphub/harpSpatial", "develop")
install_github("harpIO/harpCore","develop")

###install_github("pollyaschm/harpSpatial", "ACCORD_VS_202406")
```
**Note:**

When installing install_github("harphub/harpSpatial", "ACCORD_VS_202507"), do not update other packages (option 3: None) 

```
Changes available in harpIO/develop

harpIO:
changes to read_hdf5()
such that DMI data does not come out inverted

harpSpatial (master)
Now includes changes introduced by Polly in former VS
Installing develop version here though...

```

**Proceed with the installation of the remaining packages.**

```
install_github("harphub/Rgrib2")
```

Install hdf5 (needed when you need to access hdf files, e.g. DMI radar products)

When working on **ATOS**, follow the instructions below.


```
Exit R (ctrl-D or exit)

mkdir ~/.R/Makevars
in Makevars add: PKG_LIBS = $(HDF5_LIB)

Start R session

install.packages("hdf5r")

Exit R (ctrl-D or exit)

After successful installation remove the `Makevars` file, as this might interfere with the installation of other packages.
remove ~/.R/Makevars 
```
If working in any other machine, this step should not be necessary.
Simply do:
```
install.packages("hdf5r")
```

**Continue** installing packages

```
Start R session

install.packages("reticulate")
install.packages("here")
install.packages("tidyverse")
```

When using **conda environment**, there might be a problem with a wrong libtiff.so.* when trying to install tidyverse.

(libtiff.so.5: cannot open shared object file - while libtiff.so.6 is available). In this case "ragg" cannot be installed when trying to install "tidyverse". This can be tricked like this:

```
cd /path_to_conda_env/lib/
ln -s libtiff.so.5 libtiff.so.6
```
Then install.packages("tidyverse") should work.

**Continue** installing packages
```
install.packages("dplyr")
install.packages("ggpubr")
install.packages("RColorBrewer")
install.packages("ggplotify")
install.packages("patchwork")
```

In case you need to handle netcdf files

```
module load netcdf4
Start R session
install.packages("ncdf4")
```


To update `renv.lock` after installation use:

```
renv::settings$snapshot.type("all")
renv::snapshot()
```

## Additional dependencies 

### When using reticulate to interface with python
(More info in reading_functions.py)

We use:

```
module load python3 (python3/3.11.8-01)
module load ecmwf-toolbox (ecmwf-toolbox/2024.04.0.0, to have available the python metview interface)  
```

We need a virtual enviroment where satpy is installed (satpy is used to ahndle MSG data)

Install your own python env:

```
python -m venv satpy
pip install satpy
```

or you might source /perm/miag/venvs/satpy/bin/python3 ('rx' permissions for accord group)

The path to your python is explicitly mentioned in the definition files for reticulate examples.
Change the path to python accordingly. 

### Note on "non standard" grib parameters

```
When reading grib files which contain unknown parameters, as the precipitation parameters of DEODE:

* Modify the ECCODE grib definitions accordingly (see definitions/).

* Add ECCODE DEFINITION_PATH to your .bashrc:
  export ECCODES_DEFINITION_PATH=<path_to_definitions>/definitions:ECCODES_DEFINITION_PATH
```



