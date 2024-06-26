# Installation instructions to use the R and python libraries
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

If this happens, if you already have a personal token, this must be added to the file '.Renviron' in your home.

touch ~/.Renviron and edit to add GITHUB_PAT=thetokenabove

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


## Steps followed to install HARP in atos

In order to use the harp libraries it is recommended
to use a local environment created with the `renv` library.
The local environment can be created from scratch using the instructions below.
Remember to replace `{path_to_oper-harp_repository}` by the local
path where this repository was cloned.

**Create renv from scratch**
```bash
module load R/4.2.2 # all tests done with this version of R
module load ecmwf-toolbox #this one is necessary to install Rgrib2 dependencies
module load proj
```

```r
R
library(renv)
renv::init()
Exit R (ctrl-D or exit)

```
Once the renv environment is initiated, enter R again and install
the following libraries

```r
R
renv::install("remotes")
library(remotes)
install_github("harphub/harp")
library(harp)
install.packages("reticulate")
install_github("harphub/Rgrib2")

```

In order to install the `hdf5r` and `ncdf4` libraries needed to read hdf5 and netcdf
files, follow these steps

#### for netcdf
```bash
module load netcdf
R
install.packages("ncdf4")
```

#### for hdf5r
When working in atos, follow the following two extra steps.
If working any other machine, this step is not needed.

For the `hdf5r` library, R is does not add the correct link flags for setting rpath
when building the library. In order to avoid this issue, first create a file
` ~/.R/Makevars` with this line
```
PKG_LIBS = $(HDF5_LIB)

```
Then install `hdf5r` in R:
```bash
module load hdf5
R
install.packages("hdf5r")
```
When this is done, remove the `Makevars` file, as this might interfere with the 
installation of other packages.

After starting R in this location you should see a message like this, pointing
to a local R installation and not the standard `$HOME/x86_64-pc-linux-gnu-library`:

```
- Project '/etc/ecmwf/nfs/dh1_perm_b/miag/ACCORD_VS/R/harp_local_installation' loaded. [renv 1.0.3]
[Previously saved workspace restored]

```

To update `renv.lock` after installation use:
```
renv::settings$snapshot.type("all")
renv::snapshot()

```

### Installing local python environment
Some of the tests included here use the local python3 module and a local python environment
#### TODO Ask Fabrizio: provided venv file for recreating python environment that includes sat-data processing libraries


### Installing the modified libraries `harpSpatial` and `harpIO`

In order to install the modified libraries `harpSpatial` and `harpIO` follow these instructions
after installing the original harp libraries.

```bash
R
remove.packages("harpSpatial")  
remove.packages("harpIO")
install.packages("{path_to_oper-harp_repository}/ACCORD_VS_202406/harpSpatial_branch_DMI_GeoSphere.tar.gz", repos=NULL, type="source")
install.packages("{path_to_oper-harp_repository}/ACCORD_VS_202406/harpIO_branch_DMI_GeoSphere.tar.gz", repos=NULL, type="source")
```
