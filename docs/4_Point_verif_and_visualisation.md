# Point verification

Once the configuration file is set and the sqlite tables are created, the point verification can be carried out. 

## set_params.R

This parameter list file is used to specify the parameters considered by the `point_verif.R` script and their associated options, in particular:
- **scale_fcst**: Forecast scaling (e.g. Kelvin to degress).
- **scale_obs**: Observation scaling.
- **thresholds**: Thresholds used when computing threshold skill scores.
- **obsmin/max_val**: Max/min observation values allowed.
- **fctmax_val**: Max forecast values allowed (experimental).
- **error_sd**: Number of standard deviations used in `harpPoint::check_obs_against_fcst`. 

Typically this file does not need to be changed. By default `point_verif.R` reads parameter options from this file, but a custom parameter file can also be used by passing the `-params_file` option to `point_verif.R`. 

## point_verif.R

### Inputs and usage 

This script takes the following command line inputs (required arguments in **bold**, optional arguments are in *italics*):

- **-config_file**: The config file in the `config_files` directory (no default).
- **-start_date**: The first forecast cycle to process (in YYYYMMDDHH format, no default).
- **-end_date**: The last forecast cycle to process (in YYYYMMDDHH format, no default).
- *-params_file*: The parameter list file containing parameter scalings, thresholds, etc. (default="verification/set_params.R").
- *-params_list*: Which parameters for verify (default="ALL"). This should be a comma separated string of parameters, for example "T2m,S10m,T,S". These parameters should exist in the parameter list file, otherwise they will be skipped. If `params_list` is not specified, all parameters in the parameter list file are considered in the verification (this is not recommended in general).
- *-mod_def_rds*: A logical flag to prepend the project name to harp's default rds filenames (default=FALSE). Not generally required.
- *-add_proj_png*: A logical flag to prepend the project name to the default png filenames (default=FALSE). Not generally required. 
- *-rolling_verif*: A logcial flag to indicate "rolling" verification (default=FALSE). If TRUE, rolling verification will produce a reduced set of png files and will not produce rds files or scorecards. Generally rolling verification is restricted to a "short" (e.g. 7 days) near-real time period. This option is not compatible with `gen_sc_only=TRUE`. 
- *-gen_sc_only*: A logical flag to run scorecard generation (and plotting) only (default=FALSE). This may be useful in cases where point verification results have already been generated. This option is not compatible with `rolling_verif=TRUE`. 
- *-use_fixed_dates*: A logical flag to use the input `start_date` and `end_date` when naming the directories and png files associated with this verification (default=TRUE). If set to FALSE, the data generated will use start and end dates corresponding to the first and last `fcst_dttm`, respectively, used in the verification for a given parameter. Therefore if set to FALSE, data may be stored in different directories for different parameters if the first and last `fcst_dttm` differs (this can happen in particular for precipitation). This option does not change the start/end dates in the rds filenames. 

Alternatively, the script can be sourced directly from within R/RStudio. In this case, the arguments will be read from the `verification/source_options.R` file. This interactive mode is useful for interrogating the data. 

Typical usage:
``` 
./point_verif.R -config_file config_files/config_det_example.yml -start_date YYYYMMDDHH -end_date YYYYMMDDHH -params_file verification/set_params.R -params_list "T2m,Q"
```
This will run the verifcation for using cycles from `start_date` to `end_date` in steps of `verif:by_step` for parameters T2m and Q.

### QC

The following QC checks of the forecast and observation data are carried out by the script:

- Forcast values above `fctmax_val` are removed (if this variable is set).
- Observation values above/below the `obsmin_val`/`obsmax_val` are removed (if these variables are set).
- `harpPoint::check_obs_against_fcst` is run to remove observations which are more than `error_sd` standard devations away from the forecast.
- Station report frequency is computed and stations in the bottom 1% are removed. This acts to discard stations with very few observation reports. 

### Verification groups

By default the script assumes the following verificaiton for surface and upper-air variables (see `harpCore::make_verif_groups()` for more information):
- Surface: Data is grouped by `fcst_cycle` and the SID grouping specified by the `verif:domains` option in the config file (this is stored under the variable `station_group` in the scripts). Verification scores are then computed as a function of leadtime, valid date, and valid hour.
- Upper air: Data is group by the SID grouping specified by the `verif:domains` option and scores computed as a function of leadtime and valid hour. 

### Output

`point_verif.R` will produce standard harp `.rds` files which contain the full suite of verification scores available in harp by default. These files will be stored in:
```
{verif:verif_path}/{verif:project_name}/harpPointVerif.harp.{parameter}.harp.{start_date}-{end_date}.harp.{forecast_model_1}.model.{forecast_model_2).model...{forecast_model_N}.rds
```
Typically the filenames for the harp rds files should not be changed as the harp shiny app assumes a set format. Note that while the filenames do not contain information about the `verif:domains` considered, the domain selection is included in the rds files under the `station_group` variable. 

If `post:create_png: TRUE`, then a suite of standard verification scores are plotted as png files for local visualisation. These local files will also include plots which are not available in harp's shiny app, such as forecast timeseries and station bias/rmse maps. These files should appear in:
```
{post:plot_output}/{verif:project_name}/{start_date}-{end_date}/long_file_names.png 
```
(the filenames are somewhat convoluted and should not be changed as a strict structure is assumed in the local shiny app). See `fn_png_name` in `fn_plot_helpers.R` for more information on the convention used.

If you are also generating a scorecard (`scorecards:create_scrd: TRUE`), then all the underlying scorecard data (for each domain) will be saved to:
```
{verif:verif_path}/{verif:project_name}/harpScData-{start_date}-{end_date}-{scorecards:ref_model}-{scorecards:fcst_model).rds
```
Scorecard plots will also be available in
```
{post:plot_output}/{verif:project_name}/{start_date}-{end_date}/*scard*.png
```
The scorecard data can be passed to `fn_scorecard_signif.R` for visualisation if desired. If `scorecards:plot_signif: TRUE`, these images will be available in the same directory (search for `*sdiff*.png`).

## Visualisation 

The rds files can be visualised by using harp's built-in shiny app:
```
shiny_plot_point_verif(start_dir={verif:verif_path}/{verif:project_name})
```
For example, when loading in a surface parameter, you should see groups for `station_group` and `fcst_cycle` and different options for the "Time" axis.

To browse the png images easily, a simple shiny app is provided with a similar interface to the "monitor" tool. This shiny app is defined in `visualization/visapp/app.R`. Some small edits to `app.R` may need to be taken in order to run the app locally, as detailed below.

1. The variable `img_dir` is set to NULL by default, which just looks for files in the `sample_images` directory provided in visapp. Set `img_dir <- {plot_output}` instead to point to the images available there after verification.
2. Optional: By default the variable `app_dir` is set to `here("visualization/visapp")` and this should correctly point to the `visapp` directory provided that `here()` starts in this repository. If this is not the case for some reason, alternatively you can hardcode `img_dir` to the location of the `visapp` directory on your system. 
3. Optional: Add your project name to `all_experiment_names.R` in the `visapp` directory, following the convention given, in order to change the display name of the project in the app. If your project is not listed in `all_experiment_names.R`, it will still appear in the app but with a default display name.
4. Optional: There is a `smr_ind` variable which changes the date display in the app to cater for dates categrorsied into "Rolling", "Monthly", and "Seasonal" periods. This assumes that data in `plot_output` is stored under directories following the strcuture `Monthly/project_name/start_date-end_date`, and the same for `Rolling` and `Seasonal`. This switch can generally be left as FALSE. 

Note that you can have images for multiple projects in `plot_output` and switch easily between projects in the app. The app will read all directories which exist in `plot_output`, and therefore you need to remove the `project_name` directory from `plot_output` in order to remove the project from the app. The app can then be launched by opening `app.R` in Rstudio and hitting "Run App". From a terminal, cd into the visualization directory, open R, and use:

``` r
library(shiny)
runApp("visapp")

```


