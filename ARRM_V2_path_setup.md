# ARRM_V2 Path Setup


The main **ARRM_V2** programs are in the base ARRM_V2 folder;  various additional **ARRM_V2** functions are in several subfolders.

To run **ARRM_V2**, make sure you have the following paths in your MATLAB search path:

* *ARRM_V2_base*/util_general
* *ARRM_V2_base*/util_nc
* *ARRM_V2_base*/ncdf
* *ARRM_V2_base*/ARRM_V2_subs
* *ARRM_V2_base*

where *ARRM_V2_base* is the base folder for the **ARRM_V2** code


Running the function `ARRM_V2_setpath()` will attempt add them to the top of the current searchpath if they are not already present.  Note that it will not save the path;  to do so, run `savepath` .
`ARRM_V2_setpath()` is called at the beginning of most of the **ARRM_V2** main routines.

| Folder Name | Contents |
| ----------- | ---------|
|ARRM_V2   | ARRM_V2 main routines  |
|ARRM_V2_subs   | core subfunctions   |
|ncdf   | ncdf class (handle class based on MATLAB finfo struct)  |
|util_nc   | additional routines for working with netcdf files   |
|util_general   | general utility routines (date/time, etc.)  |

You do **not** need the following folders in your MATLAB path:

| Folder Name | Contents |
| ----------- | ---------|
| ARRM_V2_base/dev_tools | varous routines for testing development code |
| ARRM_V2_base/old_code | old code - only available in some branches. **not in master** |
| ARRM_V2_base/utils | useful (maybe) utilities related to ARRM_V2 |

