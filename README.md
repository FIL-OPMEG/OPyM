# OPyM
**work in progress**


 Python library to read OPM data into MNE-Python
 
 
### Usage - FIL BIDS format

To use the import function, it is assumed that the the data is organised in a BIDS-esque format, such that the files are named: 

```
data_folder
│
├───[FILE_PREFIX]_meg.bin             [REQUIRED]
├───[FILE_PREFIX]_meg.json            [REQUIRED]
├───[FILE_PREFIX]_channels.tsv        [REQUIRED]
├───[FILE_PREFIX]_positions.tsv       [OPTIONAL]
└───[FILE_PREFIX]_coordsystem.json    [OPTIONAL]

```
the code will try and guess the name of the files based on the name of the `.bin` file supplied. 

To import and use the funciton.

```python

from opym import read_raw_ucl

raw = read_raw_ucl('/path/to/data_folder/[FILE_PREFIX]_meg.bin')

```

### Usage - Cerca format

**warning:** Only sensor level import currently works, no position informaiton available... yet!

To use the import function, it is assumed that the the data is organised in a BIDS-esque format, such that the files are named: 

```
data_folder
│
├───[FILE_PREFIX].cMEG                [REQUIRED]
└───[FILE_PREFIX]_SessionInfo.txt     [REQUIRED]

```
the code will try and guess the name of the files based on the name of the `.bin` file supplied. 

To import and use the funciton.

```python

from opym import read_raw_cerca

raw = read_raw_cerca('/path/to/data_folder/[FILE_PREFIX].cMEG')

```
