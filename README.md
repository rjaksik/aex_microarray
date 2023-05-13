# aex_microarray
Functions used to automatically download and process data from the ArrayExpress database

This is a set of applications designed to test various hypothesis on extremely large number of samples from various microarray experiments. Apart from methods that automatically download microarray data from ArrayExpress database for a selected platform or tools which gather some basic sample data it is designed to harness the possibilities of high performance computers in order to conduct a very efficient analysis using some of the pre-prepared scripts or any other script in Python or R provided by the user. This package includes the following applications:

* **get_aex_data.py** - downloads the entire dataset from the ArrayExpress database for a specific platform
* **get_cel_dates.py** - ectracts sample processing dates (yy/mm/dd) out of each cell file in specified directory
* **get_cel_tm.py** - calculates Tm plot values for each sample in directory based on specified pgc_file
* **get_gc_ctrl.py** - calculates difference in probeset signal changes between two samples and their average slope difference + mean intensity FC, the statistics are averaged between probesest of similar GC
* **get_set_std.py** - calculates inter-proeset variance for each sample based on specified pgc_file

The example usage of each application can be viewed by running the application without any parameters.
