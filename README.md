# Privacy Preserving Genomic Beacon Data 

This is the original code (with some modification from me) of the [Bustamante Attack on Genomic Beacon Service](http://www.cell.com/ajhg/abstract/S0002-9297(15)00374-2). The code was run using qsub and scripted in R. 

The two privacy preserving solutions are also provided in the scripts folder. The analysis folder consists the make file which is not self-explanatory! Please email me or fire up an issue which I should address (no promises dear).

However I am not adding all the original codes (only the simulated data portions) as well as the data can be found at [Stanford Repo](https://purl.stanford.edu/ps677nj9353). You can also generate these large data with the makefiles, don't worry (thanks to the author). 

For the PGP individual, you can download data from https://my.pgp-hms.org/profile/hu48C4EB - the link to genomic data is http://evidence.personalgenomes.org/genome_download.php?download_genome_id=4f90693e77e7c9b39c9e6dddd43bc95934c3f654&download_nickname=CGI+sample+GS01239-DNA_D01+from+PGP+sample+ . The filename is different from the one in the makefile (but the content matches), so you'll have to rename the unzipped file to CGI_sample_GS01239-DNA_D01_from_PGP_sample_.out  for things to work in the makefile. However I did not use these. 

## Really thankful to Suyash Satyendra Shringarpure for sharing this.

![alt text][logo]

[logo]: http://www.cs.umanitoba.ca/~bruce/umanlogotight.jpg "University of Manitoba"
	
