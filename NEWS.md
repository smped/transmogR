# transmogR 0.1.0

* Initial github commit

# transmogR 0.1.1

* Added mogrify functions

# transmogR 0.1.2

* Added support for working directly from VCF files
* Added SNV, Insertion & Deletion calling (`calvInDel`)

# transmogR 0.1.3

* Added UpSet plot for showing variants by gene/transcript
* Added overlapsByVar

# transmogR 0.1.4

* Changed function names
* Added sjFromExons

# transmogR 0.1.5

* Major Bugfix in mogrifyTranscriptome when adding var_labels
* Passed ... to `mclapply()` where relevant
* Moved extraChIPs to Suggests & stopped chopping the output from `sjFromExons()`
* Allowed `sjFromExons()` to return a `GInteractions` object

# transmogR 0.1.6

* Added vignette

# transmogR 0.1.7

* Renamed mogrify*() to `transmogrify()` and `genomogrify()`
* Added masking to `genomogrify()`
