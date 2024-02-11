# MELT
Mobile Element Locator Tool (MELT) as a Snakemake workflow. 

To generate a mobile element insertion report (WGS only), MELT installation is required and some paths must be added to config_hpf.yaml:
- Download MELT from https://melt.igs.umaryland.edu/downloads.php.
- Unpack the .tar.gz file: tar zxf MELTvX.X.tar.gz  This should create a MELTvX.X directory in your current directory.
- In config_hpf.yaml, add the path to the MELTvX.X directory to config[“tools”][”melt”] .
- Generate a transposon reference text file containing a list of full paths to the mobile element references. For example: ls <full_path_to>/MELTv2.2.2/me_refs/1KGP_Hg19/*_MELT.zip > transposon_file_list.txt 
