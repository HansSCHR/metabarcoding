#!/bin/bash
#SCHRIEKE Hans 
#Dada2 on R 


module load system/conda/5.1.0
echo "conda loaded"


read -p "Do you want to install a dada2 environment with conda ? [y/n] " -n 1
if [[ $REPLY =~ ^[Yy]$ ]]
then  
	# dada2 environment installation + run the script
	conda create -n dada2 zlib=1.2.8
	source activate dada2
	conda install -c r r-base=3.5
	conda install gxx_linux-64
	conda install -c bioconda -c conda-forge bioconductor-shortread
	conda install libiconv
	Rscript dada2.R

else [[ $REPLY =~ ^[Nn]$ ]]
	# Just run the script
	source activate dada2
	Rscript dada2.R 
fi
