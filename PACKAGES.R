# Run the following lines for packages installation
# This could take a few minutes
# Please install the packages one by one

# 1- EBImage package 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EBImage")
# When the a/s/n option appears, select a

# 2- R.matlab package
install.packages("R.matlab")

# 3- mnormt package
install.packages("mnormt")    

# 4- mixtools package
install.packages("mixtools")
