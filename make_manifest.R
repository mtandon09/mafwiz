#### writeManifest() displays these warnings for Bioconductor pages
# 
# Warning: 
# * May be unable to deploy package dependency 'AnnotationDbi'; could not determine a repository URL for the
# source 'Bioconductor'.
#
#### This also causes deployment on the connect server to fail with this error
#
# Performing manifest.json to packrat transformation.
# Cannot create packrat data from manifest.json: Cannot transform manifest into packrat.lock: package MutationalPatterns does not define Repository
#
#### This solution adds the bioconductor links to the 'repos' global option
#### Note this does not affect package installation normally because CRAN also checks Bioconductor behind the scenes
#### But for full reproducibility and use with RStudio Connect, we need to provide them both explicitly
#### More here: https://github.com/rstudio/rstudio/issues/6130

# The default result only includes CRAN
# r <- getOption("repos")

# Biocmanager apparently reports both CRAN and bioc urls
bioc <- BiocManager::repositories()

# Set new repos
options(repos = bioc)

# Write manifest (should produce NO warnings)
rsconnect::writeManifest()
