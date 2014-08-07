# in case of using git repository
# git clone git@github.com:wikiselev/$1 $1
# cd $1
# git checkout -b gh-pages origin/gh-pages

# in case of just creating a folder
rm -r $2
mkdir $2
rm -r ~/public_html/misc/$2
mkdir ~/public_html/misc/$2
/homes/kiselev/R-2.15.2/bin/Rscript rmd_publish.R $1 $2