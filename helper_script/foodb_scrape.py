# FooDB Scraping

# https://foodb.ca/public/system/downloads/foodb_2020_4_7_csv.tar.gz This is the web file to scrape

import tarfile
import urllib.request

FooDB_file=urllib.request.urlopen("https://foodb.ca/public/system/downloads/foodb_2020_4_7_csv.tar.gz")
print("downloaded")

tarfile.open(name=None, mode='r:gz', fileobj=FooDB_file, bufsize=10240)

