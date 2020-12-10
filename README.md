# DRUMLR

Drug ranking using machine learning (DRUML) systematically predicts the efficacy of anti-cancer drugs
Henry Gerdes 1, Pedro Casado 1, Arran Dokal 1, Maruan Hijazi 1, #, Nosheen Akhtar1, 3, Ruth Osuntola 4, Vinothini Rajeeve 4, Jude Fitzgibbon 5, Jon Travers 6, David Britton 1,2, Shirin Khorsandi 7 & Pedro R. Cutillas 1,4,8*

1 Cell Signalling & Proteomics Group, Centre for Genomics & Computational Biology, Barts Cancer Institute, Queen Mary University of London, Charterhouse Square, London, EC1M 6BQ, United Kingdom 2 Current address: Kinomica Ltd, Alderley Park, Alderley Edge, Macclesfield SK10 4TG, United Kingdom 3 Department of Biological Sciences, National University of Medical Sciences, Rawalpindi, Pakistan 4 Mass spectrometry Laboratory, Barts Cancer Institute, Queen Mary University of London, Charterhouse Square, London, EC1M 6BQ, United Kingdom 5 Personalised Medicine Group, Centre for Genomics & Computational Biology, Barts Cancer Institute, Queen Mary University of London, Charterhouse Square, London, EC1M 6BQ, United Kingdom 6 Astra Zeneca Ltd, 1 Francis Crick Avenue, Cambridge Biomedical Campus, Cambridge, CB2 0AA, United Kingdom 7 Kings College London, Denmark Hill, Brixton, London SE5 9RS, United Kingdom 8 The Alan Turing Institute, The British Library, 2QR, 96 Euston Rd, London NW1 2DB, United Kingdom

## Copyright
The copyright holder for these data is the author. This resource is made available under a Creative Commons Attribution-NonCommercial-NoDerivatives CC-BY-NC-ND 4.0 International license.

## How to install
### for private use
download the DRUMLR_0.1.0.tar.gz file and run the following code, replacing path with your local path to the downloaded DRUMLR_0.1.0.tar.gz (usually "C:/Users/Username/Downloads/DRUMLRv01_0.1.0.tar.gz")

install.packages(path_to_, 
                 repos = NULL, 
                 type ="source")

### when the package is public
To install the package install the devtools package and run:
devtools::install_github(CutillasLab/DRUMLR)

## Code dependencies 
dplyr (1.0.2),
foreach (1.5.1),
doParallel (1.0.16),
limma (3.42.2), 
caret (6.0-86),
h2o (3.32.0.1) *,
Cubist (0.2.3),
pls(2.7-3),
glmnet (4.0-2),
kernlab (0.9-29),
randomForest (4.6-14)

* this package requires up to date java and H2O packages to operate, however, DRUMLR is written with backwards compatibility for H2O. 
