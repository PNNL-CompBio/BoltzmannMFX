This directory contains rich text files with the documentation for BoltzmannMFX.
The Github repository is configured to automatically rebuild this documentation
every time `main` is commited to and push the result to the web at:

     https://pnnl-compbio.github.io/BoltzmannMFX/

To modify these pages, all you need to do is commit something to main with your
changes.

To build this documentation on your local machine, you will need python and a few required packages.
To install them, do

   pip install -r requirements.txt

Then you can type `make html` in this directory to build the docs.
The resulting html files will be stored in `doc/build/html`.