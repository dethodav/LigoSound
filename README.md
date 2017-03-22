# LigoSound

Ligo Sound is a package to easily produce audio files from Ligo Data. Included are executables to set up ASD directories, generate single audio files, and batch groups of times into audio files. THis library is reliant upon gwpy and gwdetchar, located at 

* [gwpy](https://github.com/gwpy/gwpy)
* [gwdetchar](https://github.com/ligovirgo/gwdetchar)

The executables to create the files themselves (wsf and wsf-batch) are heavily based on similar files in gwdetchar (wdq and wdq-batch), written by Duncan Macleod.

If you are a member of the LIGO Scientific Collaboration, a virtualenv is available for you to use on the LIGO Data Grid, contianing all the dependencies needed for LigoSounds.

If you are using a bash shell, you can access this environment via 

`source /home/detchar/opt/gwpysoft-2.7/bin/activate`
