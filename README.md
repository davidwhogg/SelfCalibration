# SelfCalibration

A project on survey design for precise self-calibration.

### Authors:

* Rory Holmes (MPIA)
* David W. Hogg (NYU)

### License:

**Copyright 2012 the authors.  All rights reserved.**

If you want to license this code for use or re-use, get in touch.

### Contributors:

* Hans-Walter Rix (MPIA)

### Dependencies:

* numpy
* matplotlib
* scipy.optimize

### Running the code:

You can make everything happen by running `python something.py` at the
command line in the top-level directory.  Not sure what `something.py`
is.  Rory?

### Known issues:

* The above `something.py` is not right, obviously!
* The true flat-field is too peaked in the center.
* The `dev/` subdirectory is really a branch of the code; needs
  merging.
* Code-generated `.txt` files are checked in.

### Migration from svn:

Hogg migrated this from `svn` with something like

    cd
    git svn clone svn+ssh://astrometry.net/svn/trunk/projects/euclid/ --no-minimize-url --authors-file ~/authors
    cd euclid
    git remote add origin git@github.com:davidwhogg/SelfCalibration.git
    git push origin master
