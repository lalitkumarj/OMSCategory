Assumptions: 
    *you have a fresh copy of sage (5.6)
    *the "working_copy" folder of our git repository is located at /path/to/working_copy/
    *(recommended) you've cloned your main branch of sage and this clone is your active branch

# 1. Make sure assumptions are true, in paritular you should have a copy of our git repository on your computer

# 2. Apply patches to the Sage library:
     sage -sh
     cd $SAGE_ROOT/devel/sage
     hg qimport -P /path/to/working_copy/trac_13949.patch
     hg qimport -P /path/to/working_copy/changes_to_sagelib.patch
     hg qimport -P /path/to/working_copy/changes_to_sagelib_fam.patch

# 3. Make symlinks:
     sage -sh
     cd $SAGE_ROOT/devel/sage/sage/modular/
     ln -s /path/to/working_copy/modular/btquotients .
     ln -s /path/to/working_copy/modular/pollack_stevens .
     cd $SAGE_ROOT/devel/sage/sage/modular/overconvergent/
     ln -s /path/to/working_copy/modular/overconvergent/pollack .

# 4. Rebuild:
     sage -b