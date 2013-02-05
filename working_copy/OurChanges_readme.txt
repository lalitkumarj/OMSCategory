Additional changes:
sage -sh
cd $SAGE_ROOT/devel/sage/sage/categories/	#make sure you are in the branch you want
ln -s /path/to/OMS/categories/modsym_space_category.py .
ln -s /path/to/OMS/categories/modsym_coefficient_module_category.py .

In sage.modular.element.py: comment out 2 or 3 line with ps_modsym_...

In sage.schemes.elliptic_curves.ell_rational_field.py: comment out 2 instances of ps_modsym_from_elliptic_curve

In sage.modular.all.py: comment out from btquotients.all import *