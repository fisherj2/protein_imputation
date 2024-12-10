
#the purpose of this script is to work around dependency problems with scipenn. By default the conda install from yml is insufficient.
#it will check package versions and update to the version required if needed.

check=$(python -c "
#check for packages and install if needed
import sys
import subprocess
import pkg_resources

required = {'scipenn'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed

if missing:
    print('True')
else:
    print('False')")
    

if [ $check = "True" ]; then
  #this nonsense is needed because scipenn supposedly requires numba==0.50.0 but throws errors unless numba is upgraded
  pip install llvmlite==0.36 --ignore-installed
  pip install sciPENN --ignore-installed
  
  #fix permissions
  chmod 777 -R $CONDA_PREFIX
  conda install numba==0.53.0 
fi




