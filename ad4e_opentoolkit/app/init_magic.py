import os, sys
import shutil
from ad4e_opentoolkit.app.global_var_lib import _repo_dir

destination = None


def init_magic():
    try:
        if len(sys.argv) > 2:
            print("error: only one Argument(Ipython Profile) Accepted")
        if len(sys.argv) > 1 and sys.argv[1].upper() == "HELP":
            print('To copy the magic file to default ipython profile startup directory type just "init_magic "')
            print(
                'To copy the magic file to a specific  ipython profile startup directory type just "init_magic <myprofile>"'
            )
            print('To copy the magic file to current directory type "init_magic ."')
        else:
            inp = None
            if len(sys.argv) == 1:
                destination = os.path.expanduser("~/.ipython/profile_default/startup/openad.py")

            else:
                inp = sys.argv[1]
                destination = os.path.expanduser(f"~/.ipython/profile_{inp}/startup/openad.py")

            if inp == ".":
                destination = os.getcwd() + "/openad.py"

            if os.path.exists(destination):
                os.remove(destination)
            shutil.copyfile(_repo_dir + "/magic/openad.py", destination)
            print("Successully copied openad.py to " + destination)
    except Exception as e:
        print("Error: Unable to copy magic file, destination may not exist")
        print(destination)
        print(e)
