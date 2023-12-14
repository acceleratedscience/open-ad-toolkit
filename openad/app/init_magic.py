import os, sys
import shutil
from openad.app.global_var_lib import _repo_dir


def init_magic():
    """Installs the magic command library into the ipykernel default director or where specified"""
    destination = None
    try:
        if len(sys.argv) > 2:
            print("error: only one Argument(Ipython Profile) Accepted")
        if len(sys.argv) > 1 and sys.argv[1].upper() == "HELP":
            print('To copy the magic file to default ipython profile startup directory type just "init_magic "')
            print(
                'To copy the magic file to a specific ipython profile startup directory\
                      type just "init_magic <myprofile>"'
            )
            print('To copy the magic file to current directory type "init_magic ."')
        else:
            inp = None
            default_path = os.path.expanduser("~/.ipython/profile_default/startup")
            if not os.path.exists(default_path):
                os.mkdir(default_path)
            if len(sys.argv) == 1:
                destination = os.path.expanduser("~/.ipython/profile_default/startup/openad_magic.py")
                destination_old = os.path.expanduser(f"~/.ipython/profile_default/startup/openad.py")
            else:
                inp = sys.argv[1]
                destination = os.path.expanduser(f"~/.ipython/profile_{inp}/startup/openad_magic.py")
                destination_old = os.path.expanduser(f"~/.ipython/profile_{inp}/startup/openad.py")

            if inp == ".":
                destination = os.getcwd() + "/openad_magic.py"

            if os.path.exists(destination):
                os.remove(destination)
            if os.path.exists(destination_old):
                os.remove(destination_old)
            shutil.copyfile(_repo_dir + "/magic/openad_magic.py", destination)
            print("Successully copied openad_magic.py to " + destination)
    except Exception as e:
        print("Error: Unable to copy magic file, destination may not exist")
        print(destination)
        print(e)
