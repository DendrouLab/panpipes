import subprocess

def is_package_installed(package_name):
    try:
        subprocess.check_output(["dpkg-query", "-W", package_name])
        return True
    except subprocess.CalledProcessError:
        return False

# Usage:
# replace "time" with name of the package

if is_package_installed("time"):
    print("Package is installed.")
else:
    print("Package is not installed.")
