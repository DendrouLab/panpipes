import subprocess
import panpipes.R_scripts

def install_r_dependencies():
    r_path=panpipes.R_scripts.__path__
    subprocess.check_process(["Rscript", r_path + "install_R_libs.R"])

def main():
    install_r_dependencies()

if __name__ == "__main__":
    main()