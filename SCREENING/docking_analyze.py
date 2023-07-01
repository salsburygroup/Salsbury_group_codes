import os
import sys
import shutil
import glob

def run_scripts(screen_type, project_name):
    # Store the original directory
    original_dir = os.getcwd()

    # Define the working directory based on the screen_type and project_name
    work_dir = f'/deac/phy/salsburyGrp/autodock/{project_name}/screen_{screen_type}'

    # Change to the working directory
    os.chdir(work_dir)

    # Run the Python scripts in the working directory
    os.system('python /deac/phy/salsburyGrp/group_python/SCREENING/find_tophits.py')
    os.system('python /deac/phy/salsburyGrp/group_python/SCREENING/sort_diversity_kmeans.py 1000 10')
    os.system('python /deac/phy/salsburyGrp/group_python/SCREENING/sort_diversity_kmeans.py 10 10')
    os.system(f'python /deac/phy/salsburyGrp/group_python/SCREENING/make_report.py most_diverse_10_from10_kmeans top10_{project_name}_{screen_type} /deac/phy/salsburyGrp/autodock/ZINC_libraries_06232023/ZINC_clean_instock_leadlike.smi')
    os.system(f'python /deac/phy/salsburyGrp/group_python/SCREENING/make_report.py most_diverse_10_from1000_kmeans top10_diversity_{project_name}_{screen_type} /deac/phy/salsburyGrp/autodock/ZINC_libraries_06232023/ZINC_clean_instock_leadlike.smi')

    # Copy files starting with 'top10_' back to the original directory
    for file in glob.glob("top10_*"):
        shutil.copy(file, original_dir)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python run_scripts.py <screen_type> <project_name>")
        sys.exit(1)

    screen_type = sys.argv[1]
    project_name = sys.argv[2]
    run_scripts(screen_type, project_name)

