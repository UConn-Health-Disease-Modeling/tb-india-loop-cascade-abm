#
# This is a user-modifiable Python file designed to be a set of simple input file and directory settings that you can choose and change.
#

# The script is going to use this to store the downloaded schema file. Create 'download' directory or change to your preferred (existing) location.
schema_file="download/schema.json"
# The script is going to use this to store the downloaded Eradication binary. Create 'download' directory or change to your preferred (existing) location.
eradication_path="download/Eradication"  # for linux
# eradication_path="download/Eradication.exe"  # for windows
# Create 'Assets' directory or change to a path you prefer. idmtools will upload files found here.
assets_input_dir= "Assets"
plugins_folder = "download/reporter_plugins"

my_ep4_assets = None  # or ['dtk_post_process.py']
ep4="EP4"
ep4_path = "python_scripts"   # path to pre/post_process
requirements = "./requirements.txt"
sif="dtk_centos.id"
#males="data/Malawi_male_mortality.csv"
#females="data/Malawi_female_mortality.csv"
long_term_migration_rates_path="ltm.csv"
india_mortality_data="data/India_mortality_1990_to_2017.csv"
small_pop_file="twonodes.csv"


# ep4="EP4"
# R_script_dir="EP4/R"
# requirements = "./requirements.txt"
# sif="sif.id"
# small_pop_file="twonodes.csv"