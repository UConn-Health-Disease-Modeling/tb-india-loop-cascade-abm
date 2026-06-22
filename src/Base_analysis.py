# Example CSVAnalyzer for any experiment
# In this example, we will demonstrate how to use a CSVAnalyzer to analyze csv files for experiments

# First, import some necessary system and idmtools packages.
from logging import getLogger

from idmtools.analysis.analyze_manager import AnalyzeManager
from idmtools.core import ItemType
from idmtools.core.platform_factory import Platform

from analyzer_dev.CalibSites_SA import New_SouthAfricaCalibSite
#from CalibSites import SouthAfricaCalibSite

if __name__ == '__main__':
    logger = getLogger()
    with Platform('CALCULON') as platform:


        analyzers = New_SouthAfricaCalibSite().get_analyzers(input_vars = ['Base_Infectivity_Constant']) #is a list with the getter
        analyzers = New_SouthAfricaCalibSite().get_analyzers()  # all variables
        analyzers = New_SouthAfricaCalibSite().get_analyzers(lik_type=None) #turn off likelihood to just get csvs for plotting later etc

        # Set the experiment id you want to analyze
        experiment_test = '64db2357-8984-ec11-a9f3-9440c9be2c51'
        #set a dir to save in locally
        dir = 'test'

        # Specify the id Type, in this case an Experiment on COMPS
        manager = AnalyzeManager(partial_analyze_ok=True, ids=[(experiment_test, ItemType.EXPERIMENT)],
                                 working_dir=dir, force_manager_working_directory=False,
                                 analyzers=analyzers)
        manager.analyze()