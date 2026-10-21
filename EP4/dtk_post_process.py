#!/usr/bin/python

import os
import sys
import json
import emod_api.config.dtk_post_process_adhocevents as dpp
import emod_api.spatialreports.spatial as sr
import pandas as pd

def application( output_path ):
    dpp.application(output_path)
    os.system( "Rscript Assets/python/preprocess.R" )
    #os.remove( "output/ReportEventRecorder.csv" )