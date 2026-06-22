"""
TB India Base Simulation Runner
================================
Clean, standalone script to run the TB India EMOD simulation on COMPS.

BEFORE RUNNING: Set up the environment with setup_env.bat
Then activate:  conda activate tb_india

Usage:
    python run_base_simulation.py              # Run base case (5 replicates)
    python run_base_simulation.py --calibrate  # Run calibration mode
    python run_base_simulation.py --dry-run    # Build config/campaign only (no COMPS submission)
"""

import sys
import os
import json
import copy
import pathlib
import argparse
import numpy as np
from functools import partial

# ── Parse arguments ──────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description='TB India EMOD Simulation')
parser.add_argument('--calibrate', action='store_true', help='Run in calibration mode')
parser.add_argument('--dry-run', action='store_true', help='Build config only, do not submit to COMPS')
parser.add_argument('--pop', type=int, default=1000, help='Initial population size (default: 1000)')
parser.add_argument('--replicates', type=int, default=5, help='Number of replicates (default: 5)')
parser.add_argument('--platform', type=str, default='CALCULON', help='COMPS platform (default: CALCULON)')
args = parser.parse_args()

# ── Change to project directory ──────────────────────────────────────────
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(PROJECT_DIR)
sys.path.insert(0, PROJECT_DIR)

# ── Local imports ────────────────────────────────────────────────────────
import params
import set_config
import manifest
from analyzer_dev.CalibSites_India import IndiaCalibSite
from add_tb_drug_type import add_tb_drug_type, add_tb_drug

# ── EMOD imports ─────────────────────────────────────────────────────────
from emodpy import emod_task
from emodpy.emod_task import EMODTask
from emodpy_tbhiv.interventions import purge_campaign_event
import emod_tbhiv.bootstrap as dtk

from idmtools.core.platform_factory import Platform

# ── Verify environment ───────────────────────────────────────────────────
import emodpy
import emod_api
print(f"Python:        {sys.version}")
print(f"emodpy:        {emodpy.__version__}")
print(f"emod-api:      {emod_api.__version__}")
print(f"Project dir:   {PROJECT_DIR}")
print(f"Platform:      {args.platform}")
print(f"Mode:          {'calibration' if args.calibrate else 'base case'}")
print(f"Replicates:    {args.replicates}")
print(f"Population:    {args.pop}")
print()


# ══════════════════════════════════════════════════════════════════════════
# 1. CONFIG BUILDER
# ══════════════════════════════════════════════════════════════════════════

def set_param_fn(config):
    """Callback to set EMOD config parameters."""
    config = set_config.set_config(config)
    config.parameters.TB_Enable_Exogenous = 1
    config.parameters.Enable_Strain_Tracking = 1
    config.parameters.Death_Rate_Dependence = "NONDISEASE_MORTALITY_BY_YEAR_AND_AGE_FOR_EACH_GENDER"
    return config


def set_duration_no_hiv(task, death_fraction, duration):
    """Convert from years to days, set cure/mortality rates."""
    duration *= 365.0
    loc_rate_recover = (1.0 - death_fraction) / duration
    loc_rate_mort = death_fraction / (1.0 - death_fraction) * loc_rate_recover

    task.config.parameters.TB_Active_Cure_Rate = loc_rate_recover
    task.config.parameters.TB_Active_Mortality_Rate = loc_rate_mort

    return {'TB_Active_Cure_Rate': loc_rate_recover,
            'TB_Active_Mortality_Rate': loc_rate_mort}


def set_duration_hiv(task, death_fraction, duration):
    """Set duration and death fraction for HIV+ TB."""
    duration *= 365.0
    loc_rate_recover = (1.0 - death_fraction) / duration
    loc_rate_mort = death_fraction / (1.0 - death_fraction) * loc_rate_recover

    task.config.parameters.TB_Active_Cure_Rate_HIV = loc_rate_recover
    task.config.parameters.TB_Active_Mortality_Rate_HIV = loc_rate_mort

    return {'TB_Active_Cure_Rate_HIV': loc_rate_recover,
            'TB_Active_Mortality_Rate_HIV': loc_rate_mort}


def set_primary_hiv_pro(task, pro):
    """Scale TB_CD4_Primary_Progression by factor."""
    tmp_loc = []
    for i, val in enumerate(task.config.parameters.TB_CD4_Primary_Progression):
        tmp_loc += [pro * val]
    task.config.parameters.TB_CD4_Primary_Progression = tmp_loc
    return {'TB_CD4_Primary_Progression': tmp_loc}


def update_more_config(task):
    """Additional config parameter updates."""
    val_dur = 2.5
    val_frac = 0.7
    set_duration_no_hiv(task=task, death_fraction=val_frac, duration=val_dur)
    task.config.parameters.TB_Presymptomatic_Rate = 1.0 / (0.1 * 365)
    task.config.parameters.TB_Presymptomatic_Cure_Rate = 0
    # Note: HIV co-infection mortality handled by CoInfection_Mortality_Rate_Off/On_ART in set_config.py
    # set_duration_hiv and set_primary_hiv_pro use params not in this EMOD schema version
    task.config.parameters.Base_Individual_Sample_Rate = 0.05
    task.config.parameters.x_Base_Population = 5
    task.config.parameters.Simulation_Duration = (
        params.burn_initial + params.burn_predots + params.To_end_from_DOTS
    )
    task.config.parameters.x_Birth = 0.850


# ══════════════════════════════════════════════════════════════════════════
# 2. CAMPAIGN BUILDER
# ══════════════════════════════════════════════════════════════════════════

def build_camp(probability_per_step_slow=0.0208, probability_per_step_fast=0.05247, personal_hc_budget=5):
    """Build campaign with care-seeking, diagnostics, and treatment interventions."""

    import emod_api.campaign as camp
    camp.schema_path = manifest.schema_file
    import emod_api.interventions.outbreak as ob
    import emod_api.interventions.common as common
    import emodpy_tbhiv.interventions.hsb as hsb
    import emodpy_tbhiv.interventions.art as art
    import emodpy_tbhiv.interventions.diag_treat_neg as dtn
    import emodpy_tbhiv.interventions.tbhiv_treat as treat
    import emodpy_tbhiv.interventions.ramp_dtn as ramp
    import emodpy_tbhiv.interventions.resist_diag as rd
    import emodpy_tbhiv.interventions.hiv_diag as hd
    import emodpy_tbhiv.interventions.hiv_seeding as hs
    import emodpy_tbhiv.interventions.cd4diag as cd4
    import emodpy_tbhiv.interventions.cash as cash

    common.old_adhoc_trigger_style = False

    # Initial TB outbreak seeding
    seed = ob.seed_by_coverage(camp, timestep=40, coverage=0.05)
    camp.add(seed, first=True)
    purge_campaign_event(seed)

    # ── Health Seeking Behavior ──────────────────────────────────────────

    # Low quality care - pre-DOTS era
    camp.add(hsb.HSB(
        camp=camp,
        trigger_treatment_list=['TBActivation', 'TBFailedDrugRegimen',
                                'TBTestDefault', 'TBTestNegative', 'TBMDRTestDefault'],
        output_event='TBTestDOTSLow',
        probability_per_step=probability_per_step_slow,
        start_day=params.burn_initial,
        duration=params.burn_predots,
        property_restrictions_list=['QualityOfCare:Low'],
        fixed_initial_delay=7
    ))

    # Low quality care - DOTS era onward
    camp.add(hsb.HSB(
        camp=camp,
        trigger_treatment_list=['TBActivation', 'TBFailedDrugRegimen',
                                'TBTestDefault', 'TBTestNegative', 'TBMDRTestDefault'],
        output_event='TBTestDOTSLow',
        probability_per_step=probability_per_step_slow,
        start_day=params.burn_predots + params.burn_initial + 1,
        duration=-1,
        property_restrictions_list=['QualityOfCare:Low'],
        fixed_initial_delay=7
    ))

    # Low quality - treatment seeking after positive test
    camp.add(hsb.HSB(
        camp=camp,
        trigger_treatment_list=['TBTestPositive', 'TBDS_Positive'],
        output_event='TB_Tx',
        probability_per_step=probability_per_step_slow,
        start_day=params.burn_initial,
        duration=-1,
        property_restrictions_list=['QualityOfCare:Low'],
        fixed_initial_delay=7
    ))

    # High quality care - pre-DOTS era
    camp.add(hsb.HSB(
        camp=camp,
        trigger_treatment_list=['TBActivation', 'TBFailedDrugRegimen',
                                'TBTestDefault', 'TBTestNegative', 'TBMDRTestDefault'],
        output_event='TBTestDOTSHigh',
        probability_per_step=probability_per_step_slow * 2,
        start_day=params.burn_initial,
        duration=params.burn_predots,
        property_restrictions_list=['QualityOfCare:High'],
        fixed_initial_delay=7
    ))

    # High quality care - DOTS era onward
    camp.add(hsb.HSB(
        camp=camp,
        trigger_treatment_list=['TBActivation', 'TBFailedDrugRegimen',
                                'TBTestDefault', 'TBTestNegative', 'TBMDRTestDefault'],
        output_event='TBTestDOTSHigh',
        probability_per_step=probability_per_step_slow * 2,
        start_day=params.burn_predots + params.burn_initial + 1,
        duration=-1,
        property_restrictions_list=['QualityOfCare:High'],
        fixed_initial_delay=7
    ))

    # High quality - treatment seeking after positive test
    camp.add(hsb.HSB(
        camp=camp,
        trigger_treatment_list=['TBTestPositive', 'TBDS_Positive'],
        output_event='TB_Tx',
        probability_per_step=probability_per_step_slow * 2,
        start_day=params.burn_initial,
        duration=-1,
        property_restrictions_list=['QualityOfCare:High'],
        fixed_initial_delay=7
    ))

    # Stop testing loop when treated or cleared
    stop_testing = hsb.HSB(
        camp=camp, start_day=1,
        trigger_treatment_list=["TbDrugDistributed", "SymptomaticCleared"],
        output_event="Got_Better_Stop_Testing_Loop",
        probability_per_step=1.0
    )
    stop_testing.Event_Coordinator_Config.Intervention_Config.\
        Actual_IndividualIntervention_Config.Enable_Intervention_Replacement = 1
    camp.add(stop_testing)

    # ── Epoch 1: Pre-DOTS diagnostics & treatment ────────────────────────

    camp.add(dtn.DiagnosticTreatNeg(
        camp=camp, trigger_treatment_list=['TBTestDOTSHigh'],
        base_sensitivity_smearpos=params.sens_smear_pos_pre_GH,
        base_sensitivity_smearneg=params.sens_smear_neg_pre_GH,
        treatment_fraction=0.6, start_day=params.burn_initial,
        duration=params.burn_predots,
        property_restrictions_list=['QualityOfCare:High']
    ))
    camp.add(dtn.DiagnosticTreatNeg(
        camp=camp, trigger_treatment_list=['TBTestDOTSLow'],
        base_sensitivity_smearpos=params.sens_smear_pos_pre_GL,
        base_sensitivity_smearneg=params.sens_smear_neg_pre_GL,
        treatment_fraction=0.6, start_day=params.burn_initial,
        duration=params.burn_predots,
        property_restrictions_list=['QualityOfCare:Low']
    ))
    camp.add(treat.TBHIVDrugTreatment(
        camp=camp, trigger_treatment_list=['TBTestPositive'],
        drug_name='PreDOTSHigh', latent_multiplier=0,
        start_day=params.burn_initial, duration=params.burn_predots,
        property_restrictions_list=['QualityOfCare:High']
    ))
    camp.add(treat.TBHIVDrugTreatment(
        camp=camp, trigger_treatment_list=['TBTestPositive'],
        drug_name='PreDOTSLow', latent_multiplier=0,
        start_day=params.burn_initial, duration=params.burn_predots,
        property_restrictions_list=['QualityOfCare:Low']
    ))

    # ── Epoch 2: DOTS era (pre-GeneXpert) diagnostics & treatment ────────

    dots_to_gx = params.genexpert_introduction - params.dots_start

    camp.add(dtn.DiagnosticTreatNeg(
        camp=camp, trigger_treatment_list=['TBTestDOTSHigh'],
        base_sensitivity_smearpos=params.sens_smear_pos_pre_GH,
        base_sensitivity_smearneg=params.sens_smear_neg_pre_GH,
        treatment_fraction=0.7, start_day=params.dots_start,
        duration=dots_to_gx,
        property_restrictions_list=['QualityOfCare:High']
    ))
    camp.add(dtn.DiagnosticTreatNeg(
        camp=camp, trigger_treatment_list=['TBTestDOTSLow'],
        base_sensitivity_smearpos=params.sens_smear_pos_pre_GL,
        base_sensitivity_smearneg=params.sens_smear_neg_pre_GL,
        treatment_fraction=0.7, start_day=params.dots_start,
        duration=dots_to_gx,
        property_restrictions_list=['QualityOfCare:Low']
    ))
    for qoc, drug in [('High', 'DOTSHQ'), ('Low', 'DOTSLQ')]:
        for triggers in [['TBTestPositive', 'TBDS_Positive'], ['TB_Tx']]:
            camp.add(treat.TBHIVDrugTreatment(
                camp=camp, trigger_treatment_list=triggers,
                drug_name=drug, latent_multiplier=0,
                start_day=params.dots_start, duration=dots_to_gx,
                property_restrictions_list=[f'QualityOfCare:{qoc}']
            ))
    camp.add(treat.TBHIVDrugTreatment(
        camp=camp, trigger_treatment_list=['TBMDRTestPositive'],
        drug_name='DOTSMDR', latent_multiplier=0,
        start_day=params.dots_start, duration=dots_to_gx,
        property_restrictions_list=[]
    ))

    # ── Epoch 3: GeneXpert era diagnostics & treatment ───────────────────

    # Ramped GeneXpert diagnostics
    camp.add(ramp.RampDTN(
        camp=camp, trigger_treatment_list=['TBTestDOTSHigh'],
        ramp_time=params.length_gene_xpert_ramp,
        base_sensitivity_smearpos=params.sens_smear_pos_GXH,
        base_sensitivity_smearneg=params.sens_smear_neg_GXH,
        base_sensitivity_smearpos2=params.sens_smear_pos_pre_GH,
        base_sensitivity_smearneg2=params.sens_smear_neg_pre_GH,
        treatment_fraction=0.65,   # was 0.8 — ~35% pre-Tx LTFU in high QoC
        pos_event='ProviderOrdersTBTest', pos_event2='ProviderTestNoR',
        start_day=params.genexpert_introduction, duration=-1,
        property_restrictions_list=['QualityOfCare:High']
    ))
    camp.add(ramp.RampDTN(
        camp=camp, trigger_treatment_list=['TBTestDOTSLow'],
        ramp_time=params.length_gene_xpert_ramp,
        base_sensitivity_smearpos=params.sens_smear_pos_GXL,
        base_sensitivity_smearneg=params.sens_smear_neg_GXL,
        base_sensitivity_smearpos2=params.sens_smear_pos_pre_GL,
        base_sensitivity_smearneg2=params.sens_smear_neg_pre_GL,
        treatment_fraction=0.30,   # was 0.8 — ~70% LTFU in private sector / low QoC
        pos_event='ProviderOrdersTBTest', pos_event2='ProviderTestNoR',
        start_day=params.genexpert_introduction, duration=-1,
        property_restrictions_list=['QualityOfCare:Low']
    ))

    # Resistance diagnostics
    for qoc, sens in [('High', params.sens_resistance_H), ('Low', params.sens_resistance_L)]:
        for trigger in ['ProviderOrdersTBTest', 'ProviderTestNoR']:
            camp.add(rd.ResistanceDiagnostic(
                camp=camp, trigger_treatment_list=[trigger],
                sensitivity=sens, specificity=params.specificity_resistance,
                neg_event='TBDS_Positive',
                treatment_fraction=0.5, treatment_fraction_negative_test=1.0,
                start_day=params.genexpert_introduction, duration=-1,
                property_restrictions_list=[f'QualityOfCare:{qoc}']
            ))

    # GeneXpert era treatment
    for qoc, drug in [('High', 'DOTSHQ'), ('Low', 'DOTSLQ')]:
        for triggers in [['TBTestPositive', 'TBDS_Positive'], ['TB_Tx']]:
            camp.add(treat.TBHIVDrugTreatment(
                camp=camp, trigger_treatment_list=triggers,
                drug_name=drug, latent_multiplier=0,
                start_day=params.genexpert_introduction, duration=-1,
                property_restrictions_list=[f'QualityOfCare:{qoc}']
            ))
    camp.add(treat.TBHIVDrugTreatment(
        camp=camp, trigger_treatment_list=['TBMDRTestPositive'],
        drug_name='DOTSMDR', latent_multiplier=0,
        start_day=params.genexpert_introduction, duration=-1,
        property_restrictions_list=[]
    ))

    # Cash on infection (for economic tracking)
    # Note: cash.triggered_cash was removed in emodpy-tbhiv 1.3.1; use cash.cash() instead
    cash_on_infection = cash.cash(camp, amount=personal_hc_budget)
    camp.add(cash_on_infection)

    return camp


# ══════════════════════════════════════════════════════════════════════════
# 3. DEMOGRAPHICS BUILDER
# ══════════════════════════════════════════════════════════════════════════

def build_mi_demog(hint_group_a_bi=3.1, hint_group_b_bi=2.9, hint_group_c_bi=3.5):
    """Build India-like demographics with quality-of-care stratification."""
    import emodpy_tbhiv.demographics.TBHIVDemographics as Demographics

    demog = Demographics.from_template_node(lat=0, lon=0, pop=args.pop, name=1, forced_id=1)

    # Quality of care: 15% high (urban/public), 85% low (rural/private)
    demog.AddIndividualPropertyAndHINT(
        "QualityOfCare", ["High", "Low", "None"], [0.15, 0.85, 0]
    )

    # Geographic prevalence groups with HINT transmission matrix
    demog.AddIndividualPropertyAndHINT(
        "Geographic", ["Low", "Mid", "High"], [0.333, 0.333, 0.334],
        [[hint_group_a_bi, 1, 1],
         [1, hint_group_b_bi, 1],
         [1, 1, hint_group_c_bi]]
    )

    # Age-dependent transmission (PREM India contact matrix)
    demog.AddAgeDependentTransmission(
        Age_Bin_Edges_In_Years=[0, 20, 40, 60, -1],
        TransmissionMatrix=[
            [3.5, 1.1, 0.5, 0.1],
            [1.1, 2.0, 0.8, 0.1],
            [0.8, 1.1, 0.8, 0.1],
            [0.4, 0.4, 0.3, 0.1]
        ]
    )

    # India-specific mortality and fertility
    demog.SetMortalityOverTimeFromData(manifest.india_mortality_data, params.base_year)
    demog.SetFertilityOverTimeFromParams(
        years_region1=110, years_region2=60,
        start_rate=0.2, inflection_rate=0.17, end_rate=0.07
    )

    return demog


# ══════════════════════════════════════════════════════════════════════════
# 4. SWEEP FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════

def update_sim_random_seed(simulation, value):
    simulation.task.config.parameters.Run_Number = value
    return {'Run_Number': value}


def map_sample_to_model_input(simulation, sample):
    """Map calibration subsample parameters to model inputs."""
    sample = copy.deepcopy(sample)

    # Care seeking
    fast_days = sample.pop('Care Seeking Fast Duration Days')
    slow_days = sample.pop('Care Seeking Slow Duration Days')
    set_seeking(simulation.task, slow_days, fast_days)

    # Duration and mortality
    dur = sample.pop('TB Duration')
    frac = sample.pop('Death Fraction')
    set_duration_no_hiv(simulation.task, death_fraction=frac, duration=dur)

    # Presymptomatic duration
    presymp_yrs = sample.pop('TB Presymptomatic Duration Years', 0.1)
    simulation.task.config.parameters.TB_Presymptomatic_Rate = 1.0 / (presymp_yrs * 365)

    # TBHIV duration
    hiv_dur = sample.pop('TBHIV Duration')
    set_duration_hiv(simulation.task, death_fraction=frac, duration=hiv_dur)

    # HIV multiplier
    hiv_pro = sample.pop('Primary HIV multiplier', 1.0)
    set_primary_hiv_pro(simulation.task, hiv_pro)

    # ART factor
    art_factor = sample.pop('ART Factor', 1.0)

    # CD4 strata
    cd4_vals = [
        sample.pop('CD4_aq_below_200_rel', 54.55),
        sample.pop('CD4_aq_below_200_rel', 54.55),
        sample.pop('CD4_aq_200_300_rel', 20.434),
        sample.pop('CD4_aq_300_400_rel', 3.0),
        sample.pop('CD4_aq_400_500_rel', 2.0),
        sample.pop('CD4_aq_above_500', 1.0),
        sample.pop('CD4_aq_above_500', 1.0),
    ]

    # Relative infectiousness
    sample.pop('Relative TB CD4 Infectiousness', 0.1)

    # Map remaining directly to config params
    for key, value in sample.items():
        if hasattr(simulation.task.config.parameters, key):
            setattr(simulation.task.config.parameters, key, value)

    return sample


def set_seeking(task, duration_slow, duration_fast):
    """Update campaign with new care-seeking durations."""
    loc_rate_slow = 1.0 / duration_slow
    loc_rate_fast = 1.0 / duration_fast
    build_camp_partial = partial(
        build_camp,
        probability_per_step_slow=loc_rate_slow,
        probability_per_step_fast=loc_rate_fast
    )
    task.create_campaign_from_callback(build_camp_partial)
    return {'LowSeek_Duration': duration_slow, 'Highseek_Duration': duration_fast}


# ══════════════════════════════════════════════════════════════════════════
# 5. DRUG CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════

def configure_drugs(task):
    """Set up all TB drug types."""
    import emodpy_tbhiv.standard_drugs as sd
    sd.bootstrap(task.config)

    add_tb_drug(task=task, drug_name='ACFDOTS')

    add_tb_drug_type(task=task, drug_name='DOTSMDR',
                     duration=270, cure_proportion=0.7, death_proportion=0.03,
                     resistance_proportion=0, inactivation_proportion=0.10,
                     relapse_proportion=0.02)

    add_tb_drug_type(task=task, drug_name='PreDOTSHigh',
                     duration=180, cure_proportion=0.7, death_proportion=0.03,
                     resistance_proportion=0, inactivation_proportion=0.10,
                     relapse_proportion=0.02)

    add_tb_drug_type(task=task, drug_name='PreDOTSLow',
                     duration=180, cure_proportion=0.6, death_proportion=0.03,
                     resistance_proportion=0, inactivation_proportion=0.10,
                     relapse_proportion=0.02)

    add_tb_drug_type(task=task, drug_name='Universal',
                     duration=90, cure_proportion=0.65, death_proportion=0.03,
                     resistance_proportion=0.02, inactivation_proportion=0.10,
                     relapse_proportion=0.02)


# ══════════════════════════════════════════════════════════════════════════
# 6. MAIN RUN FUNCTION
# ══════════════════════════════════════════════════════════════════════════

def run_simulation():
    """Build and run the TB India experiment on COMPS."""

    # Download EMOD executable
    print("Setting up EMOD binaries...")
    dtk.setup(pathlib.Path(manifest.eradication_path).parent)

    # Connect to COMPS platform
    if not args.dry_run:
        print(f"Connecting to {args.platform}...")
        platform = Platform(args.platform, priority='AboveNormal', node_group='idm_48cores')

    # Create EMODTask
    print("Creating EMODTask...")
    emod_task.dev_mode = True

    task = EMODTask.from_default2(
        eradication_path=manifest.eradication_path,
        campaign_builder=build_camp,
        schema_path=manifest.schema_file,
        param_custom_cb=set_param_fn,
        ep4_path=manifest.ep4_path,
        demog_builder=build_mi_demog,
        config_path='my_config.json'
    )

    # Add overlay assets
    print("Adding asset files...")
    task.common_assets.add_directory(assets_directory=manifest.assets_input_dir)

    # Additional config
    update_more_config(task)

    # Configure drugs
    configure_drugs(task)

    # Set singularity container
    task.set_sif(manifest.sif)

    if args.dry_run:
        print("\n=== DRY RUN: Config and campaign built successfully ===")
        print("Config parameters set. Campaign events created.")
        print("No COMPS submission.")
        return None

    if args.calibrate:
        # Calibration mode
        from idmtools_calibra import calib_manager
        resist = 0.0
        from add_tb_drug_type import add_tb_drug_type as add_drugs_calib_fn
        calib_manager.platform = platform

        sites = [IndiaCalibSite()]
        # ... calibration logic would go here
        print("Calibration mode - see original notebook for full calibration setup")
        return None
    else:
        # Base case run
        from idmtools.entities.templated_simulation import TemplatedSimulations
        from idmtools.entities.experiment import Experiment
        from idmtools.builders import SimulationBuilder

        # Fixed parameter values (from calibration results)
        subsample = {
            'ART Factor': 1.0,
            'Base_Infectivity_Constant': 0.017,
            'CD4_aq_200_300_rel': 20.434,
            'CD4_aq_300_400_rel': 3.0,
            'CD4_aq_400_500_rel': 2.0,
            'CD4_aq_above_500': 1.0,
            'CD4_aq_below_200_rel': 54.55,
            'Care Seeking Fast Duration Days': 30,
            'Care Seeking Slow Duration Days': 60,
            'Death Fraction': 0.7,
            'Post_Infection_Acquisition_Multiplier': 0.8,
            'Primary HIV multiplier': 1.0,
            'Relative TB CD4 Infectiousness': 0.1,
            'TB Duration': 1.8,
            'TB Presymptomatic Duration Years': 0.1,
            'TBHIV Duration': 1.82,
            'TB_Active_Presymptomatic_Infectivity_Multiplier': 0.1660,
            'TB_Fast_Progressor_Fraction_Adult': 0.19,
            'TB_Slow_Progressor_Rate': 1.65e-05,
        }

        # Build sweep
        builder = SimulationBuilder()
        builder.add_sweep_definition(update_sim_random_seed, range(args.replicates))
        builder.add_sweep_definition(map_sample_to_model_input, [subsample])

        # Create and run experiment
        experiment = Experiment.from_builder(
            builders=builder, base_task=task,
            name="TB India_Base_25yr"
        )

        print(f"Submitting experiment to {args.platform}...")
        experiment.run(wait_until_done=True, platform=platform)

        # Handle completion
        task.handle_experiment_completion(experiment)
        print(f"\nExperiment {experiment.uid} succeeded!")

        # Download results
        EMODTask.get_file_from_comps(experiment.uid, [
            "ReportEventRecorder.csv",
            "PropertyReportTB.json"
        ])

        EMODTask.cache_experiment_metadata_in_sql(experiment.uid)

        print(f"\nExperiment ID: {experiment.uid}")
        print("Results downloaded to local directory.")
        return experiment.uid


# ══════════════════════════════════════════════════════════════════════════
# ENTRY POINT
# ══════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    experiment_id = run_simulation()
    if experiment_id:
        print(f"\nDone! Experiment ID: {experiment_id}")
        # Save for later analysis
        with open("COMPS_ID", "w") as f:
            f.write(str(experiment_id))
