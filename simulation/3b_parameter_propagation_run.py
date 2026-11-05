"""Propagation run: the production task, swept over 20 NROY parameter sets.

Generated from cell 25 of 3_production_scenarios.ipynb; the task-building code
is that cell's, unchanged.  Only the sweep differs - see the comment block
inside the function.

Usage, in a new notebook cell after cell 25 (do NOT run cell 26):

    exec(open("3b_parameter_propagation_run.py").read())
    uid = run_parameter_propagation()

Set REPLICATES_PROPAGATION before calling if 10 is not wanted.
"""
REPLICATES_PROPAGATION = 10

def run_parameter_propagation():
    from idmtools.entities.experiment import Experiment
    from idmtools.builders import SimulationBuilder

    # Create EMODTask
    print('Creating EMODTask...')
    emod_task.dev_mode = False

    task = EMODTask.from_default2(
        eradication_path=manifest.eradication_path,
        campaign_builder=build_camp,
        schema_path=manifest.schema_file,
        param_custom_cb=set_param_fn,
        ep4_path=manifest.ep4_path,
        demog_builder=build_mi_demog,
        config_path='my_config.json'
    )

    print('Adding asset files...')
    task.common_assets.add_directory(assets_directory=manifest.assets_input_dir)

    # Additional config
    update_more_config(task)

    # Bootstrap standard drugs
    import emodpy_tbhiv.standard_drugs as sd
    sd.bootstrap(task.config)

    # Custom drug definitions
    add_tb_drug(task=task, drug_name='ACFDOTS')
    add_tb_drug_type(task=task, drug_name='DOTSMDR', duration=270,
                     cure_proportion=0.7, death_proportion=0.03,
                     resistance_proportion=0, inactivation_proportion=0.10,
                     relapse_proportion=0.02)
    add_tb_drug_type(task=task, drug_name='PreDOTSHigh', duration=180,
                     cure_proportion=0.7, death_proportion=0.03,
                     resistance_proportion=0, inactivation_proportion=0.10,
                     relapse_proportion=0.02)
    add_tb_drug_type(task=task, drug_name='PreDOTSLow', duration=180,
                     cure_proportion=0.6, death_proportion=0.03,
                     resistance_proportion=0, inactivation_proportion=0.10,
                     relapse_proportion=0.02)
    add_tb_drug_type(task=task, drug_name='Universal', duration=90,
                     cure_proportion=0.65, death_proportion=0.03,
                     resistance_proportion=0.02, inactivation_proportion=0.10,
                     relapse_proportion=0.02)

    # Singularity container
    task.set_sif(manifest.sif)

    # ── Build sweep ───────────────────────────────────────────────
    subsample = {
        # ── Production point: round-2 combined fit, wave2 design 4 (experiment 496ded25) ──
        # max I = 0.27 on the six burden targets; 2024 treatment fraction 0.746,
        # care-seeking x1.56, infectivity x0.71 -- a LOOSE present-day cascade,
        # which is what leaves headroom for the Dx and Tx scenarios.
        # Keeps every key the calibration sweeps also saw; only the ten design
        # columns are replaced.
        'ART Factor': 1.0,
        'Base_Infectivity_Constant': 0.0419561921123124,   # round-2 calibration (wave2 design 4)
        'Care Seeking Fast Duration Days': 14.170451537838758,   # round-2 calibration (wave2 design 4)
        'Care Seeking Slow Duration Days': 157.77075650568554,   # round-2 calibration (wave2 design 4)
        'CD4_aq_200_300_rel': 20.434,
        'CD4_aq_300_400_rel': 3.0,
        'CD4_aq_400_500_rel': 2.0,
        'CD4_aq_above_500': 1.0,
        'CD4_aq_below_200_rel': 54.55,
        'cs_gain_per_year': 0.0164967327813104,   # round-2 calibration (wave2 design 4)
        'Death Fraction': 0.55,   # fixed; present throughout calibration, so it must stay here too
        'infectivity_decline_per_year': 0.0128311263239953,   # NEW in round 2: secular non-cascade decline (NodeInfectivityMult), 2024 multiplier x0.71
        'Post_Infection_Acquisition_Multiplier': 0.3197498251614573,   # round-2 calibration (wave2 design 4)
        'Primary HIV multiplier': 1.0,
        'Relative TB CD4 Infectiousness': 0.1,
        'sens_smear_pos_GXL': 0.65,
        'TB Duration': 3.5,   # FIXED, not calibrated: applied only together with 'Death Fraction', which no calibration design carries, so the swept range never reached the model
        'TB Presymptomatic Duration Years': 0.2330436272992007,   # round-2 calibration (wave2 design 4)
        'TB_Active_Presymptomatic_Infectivity_Multiplier': 0.166,
        'TB_Fast_Progressor_Fraction_Adult': 0.0764744993503116,   # round-2 calibration (wave2 design 4)
        'TB_Slow_Progressor_Rate': 8.265850430929792e-06,   # round-2 calibration (wave2 design 4)
        'TBHIV Duration': 1.5,
        'tf_gain_per_year': 0.0017202096282671,   # round-2 calibration (wave2 design 4)
        'treatment_fraction_LQ': 0.55,
    }
    # ====================================================================
    # PROPAGATION SWEEP (replaces the single-parameter-set sweep)
    #
    # The production version does:
    #     PRODUCTION_CAMP_KWARGS = campaign_kwargs_from_sample(subsample)
    #     builder.add_sweep_definition(map_sample_to_model_input, [subsample])
    #     builder.add_sweep_definition(sweep_combined, COMBINED_SCENARIOS)
    # which is correct for one parameter set and wrong for twenty: the scenario
    # callback replaces the campaign and re-applies a single global, so the five
    # campaign-side parameters would be frozen while the config-side ones varied.
    # Sweeping the cross product instead lets each simulation carry its own.
    # ====================================================================
    import json as _json
    from functools import partial as _partial

    _FREE = ["Base_Infectivity_Constant", "TB Presymptomatic Duration Years",
             "Care Seeking Fast Duration Days", "Care Seeking Slow Duration Days",
             "TB_Slow_Progressor_Rate", "Post_Infection_Acquisition_Multiplier",
             "TB_Fast_Progressor_Fraction_Adult", "tf_gain_per_year",
             "cs_gain_per_year", "infectivity_decline_per_year"]
    _SAMPLES = _json.load(open("../analysis/derived/nroy_parameter_sets.json", encoding="utf-8"))

    def _sample_for(rec):
        """The production dict with only the ten calibrated columns replaced."""
        s = dict(subsample)
        for k in _FREE:
            s[k] = rec[k]
        return s

    def sweep_param_and_scenario(simulation, value):
        rec, scen = value
        sample = _sample_for(rec)
        acf_cov, dx_boost, tx_boost, label = scen

        # config-side parameters; this also registers a campaign, which the
        # registration below deliberately replaces, exactly as in production
        map_sample_to_model_input(simulation, sample)

        simulation.task.create_campaign_from_callback(
            _partial(build_camp,
                     acf_coverage=acf_cov,
                     dx_boost=dx_boost,
                     tx_boost=tx_boost,
                     **campaign_kwargs_from_sample(sample)))

        if tx_boost:
            add_tb_drug_type(task=simulation.task, drug_name=TX_DRUG_NAME,
                             duration=180.0, cure_proportion=TX_CURE_PROPORTION,
                             death_proportion=0.03, resistance_proportion=0,
                             inactivation_proportion=0.10, relapse_proportion=0.02,
                             mdr_cure_proportion=0.1)

        scenario_code = (int(bool(acf_cov)) * 4 + int(bool(dx_boost)) * 2
                         + int(bool(tx_boost)))
        # tags stay numeric: cache_experiment_metadata_in_sql calls float()
        return {"Scenario_Code": scenario_code,
                "Param_Set": int(rec["param_set_id"]),
                "Is_Production": int(rec["is_production"]),
                "ACF_Coverage": acf_cov,
                "Dx_Boost": int(dx_boost),
                "Tx_Boost": int(tx_boost)}

    _combos = [(rec, scen) for rec in _SAMPLES for scen in COMBINED_SCENARIOS]

    builder = SimulationBuilder()
    builder.add_sweep_definition(update_sim_random_seed, range(REPLICATES_PROPAGATION))
    builder.add_sweep_definition(sweep_param_and_scenario, _combos)

    exp_name = "TB_India_v28_NROY_propagation"
    _n = REPLICATES_PROPAGATION * len(_combos)
    print(f"  parameter sets    : {len(_SAMPLES)}")
    print(f"  scenarios         : {len(COMBINED_SCENARIOS)}")
    print(f"  replicates        : {REPLICATES_PROPAGATION}")
    print(f"  TOTAL SIMULATIONS : {_n}")

    print(f'\nExperiment: {exp_name}')

    experiment = Experiment.from_builder(
        builders=builder, base_task=task, name=exp_name)
    experiment.run(wait_until_done=True, platform=platform)
    task.handle_experiment_completion(experiment)

    print(f'\nExperiment {experiment.uid} succeeded.')

    # Download outputs
    EMODTask.get_file_from_comps(experiment.uid, [
        'InsetChart.json',
        'ReportEventRecorder.csv',
        'Report_TBHIV_ByAge.csv',
        'PropertyReportTB.json',
        'config.json'
    ])
    # v8: cache step wrapped — non-critical, fails on non-numeric tags
    try:
        EMODTask.cache_experiment_metadata_in_sql(experiment.uid)
    except Exception as e:
        print(f'[v8] cache_experiment_metadata_in_sql skipped: {e}')

    return experiment.uid
