"""Drop-in replacement for the production run cell, sweeping NROY parameter sets.

WHY THIS IS NOT A ONE-LINE CHANGE
---------------------------------
The production notebook sweeps a single parameter set and publishes its
campaign parameters to a module-level global:

    PRODUCTION_CAMP_KWARGS = campaign_kwargs_from_sample(subsample)
    builder.add_sweep_definition(map_sample_to_model_input, [subsample])
    builder.add_sweep_definition(sweep_combined, COMBINED_SCENARIOS)

`sweep_combined` *replaces* the campaign registered by the first sweep, and
re-applies the global so the secular cascade trends are not lost.  With one
parameter set that is correct.  With twenty it is not: every simulation would
receive the campaign parameters of whichever sample happened to fill the
global, so the care-seeking delays, the treatment-initiation and care-seeking
gains, and the secular infectivity decline - five of the ten free parameters,
including the three the calibration constrains least - would be frozen while
the config-side parameters varied.  Nothing would fail; the run would simply
not be the propagation it claims to be.

The fix is to stop treating the parameter set and the scenario as independent
sweeps.  One sweep over their cross-product registers the campaign once,
carrying that simulation's own campaign parameters together with its scenario.

USAGE
-----
Run the production notebook as far as the task is built, then run this cell in
place of the sweep block at the end of run_experiment().  It expects the
notebook's own objects: build_camp, map_sample_to_model_input,
campaign_kwargs_from_sample, add_tb_drug_type, update_sim_random_seed,
COMBINED_SCENARIOS, TX_DRUG_NAME, TX_CURE_PROPORTION, and the production
`subsample` dict.
"""
import json
from functools import partial
from pathlib import Path

from idmtools.builders import SimulationBuilder
from idmtools.entities.experiment import Experiment

HERE = Path(__file__).resolve().parent
DATA = HERE / "derived"
SAMPLES = json.load(open(DATA / "nroy_parameter_sets.json", encoding="utf-8"))
FREE = ["Base_Infectivity_Constant", "TB Presymptomatic Duration Years",
        "Care Seeking Fast Duration Days", "Care Seeking Slow Duration Days",
        "TB_Slow_Progressor_Rate", "Post_Infection_Acquisition_Multiplier",
        "TB_Fast_Progressor_Fraction_Adult", "tf_gain_per_year",
        "cs_gain_per_year", "infectivity_decline_per_year"]
REPLICATES_PROPAGATION = 10


def build_propagation_experiment(task, subsample, build_camp,
                                 map_sample_to_model_input,
                                 campaign_kwargs_from_sample, add_tb_drug_type,
                                 update_sim_random_seed, COMBINED_SCENARIOS,
                                 TX_DRUG_NAME, TX_CURE_PROPORTION,
                                 replicates=REPLICATES_PROPAGATION):
    """Cross-product sweep over NROY parameter sets and intervention scenarios."""

    def sample_for(rec):
        """The production dict with only the ten calibrated columns replaced.

        Every other key - the HIV block, Death Fraction, TB Duration - is left
        at its production value, because those were fixed during calibration
        and are not part of the retained region.
        """
        s = dict(subsample)
        for k in FREE:
            s[k] = rec[k]
        return s

    combos = [(rec, scen) for rec in SAMPLES for scen in COMBINED_SCENARIOS]

    def sweep_param_and_scenario(simulation, value):
        rec, scen = value
        sample = sample_for(rec)
        acf_cov, dx_boost, tx_boost, label = scen

        # config-side parameters; this also registers a campaign, which the
        # scenario registration below deliberately replaces, exactly as in the
        # production run
        map_sample_to_model_input(simulation, sample)

        # the campaign parameters of THIS sample, not of a global
        simulation.task.create_campaign_from_callback(
            partial(build_camp,
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
        # tags must stay numeric: cache_experiment_metadata_in_sql calls float()
        return {"Scenario_Code": scenario_code,
                "Param_Set": int(rec["param_set_id"]),
                "Is_Production": int(rec["is_production"]),
                "ACF_Coverage": acf_cov,
                "Dx_Boost": int(dx_boost),
                "Tx_Boost": int(tx_boost)}

    builder = SimulationBuilder()
    builder.add_sweep_definition(update_sim_random_seed, range(replicates))
    builder.add_sweep_definition(sweep_param_and_scenario, combos)

    n = replicates * len(combos)
    print(f"parameter sets      : {len(SAMPLES)}")
    print(f"scenarios           : {len(COMBINED_SCENARIOS)}")
    print(f"replicates          : {replicates}")
    print(f"total simulations   : {n}")
    assert n == replicates * len(SAMPLES) * len(COMBINED_SCENARIOS)

    return Experiment.from_builder(builders=builder, base_task=task,
                                   name="TB_India_v28_NROY_propagation")
