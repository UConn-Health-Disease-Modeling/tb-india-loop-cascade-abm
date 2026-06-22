from numpy import log as log


def add_tb_drug_type(task, drug_name, duration, cure_proportion, death_proportion,
                     resistance_proportion, inactivation_proportion, relapse_proportion,
                     reduced_transmit=1, reduced_acquire=1,
                     cure_proportion_hiv=None, death_proportion_hiv=None,
                     resistance_proportion_hiv=None, relapse_proportion_hiv=None,
                     inactivation_proportion_hiv=None, mdr_cure_proportion=None,
                     mdr_death_proportion=None, mdr_inactivation_proportion=None,
                     mdr_relapse_proportion=None):
    if cure_proportion_hiv is None:
        cure_proportion_hiv = cure_proportion
    if inactivation_proportion_hiv is None:
        inactivation_proportion_hiv = inactivation_proportion
    if relapse_proportion_hiv is None:
        relapse_proportion_hiv = relapse_proportion
    if death_proportion_hiv is None:
        death_proportion_hiv = death_proportion
    if resistance_proportion_hiv is None:
        resistance_proportion_hiv = resistance_proportion

    if mdr_cure_proportion is None:
        mdr_cure_proportion = cure_proportion
    if mdr_inactivation_proportion is None:
        mdr_inactivation_proportion = inactivation_proportion
    if mdr_relapse_proportion is None:
        mdr_relapse_proportion = relapse_proportion
    if mdr_death_proportion is None:
        mdr_death_proportion = death_proportion

    rate_res = converttorate(resistance_proportion, duration)
    rate_res_hiv = converttorate(resistance_proportion_hiv, duration)

    total_rate = converttorate(cure_proportion + inactivation_proportion + death_proportion + relapse_proportion,
                               duration)
    rate_cure = cure_proportion * total_rate
    rate_relapse = relapse_proportion * total_rate
    rate_death = death_proportion * total_rate
    rate_inactivation = inactivation_proportion * total_rate

    total_ratehiv = converttorate(cure_proportion_hiv + inactivation_proportion_hiv +
                                  death_proportion_hiv + relapse_proportion_hiv, duration)
    rate_curehiv = cure_proportion_hiv * total_ratehiv
    rate_relapsehiv = relapse_proportion_hiv * total_ratehiv
    rate_deathhiv = death_proportion_hiv * total_ratehiv
    rate_inactivationhiv = inactivation_proportion_hiv * total_ratehiv

    total_rate_mdr = converttorate(mdr_cure_proportion + mdr_inactivation_proportion +
                                   mdr_death_proportion + mdr_relapse_proportion, duration)
    rate_curemdr = mdr_cure_proportion * total_rate_mdr
    rate_relapsemdr = mdr_relapse_proportion * total_rate_mdr
    rate_deathmdr = mdr_death_proportion * total_rate_mdr
    rate_inactivationmdr = mdr_inactivation_proportion * total_rate_mdr

    drug_dict = {drug_name: dict(TB_Drug_Cure_Rate=rate_cure, TB_Drug_Cure_Rate_MDR=rate_curemdr,
                                 TB_Drug_Cure_Rate_HIV=rate_curehiv, TB_Drug_Inactivation_Rate=rate_inactivation,
                                 TB_Drug_Inactivation_Rate_MDR=rate_inactivationmdr,
                                 TB_Drug_Inactivation_Rate_HIV=rate_inactivationhiv,
                                 TB_Drug_Mortality_Rate=rate_death, TB_Drug_Mortality_Rate_MDR=rate_deathmdr,
                                 TB_Drug_Mortality_Rate_HIV=rate_deathhiv, TB_Drug_Relapse_Rate=rate_relapse,
                                 TB_Drug_Relapse_Rate_MDR=rate_relapsemdr, TB_Drug_Relapse_Rate_HIV=rate_relapsehiv,
                                 TB_Drug_Resistance_Rate=rate_res, TB_Drug_Resistance_Rate_HIV=rate_res_hiv,
                                 TB_Drug_Primary_Decay_Time_Constant=duration, TB_Reduced_Acquire=reduced_acquire,
                                 TB_Reduced_Transmit=reduced_transmit)
                 }

    current_drug_params = task.config.parameters.TBHIV_Drug_Params
    current_drug_params[drug_name] = drug_dict[drug_name]
    task.config.parameters.TBHIV_Drug_Params = current_drug_params

    current_drug_names = task.config.parameters.TBHIV_Drug_Types
    if drug_name not in current_drug_names:
        current_drug_names.append(drug_name)
        task.config.parameters.TBHIV_Drug_Types = current_drug_names


def converttorate(proportion, interval):
    return -log(1.0 - proportion) / interval



def add_tb_drug(task, drug_name, rate_cure=None, rate_curehiv=None, rate_curemdr=None, rate_inactivation=None,
                rate_inactivationhiv=None, rate_inactivationmdr=None, rate_death=None, rate_deathhiv=None,
                rate_deathmdr=None, duration=None, rate_relapse=None, rate_relapsehiv=None,
                rate_relapsemdr=None, rate_res=None, reduced_acquire=1, reduced_transmit=1):
    rate_cure = 0.0007701635339554948
    rate_curehiv = 0.0007701635339554948
    rate_curemdr = 0.0007701635339554948
    rate_inactivation = 0.0007701635339554948
    rate_inactivationhiv = 0.0007701635339554948
    rate_inactivationmdr = 0.0007701635339554948
    rate_death = 0.0007701635339554948
    rate_deathhiv = 0.0007701635339554948
    rate_deathmdr = 0.0007701635339554948
    duration = 90
    rate_relapse = 0.0015403270679109895
    rate_relapsehiv = 0.0015403270679109895
    rate_relapsemdr = 0.0015403270679109895
    rate_res = 0.0011706723961980698
    rate_res_hiv = 0.0011706723961980698
    reduced_acquire = 1
    reduced_transmit = 1

    drug_dict = {drug_name: dict(TB_Drug_Cure_Rate=rate_cure, TB_Drug_Cure_Rate_MDR=rate_curemdr,
                                 TB_Drug_Cure_Rate_HIV=rate_curehiv, TB_Drug_Inactivation_Rate=rate_inactivation,
                                 TB_Drug_Inactivation_Rate_MDR=rate_inactivationmdr,
                                 TB_Drug_Inactivation_Rate_HIV=rate_inactivationhiv,
                                 TB_Drug_Mortality_Rate=rate_death, TB_Drug_Mortality_Rate_MDR=rate_deathmdr,
                                 TB_Drug_Mortality_Rate_HIV=rate_deathhiv, TB_Drug_Relapse_Rate=rate_relapse,
                                 TB_Drug_Relapse_Rate_MDR=rate_relapsemdr, TB_Drug_Relapse_Rate_HIV=rate_relapsehiv,
                                 TB_Drug_Resistance_Rate=rate_res, TB_Drug_Resistance_Rate_HIV=rate_res_hiv,
                                 TB_Drug_Primary_Decay_Time_Constant=duration, TB_Reduced_Acquire=reduced_acquire,
                                 TB_Reduced_Transmit=reduced_transmit)}

    current_drug_params = task.config.parameters.TBHIV_Drug_Params
    current_drug_params[drug_name] = drug_dict[drug_name]
    task.config.parameters.TBHIV_Drug_Params = current_drug_params

    current_drug_names = task.config.parameters.TBHIV_Drug_Types
    if drug_name not in current_drug_names:
        current_drug_names.append(drug_name)
        task.config.parameters.TBHIV_Drug_Types = current_drug_names
