
### 2-D Sweep Settings
replicates = 2
samples = 2
xrm_min = 0.1
xrm_max = 0.3
rural_inf_mult_min = 0.0
rural_inf_mult_max = 1.0

### Calibration Settings
guess = 0.01 # b_i_c
#guess = 1e-2
initial_res = 1.0
#search_range = 0.95
search_range = 3
#groupby_tag="TB_Fast_Progressor_Rate"
#groupby_tag="Base_Infectivity_Constant"


exp_name = "India dynamic"
nSims = 1
base_year = 1850
base_infectivity = 1.0
burn_initial = 2 * 365
burn_predots = 145 * 365     # India started DOTS 1997, scaled up to 2006

dots_start = burn_initial + burn_predots     # ~1997
hiv_epidemic_start = dots_start - 15 * 365
art_start = hiv_epidemic_start + 21 * 365
genexpert_introduction = dots_start + 15 * 365     # India CBNAAT/GeneXpert scale-up from 2012
length_gene_xpert_ramp = 8 * 365     # Scaled up from 2012-2020
To_end_from_DOTS = 41.0 * 365

vaccine_introduction = dots_start + 28 * 365     # Assume new vaccine from 2025

# passive seeking rates
high_seek = 1.0 / (3 * 30)
low_seek = 1.0 / (7 * 30)
b_offset = dots_start / 365.0 - 2000     # ie DOTS starting in 2000

# ART introduction in 2007
seek_200 = 1.0 / 30.0
seek_350_500 = 1.0 / 30.0

# sensitivities PreDOTS
sens_smear_pos_pre_GL = 0.6 #*0.5
sens_smear_neg_pre_GL = 0.6 #*0.5  # 0.3  # 0.7

sens_smear_pos_pre_GH = 0.7 #*0.5
sens_smear_neg_pre_GH = 0.7 #*0.5  # 0.3  # 0.7

sens_smear_pos_GXL = 0.65   # private sector: gradual improvement via X-ray/clinical dx
sens_smear_neg_GXL = 0.60   # private sector: improved smear-neg detection

sens_smear_pos_GXH = 0.75   # was 0.99 — GeneXpert with operational/specimen losses
sens_smear_neg_GXH = 0.65   # was 0.90 — smear-neg harder even with molecular Dx

sens_resistance_L = 1.0
sens_resistance_H = 1.0
specificity_resistance = 1.00
sens_resistance_L_clinical = 1.0
sens_resistance_H_clinical = 1.0
specificity_resistance_clinical = 1.00

CRP_Sensitivity = 0.91
CRP_Specificity = 0.59

intervention_day = dots_start + 22 * 365     # Assume new interventions start ~2022
dots_plus_15 = dots_start + 15.0 * 365.0

# Dynamic = False will set at Guess value
parameters = [
    {
        'Name': 'ART Factor',
        'Dynamic': False,
        # coinfection death scalar, 0 is like no HIV, 1 like HIV with no ART
        'Guess': 1.0,
        'Min': 0.0,
        'Max': 1.0
    },
    # {
    #     'Name': 'ART_Reactivation_Factor',
    #     # SC: This parameter is not listed in "subsample" dict
    #     'Dynamic': False,
    #     'MapTo': 'ART_Reactivation_Factor',
    #     'Guess': 1.0,
    #     'Min': 0.05,
    #     'Max': 1.0
    # },
    {
        'Name': 'Base_Infectivity_Constant',
        'Dynamic': True,
        'MapTo': 'Base_Infectivity_Constant',
        'Guess': 0.009,
        'Min': 0.005,
        'Max': 0.016
    },
    {
        'Name': 'CD4_aq_below_200_rel',
        'Dynamic': False,
        # first two
        'Guess': 54.5,
        'Min': 4.0,
        'Max': 60.0
    },
    {
        'Name': 'CD4_aq_200_300_rel',
        'Dynamic': False,
        # third
        'Guess': 4.44,
        'Min': 3.0,
        'Max': 50.0
    },
    {
        'Name': 'CD4_aq_300_400_rel',
        'Dynamic': False,
        # fourth
        'Guess': 2.0,
        'Min': 2.0,
        'Max': 50.0
    },
    {
        'Name': 'CD4_aq_400_500_rel',
        'Dynamic': False,
        # fifth
        'Guess': 1.0,
        'Min': 1.0,
        'Max': 23.0
    },
    {
        'Name': 'CD4_aq_above_500',
        'Dynamic': False,
        # six and 7th (for now 7th is really superfluous in future, but for compatibility)
        'Guess': 1.0,
        'Min': 1.0,
        'Max': 4.0
    },
    {
        'Name': 'Care Seeking Fast Duration Days',
        'Dynamic': True,
        # for now
        'Guess': 173.14,
        'Min': 30.0,
        'Max': 270.0
    },
    {
        'Name': 'Care Seeking Slow Duration Days',
        'Dynamic': True,
        # for now
        'Guess': 493.25,
        'Min': 90.0,
        'Max': 725.0
    },
    {
        'Name': 'Death Fraction',
        'Dynamic': False,
        # for now
        'Guess': 0.7,
        'Min': 0.4,
        'Max': 0.8
    },
    {
        'Name': 'Post_Infection_Acquisition_Multiplier',
        'Dynamic': True,
        'MapTo': 'Post_Infection_Acquisition_Multiplier',
        'Guess': 0.483,
        'Min': 0.3,
        'Max': 0.9
    },
    {
        'Name': 'Primary HIV multiplier',
        'Dynamic': False,
        # for now
        'Guess': 1.0,
        'Min': 1.0,
        'Max': 4.0
    },
    {
        'Name': 'Relative TB CD4 Infectiousness',
        'Dynamic': True,
        # for parsimony not CD4 or ART dependent
        'Guess': 0.1,
        'Min': 0.0,
        'Max': 0.9
    },
    {
        'Name': 'TB Duration',
        'Dynamic': True,
        # constrain to less than duration without HIV
        'Guess': 1.5,
        'Min': 1.0,
        'Max': 3.0
    },
    {
        'Name': 'TB_Presymptomatic_Duration_Years',
        'Dynamic': True,
        # 'MapTo': 'M'
        'Guess': 0.5,
        'Min': 10.0 / 365.0,
        'Max': 4.0
    },
    {
        'Name': 'TBHIV Duration',
        'Dynamic': True,
        # constrain to less than duration without HIV
        'Guess': 1.8,
        'Min': 1.0,
        'Max': 4
    },
    {
        'Name': 'TB_Active_Presymptomatic_Infectivity_Multiplier',
        'Dynamic': True,
        'MapTo': 'TB_Active_Presymptomatic_Infectivity_Multiplier',
        # This is based on a relative smear-negative of 0.34604 fixed (ok sig figs not withstanding) so
        # ranging from no infectivity to equivalent to a smear negative case at best
        'Guess': 0,
        'Min': 0.0,
        'Max': 1.0
    },
    {
        'Name': 'TB_Fast_Progressor_Fraction_Adult',
        'Dynamic': True,
        'MapTo': 'TB_Fast_Progressor_Fraction_Adult',
        'Guess': 0.15,
        'Min': 0.05,
        'Max': 0.20
    },
    {
        'Name': 'TB_Slow_Progressor_Rate',
        'Dynamic': True,
        'MapTo': 'TB_Slow_Progressor_Rate',
        'Guess': 1.75e-5,
        'Min': 0.001 / 365.0,
        'Max': 0.0075 / 365.0
    }
]