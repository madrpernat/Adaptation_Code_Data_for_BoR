# No config for experiment 1 since it only has one policy, and HiPlot seems to need at least two observations to work

configs = [
    {
        'name': 'experiment2',
        'exp_directories': ['output/borg_experiments_analysis/TieredPowCoordOps_TieredMead_experiment2_seed8'],
        'objectives': ['Powell.3490', 'Powell.WY.Release', 'Mead.1000', 'Avg.LB.Shortage'],
        'output': 'output/borg_experiments_analysis/TieredPowCoordOps_TieredMead_experiment2_seed8'
    },
    {
        'name': 'experiment3',
        'exp_directories': ['output/borg_experiments_analysis/TieredPowCoordOps_TieredMead_experiment3_seed8'],
        'objectives': ['Powell.3490', 'Powell.WY.Release', 'Mead.1000', 'Avg.LB.Shortage'],
        'output': 'output/borg_experiments_analysis/TieredPowCoordOps_TieredMead_experiment3_seed8'
    },
    {
        'name': 'experiment4',
        'exp_directories': ['output/borg_experiments_analysis/TieredPowCoordOps_TieredMead_experiment4_seed8'],
        'objectives': ['Powell.3490', 'Powell.WY.Release', 'Mead.1000', 'Avg.LB.Shortage'],
        'output': 'output/borg_experiments_analysis/TieredPowCoordOps_TieredMead_experiment4_seed8'
    },
    {
        'name': 'all_original_experiments',
        'exp_directories': ['output/borg_experiments_analysis/TieredPowCoordOps_TieredMead_experiment1_seed8',
                            'output/borg_experiments_analysis/TieredPowCoordOps_TieredMead_experiment2_seed8',
                            'output/borg_experiments_analysis/TieredPowCoordOps_TieredMead_experiment3_seed8',
                            'output/borg_experiments_analysis/TieredPowCoordOps_TieredMead_experiment4_seed8'
                            ],
        'objectives': ['Powell.3490', 'Powell.WY.Release', 'Mead.1000', 'Avg.LB.Shortage'],
        'output': 'output/borg_experiments_analysis'
    },
    {
        'name': 'experiment2_new__orig_objectives',
        'exp_directories': ['output/borg_experiments_analysis/TieredPowCoordOps_TieredMead_experiment2_seed8_NEW'],
        'objectives': ['Powell.3490', 'Powell.WY.Release', 'Mead.1000', 'Avg.LB.Shortage'],
        'output': 'output/borg_experiments_analysis/TieredPowCoordOps_TieredMead_experiment2_seed8_NEW'
    },
    {
        'name': 'experiment4_new__orig_objectives',
        'exp_directories': ['output/borg_experiments_analysis/TieredPowCoordOps_TieredMead_experiment4_seed8_NEW'],
        'objectives': ['Powell.3490', 'Powell.WY.Release', 'Mead.1000', 'Avg.LB.Shortage'],
        'output': 'output/borg_experiments_analysis/TieredPowCoordOps_TieredMead_experiment4_seed8_NEW'
    },
    {
        'name': 'experiment2_new__new_objectives',
        'exp_directories': ['output/borg_experiments_analysis/TieredPowCoordOps_TieredMead_experiment2_seed8_NEW'],
        'objectives': ['Avg.Powell.PE', 'Powell.Release.LTEMP', 'Avg.Mead.PE', 'Avg.LB.Shortage'],
        'output': 'output/borg_experiments_analysis/TieredPowCoordOps_TieredMead_experiment2_seed8_NEW'
    },
    {
        'name': 'experiment4_new__new_objectives',
        'exp_directories': ['output/borg_experiments_analysis/TieredPowCoordOps_TieredMead_experiment4_seed8_NEW'],
        'objectives': ['Avg.Powell.PE', 'Powell.Release.LTEMP', 'Avg.Mead.PE', 'Avg.LB.Shortage'],
        'output': 'output/borg_experiments_analysis/TieredPowCoordOps_TieredMead_experiment4_seed8_NEW'
    },
    {
        'name': 'all_new_experiments__original_objectives',
        'exp_directories': ['output/borg_experiments_analysis/TieredPowCoordOps_TieredMead_experiment2_seed8_NEW',
                            'output/borg_experiments_analysis/TieredPowCoordOps_TieredMead_experiment4_seed8_NEW'
                            ],
        'objectives': ['Powell.3490', 'Powell.WY.Release', 'Mead.1000', 'Avg.LB.Shortage'],
        'output': 'output/borg_experiments_analysis'
    },
    {
        'name': 'all_new_experiments__new_objectives',
        'exp_directories': ['output/borg_experiments_analysis/TieredPowCoordOps_TieredMead_experiment2_seed8_NEW',
                            'output/borg_experiments_analysis/TieredPowCoordOps_TieredMead_experiment4_seed8_NEW'
                            ],
        'objectives': ['Avg.Powell.PE', 'Powell.Release.LTEMP', 'Avg.Mead.PE', 'Avg.LB.Shortage'],
        'output': 'output/borg_experiments_analysis'
    }
]
