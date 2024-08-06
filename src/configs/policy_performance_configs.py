
objective_view_configs = [
    {
        'experiment': '2_new',
        'policy': 191,
        'objectives': ['Avg_Mead_PE', 'LB_Shortage_Volume', 'Mead_1000', 'Powell_3490', 'Powell_Release_LTEMP']
    },
    {
        'experiment': '4_new',
        'policy': 84,
        'objectives': ['Avg_Mead_PE', 'LB_Shortage_Volume', 'Mead_1000', 'Powell_3490', 'Powell_Release_LTEMP']
    }
]

acceptability_view_configs = [
    {
        'experiment': '2_new',
        'policy': 191,
        'objectives': ['Avg_Mead_PE', 'LB_Shortage_Volume'],
        'thresholds': [-1100, 2000000]
    },
    {
        'experiment': '4_new',
        'policy': 84,
        'objectives': ['Avg_Mead_PE', 'LB_Shortage_Volume'],
        'thresholds': [-1100, 2000000]
    }

]

