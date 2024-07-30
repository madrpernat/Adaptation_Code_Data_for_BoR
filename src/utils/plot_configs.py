from src.utils import ids


PLOT_TITLES = {
    ids.MEDIAN_FLOW: 'Median Annual Lees Ferry Flow [MAF]',
    ids.DEMAND: 'Annual Upper Basin Demand [MAF]',
    ids.INIT_STORAGE: 'Initial Combined Storage [Powell + Mead]',
    ids.MAX_ANNUAL_FLOW: 'Max Annual Lees Ferry Flow [MAF]',
    ids.MIN_ANNUAL_FLOW: 'Min Annual Lees Ferry Flow [MAF]',
    ids.WETTEST_10_YEAR_FLOW: 'Wettest 10-Year Average Annual Lees Ferry Flow [MAF]',
    ids.DRIEST_10_YEAR_FLOW: 'Driest 10-Year Average Annual Lees Ferry Flow [MAF]',
    ids.IQR_FLOW: 'Inter-Quartile Range of Annual Lees Ferry Flow [MAF]'
}


COLOR_SCHEMES = {
    ids.MEDIAN_FLOW: 'wet_dry',
    ids.DEMAND: 'warm_cool',
    ids.INIT_STORAGE: 'warm_cool',
    ids.MAX_ANNUAL_FLOW: 'wet_dry',
    ids.MIN_ANNUAL_FLOW: 'wet_dry',
    ids.WETTEST_10_YEAR_FLOW: 'wet_dry',
    ids.DRIEST_10_YEAR_FLOW: 'wet_dry',
    ids.IQR_FLOW: 'warm_cool'
}


INVERSE_COLORBAR = {
    ids.MEDIAN_FLOW: True,
    ids.DEMAND: False,
    ids.INIT_STORAGE: True,
    ids.MAX_ANNUAL_FLOW: True,
    ids.MIN_ANNUAL_FLOW: True,
    ids.WETTEST_10_YEAR_FLOW: True,
    ids.DRIEST_10_YEAR_FLOW: True,
    ids.IQR_FLOW: True
}
