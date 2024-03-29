import seaborn as sns

VERBOSE = 0
palette = sns.color_palette('deep')
alpha = 0.8

"""
backgroundcolor and fontcolor: matplotlib named colors, rgb, or rgba tuples.
background_score_cmap: seaborn color palette name, or a list of rgb colors.
"""

config_default = {
    'CDS':{
        'backgroundcolor':palette[0] + (alpha,),
        'css_attributes':{}
    },
    'RNA':{
        'backgroundcolor':palette[2] + (alpha,),
        'css_attributes':{}
    },
    'coordinates':{
        'css_attributes':{'height':'100%',  # this is needed in order to color the empty elements in the row
                         'width':'100%'
                         }
    },
    'start_codon':{
        'fontcolor':'green',
        'css_attributes':{'font-weight':'bold'}
    },
    'stop_codon':{
        'fontcolor':'red',
        'css_attributes':{'font-weight':'bold'}
    },
    'RBS':{
        'css_attributes':{'text-decoration':'underline'}
    },
    'peptide_found':{
        'css_attributes':{'font-weight':'bold'}
    },
    'region':{
        'fontcolor':palette[2],
        'css_attributes':{'text-decoration':'underline'}
    },
    'region2':{
        'backgroundcolor':palette[3] + (0.5,)
    },
    'region_score':{
        'background_score_name':'score',
        # 'background_score_cmap':'Blues',
        'background_score_cmap':sns.light_palette((259,80,60), as_cmap=True, input='husl'),
        'background_score_vmin':0,
        'background_score_vmax':1
    }
}