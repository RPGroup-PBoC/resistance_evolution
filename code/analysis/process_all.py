import numpy as np
import scipy.stats as st
import pandas as pd

import cmdstanpy
import arviz as az

#import iqplot
import bebi103

import bokeh.io
import bokeh.plotting
from bokeh.io import export_svgs

files = [
    "20230207_r1_plate1",
    #"20230208_r1_plate2",
    "20230302_r1_plate3",
    #"20230303_r1_plate4",
    #"20230308_r1_plate5",
    #"20230505_r2_plate1",
    #"20230509_r2_plate2",
    #"20230512_r2_plate3",
    #"20230519_r2_plate4",
    #"20230630_r2_plate5",
]

bebi103.stan.clean_cmdstan()
sm = cmdstanpy.CmdStanModel(stan_file='mixed_growth.stan')

def fit_mixed_model(file, skip_wells=[]):
    
    df = pd.read_csv("../../processing/plate_reader/" + file + "/growth_plate.csv")

    plot_list = []
    df_return = pd.DataFrame()
    for i, well in enumerate(df['well'].unique()):

        strain = df.loc[df['well'] == well, "strain"].values[0]
        pos_sel = df.loc[df['well'] == well, "pos_selection"].values[0]
        date = file.split('_')[0]

        print(strain, pos_sel)
        inds = ~np.isnan(np.log(df.loc[df['well'] == well, "OD600_norm"].values))
        data = {
            "t": df.loc[df['well'] == well, "time_min"].values[inds],
            "y": np.log(df.loc[df['well'] == well, "OD600_norm"].values)[inds],
            "N": len(df.loc[df['well'] == well, "time_min"].values[inds])
        }
        with bebi103.stan.disable_logging():
            samples = sm.sample(
                data=data,
                chains=4,
                iter_sampling=1000,
                iter_warmup=1000, 
            )
        samples = az.from_cmdstanpy(posterior=samples, posterior_predictive=['y_ppc'])
        samples_df = samples.posterior.to_dataframe()
        y_ppc = samples.posterior_predictive["y_ppc"].stack({"sample": ("chain", "draw")})
        y_ppc = y_ppc.transpose('sample', 'y_ppc_dim_0')
        df_sum = samples_df[['K', 'lambda_G', 'lambda_P', 'y0', 'f', 'sigma']].agg(np.mean)


        df_return = pd.concat([df_return, pd.DataFrame(
            data={"K": [df_sum['K']],
                  "lambda_G": [df_sum['lambda_G']],
                  "lambda_P": [df_sum['lambda_P']],
                  "y0": [df_sum['y0']],
                  "f": [df_sum['f']],
                  "sigma": [df_sum['sigma']],
                  "strain": [strain],
                  "tc": [pos_sel],
                  "date": [date]
                 })])

        p = bebi103.viz.predictive_regression(
                y_ppc,
                percentiles=[30, 50, 70, 99],
                #data=ell,
                samples_x=data['t'],
                data=np.transpose([data['t'], data['y']]),
                x_axis_label='time [min]',
                title=pos_sel
            )
        p.output_backend = "svg"
        plot_list.append(p)

    return plot_list, df_return

df_out = pd.DataFrame()
for file in files:
    print("Processing {}".format(file))
    p, df_fit = fit_mixed_model(file)
    export_svgs(bokeh.layouts.gridplot(p, ncols=2), filename="output/{}.svg".format(file))
    df_out = pd.concat([df_out, df_fit])
    df_fit.to_csv("{}_processed.csv".format(file))


    
