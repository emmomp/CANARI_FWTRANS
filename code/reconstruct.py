import sys
import glob
import calendar
sys.path.insert(0,'/users/emmomp/Python/ECCOv4-py')
import ecco_v4_py as ecco
import xarray as xr
import numpy as np
import dask
import utils as ut
import os.path
import matplotlib.pyplot as plt

ecco_grid=xr.open_dataset('~/data/orchestra/other_data/ECCO_r3_alt/ECCOv4r3_grid.nc')
rootdir='/users/emmomp/data/canari/experiments'
expt='fwd_26y'
fcname='horflux_fw_denm'
fclabel='FW Flux through Denmark Strait'

eyears=['2006','2014','2000']

conv_dir='../data_out/denm_X_ECCOclimanom_hfreq'
adj_diag_map={'adxx_qnet':['EXFqnet','oceQnet'],'adxx_tauu':['oceTAUU','EXFtauu'],'adxx_tauv':['oceTAUV','EXFtauv'],
                 'adxx_empmr':['EXFempmr','oceFWflx']
                }

ecco_convs={}
ecco_convs['all']=[]
ecco_convs['OCE']=[]
ecco_convs['EXF']=[]
ecco_convs_2d=[]
for x in adj_diag_map:
    ecco_convs['all']=ecco_convs['all']+[x+'X'+v+'_sum' for v in adj_diag_map[x]]
    ecco_convs_2d=ecco_convs_2d+[x+'X'+v for v in adj_diag_map[x]]
    ecco_convs['OCE']=ecco_convs['OCE']+[x+'X'+v+'_sum' for v in adj_diag_map[x] if v[:3]=='oce']
    ecco_convs['EXF']=ecco_convs['EXF']+[x+'X'+v+'_sum' for v in adj_diag_map[x] if v[:3]=='EXF']

oexps=['transfw_Mar_noparam_7d','transfw_Jun_noparam_7d','transfw_Sep_noparam_7d','transfw_Dec_noparam_7d']
mth=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
mthi=dict(zip(mth,list(range(1,13))))

print('Loading solution')
fc_climanom,fc_mth=ut.get_soln(fcname,rootdir+'/'+expt)

for eyear in eyears:
    print(f'Calculating reconstructions for {eyear}')
    
    cexps_full=glob.glob(f'{conv_dir}/transfw_*_7d_{eyear}')
    cexps=[exp.split('/')[-1] for exp in cexps_full]

    cexps_mdict={}
    for exp in cexps:
        cexps_mdict[exp]=mthi[exp.split('/')[-1].split('_')[1]]
    cexps_edict={m:k for k,m in cexps_mdict.items() }
    
    conv_ecco=[]
    for exp in cexps:
        print(f'Loading {exp}')
        ds_exp=[]
        for year in range(1992,2018):
            ds=xr.open_mfdataset(f'{conv_dir}/{exp}/{year}/*.nc',coords='minimal').assign_coords({'year':year})
            ds_exp.append(ds[ecco_convs['all']])
        ds_exp=xr.concat(ds_exp,'year').assign_coords({'exp':exp,'month':cexps_mdict[exp]}).chunk({'year':26,'lag_years':260})
        conv_ecco.append(ds_exp)
    conv_ecco=xr.concat(conv_ecco,'exp').sel(year=slice(1992,2018))

    with dask.config.set(**{'array.slicing.split_large_chunks': False}):
        conv_ecco=conv_ecco.sortby(conv_ecco.lag_years,ascending=False)

    plotdates=[]
    for ie,exp in enumerate(cexps):
        plotdates.append([np.datetime64(f'{conv_ecco.year[i].data}-{cexps_mdict[exp]:02.0f}-16','ns') for i in range(0,26)])

    dJpred_ecco=conv_ecco.sum(dim='lag_years').assign_coords(dates=(['exp','year'],plotdates))
    dJpred_ecco['wind_EXF']=dJpred_ecco['adxx_tauuXEXFtauu_sum']+dJpred_ecco['adxx_tauvXEXFtauv_sum']
    dJpred_ecco['wind_OCE']=dJpred_ecco['adxx_tauuXoceTAUU_sum']+dJpred_ecco['adxx_tauvXoceTAUV_sum']
    dJpred_ecco['all_OCE']=dJpred_ecco[ecco_convs['OCE']].to_array().sum('variable')
    dJpred_ecco['all_EXF']=dJpred_ecco[ecco_convs['EXF']].to_array().sum('variable')

    dJpred_ecco_cumsum=conv_ecco.cumsum('lag_years').assign_coords(dates=(['exp','year'],plotdates))
    dJpred_ecco_cumsum['wind_EXF']=dJpred_ecco_cumsum['adxx_tauuXEXFtauu_sum']+dJpred_ecco_cumsum['adxx_tauvXEXFtauv_sum']
    dJpred_ecco_cumsum['wind_OCE']=dJpred_ecco_cumsum['adxx_tauuXoceTAUU_sum']+dJpred_ecco_cumsum['adxx_tauvXoceTAUV_sum']
    dJpred_ecco_cumsum['all_OCE']=dJpred_ecco_cumsum[ecco_convs['OCE']].to_array().sum('variable')
    dJpred_ecco_cumsum['all_EXF']=dJpred_ecco_cumsum[ecco_convs['EXF']].to_array().sum('variable')
    dJpred_ecco_cumsum=dJpred_ecco_cumsum.assign_coords(conv_ecco.coords)

    dJ_vars=dJpred_ecco.squeeze().stack(yearexp=['exp','year']).sortby('dates')
    
    cum_ev=xr.open_dataset(f'../data_out/horflux_fw_denm_cumev_bylag_byvar_{eyear}.nc')
    cum_ev_bym=xr.open_dataset(f'../data_out/horflux_fw_denm_cumev_bylag_byvar_bymonth_{eyear}.nc')
    cum_ev_bym['lag_years']=cum_ev['lag_years']
    lagmax=cum_ev_bym.idxmax('lag_years').squeeze().load()
    lagmax_ds=xr.concat([lagmax[var] for var in lagmax],'var').assign_coords(var=[var for var in lagmax])
    lagmax_ds.name='lag_max'
    
    for var in ecco_convs['all']+['wind_OCE','wind_EXF','all_OCE','all_EXF']:
        print('Writing reconstructions')
        print(f'Full reconstruction, {var}')
        YY=dJ_vars[var].swap_dims({'yearexp':'dates'}).rename({'dates':'time'}).assign_coords({'eyear':eyear})
        YY.drop_vars(['yearexp','exp']).to_netcdf(f'../data_out/denstr_fwflux_4yrecon_{eyear}_{var}.nc')

        # Peak year reconstruction   
        lmax=cum_ev[var].idxmax('lag_years')
        dJpred_maxev=dJpred_ecco_cumsum[var].sel(lag_years=lmax,method='nearest').stack(yearexp=['exp','year']).sortby('dates')
        dJpred_maxev=dJpred_maxev.swap_dims({'yearexp':'dates'}).drop_vars('time').rename({'dates':'time'})
        print(f'{var} recon max {lmax:2.1f}y')

        YY=dJpred_maxev.assign_coords({'eyear':eyear})
        print(f'Peak reconstruction, {var}')
        YY.drop_vars(['yearexp','exp']).to_netcdf(f'../data_out/denstr_fwflux_peakEVrecon_{eyear}_{var}.nc')

        # Peak seasonal reconstruction      
        lmax=lagmax_ds.sortby('month').sel(var=var)
        dJpred_maxev=dJpred_ecco_cumsum[var].squeeze().swap_dims({'exp':'month'}).sortby('month').sel(lag_years=lmax,method='nearest').sortby('month')
        dJpred_maxev=dJpred_maxev.stack(yearmonth=['year','month']).drop_vars('time')

        YY=dJpred_maxev.swap_dims({'yearmonth':'dates'}).rename({'dates':'time'}).assign_coords({'eyear':eyear})
        print(f'Peak monthly reconstruction, {var}')
        YY.drop_vars(['yearmonth','exp']).to_netcdf(f'../data_out/denstr_fwflux_mthEVrecon_{eyear}_{var}.nc')
        
        print(f'Done {var}')