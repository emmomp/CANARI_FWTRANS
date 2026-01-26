import xarray as xr
import numpy as np
import sys
sys.path.insert(0,'/users/emmomp/Python/ECCOv4-py')
import ecco_v4_py as ecco

root_dir='/data/smurphs/emmomp/canari/experiments/'
grid_dir='/data/smurphs/emmomp/orchestra/grid2/'
ecco_grid=xr.open_dataset('/data/expose/ECCOv4-r4/Version4/Release4/nctiles_grid/ECCO-GRID.nc') 

section='Denmark Strait'
section_label='Den'
[section_pt1,section_pt2]=ecco.get_section_endpoints(section)
masks_loaded=False

exps=['pert_10y_tauu_NGlandJANpulseplus','pert_10y_tauu_NGlandJANpulseminus']
#exps=['fwd_10y_1996_ctrl']
#exps=['pert_10y_tauu_NAlaskaSEPpulseplus','pert_10y_tauu_NAlaskaSEPpulseminus']

for exp in exps:
    print(f'Calculating {exp} {section} transports')
    
    ds=[]
  #  for var in ['SALT','THETA','ADVx_SLT','ADVy_SLT','DFxE_SLT','DFyE_SLT','ADVx_TH','ADVy_TH','DFxE_TH','DFyE_TH','UVELMASS','VVELMASS','GM_PsiX','GM_PsiY']:
    for var in ['SALT','THETA','UVELMASS','VVELMASS','GM_PsiX','GM_PsiY']:
        ds.append(ecco.load_ecco_vars_from_mds(f'{root_dir}{exp}/diags/{var}_mon_mean',mds_grid_dir=grid_dir,mds_files=f'{var}_mon_mean',
                                               #output_freq_code='AVG_MON',model_start_datetime=np.datetime64('1996-01-01T12:00:00'),read_grid=False))
                                                output_freq_code='AVG_MON',model_start_datetime=np.datetime64('1992-01-01T12:00:00'),read_grid=False))
    ds.append(ecco_grid)
    ds=xr.merge(ds).load()
    if not masks_loaded:
        line_maskC, line_maskW, line_maskS = ecco.get_section_line_masks(section_pt1,section_pt2,ds)
        masks_loaded=True
    grid = ecco.get_llc_grid(ds)
#    full_heat_transport=ecco.calc_section_heat_trsp(ds,maskC=line_maskC,along_section=True,grid=grid)
    full_vol_transport=ecco.calc_section_vol_trsp(ds,maskC=line_maskC,along_section=True,grid=grid)
    full_fw_transport=ecco.calc_section_fw_trsp(ds,Sref=35,maskC=line_maskC,grid=grid,along_section=True)
#    full_salt_transport=ecco.calc_section_salt_trsp(ds,maskC=line_maskC,along_section=True,grid=grid)
    full_transport=xr.merge([
       # full_heat_transport,
        full_vol_transport,
        full_fw_transport,
     #   full_salt_transport
    ])
    full_transport=full_transport.assign_coords({'exp':exp})
   # full_transport['fw_trsp']=full_transport['fw_trsp_adv']+full_transport['fw_trsp_dif']
    print('Writing to file')
    full_transport.reset_index('ij').to_netcdf(f'../data_out/{exp}_{section_label}_2dtransports.nc')
print('Done')
