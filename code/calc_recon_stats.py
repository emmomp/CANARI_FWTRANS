import xarray as xr
import numpy as np
import utils as ut
from scipy import signal, stats
from scipy.stats import chi2
import glob
from inputs import DATA_DIR, CONTR_DIR, ecco_convs, eyears, EV_DIR, FCNAME
import montecarlo_function as mont

var_list = ecco_convs + ["wind_OCE", "all_OCE"]
conv_var_list = [x.removesuffix("_sum") for x in var_list]
masks_plot=['global','gland','natl','arct','norw']
rtypes = ["full", "opt", "opt_lo", "opt_hi"]


def preproc(ds):
    ds = ds.interp(lag_years=np.arange(-5, 0, 0.025))
    return ds


ds_tseries_all = xr.open_dataset(
    f"{DATA_DIR}/horflux_fw_denm_adjsens_basin_tseries_vlag.nc"
)
dJ_all = []
for eyear in eyears:
    print(eyear)
    foo = xr.open_mfdataset(
        f"{CONTR_DIR}/{eyear}/horflux_fw_denm_contr_tseries_*.nc",
        combine="nested",
        concat_dim="exp",
        preprocess=preproc,
    )
    dJ_all.append(
        foo.assign_coords({"eyear": eyear})
        .swap_dims({"exp": "month"})
        .drop_vars(["exp"])
    )
dJ_all = xr.concat(dJ_all, "eyear", coords="minimal", compat="override")


dJopt_all = []
dJful_all = []
for eyear in eyears:
    print(eyear)
    files = glob.glob(f"{CONTR_DIR}/{eyear}/horflux_fw_denm_contr_tseries_*.nc")
    cum_ev = xr.open_dataset(
        f"{EV_DIR}/{eyear}/horflux_fw_denm_cumev_bylag_byvar_{eyear}.nc"
    )
    cum_ev_bym = xr.open_dataset(
        f"{EV_DIR}/{eyear}/horflux_fw_denm_cumev_bylag_byvar_bymonth_{eyear}.nc"
    )
    cum_ev_bym["lag_years"] = cum_ev["lag_years"]
    lagmax = cum_ev_bym.idxmax("lag_years").squeeze().load()
    dJopt_year = []
    dJful_year = []
    for file in files:
        foo = xr.open_dataset(file).sel(stat="sum").squeeze()
        dJ_var = []
        for var in var_list:
            lmax = lagmax[var].sel(month=foo.month)
            dJ_var.append(
                foo[var.removesuffix("_sum")]
                .sel(lag_years=lmax, method="nearest")
                .swap_dims({"year": "dates"})
            )
        dJ_var = xr.merge(dJ_var, compat="override")
        dJopt_year.append(dJ_var)
        dJful_year.append(
            foo[conv_var_list].isel(lag_years=-1).swap_dims({"year": "dates"})
        )
    dJopt_year = (
        xr.concat(dJopt_year, "dates").sortby("dates").assign_coords({"eyear": eyear})
    )
    dJful_year = (
        xr.concat(dJful_year, "dates").sortby("dates").assign_coords({"eyear": eyear})
    )
    dJopt_all.append(dJopt_year)
    dJful_all.append(dJful_year)
dJopt_all = xr.concat(dJopt_all, "eyear", coords="minimal", compat="override")
dJful_all = xr.concat(dJful_all, "eyear", coords="minimal", compat="override")

dJ_all_gma = (dJopt_all.sel(mask="global") - dJopt_all.sel(mask="arct")).assign_coords(
    {"mask": "global-arct"}
)
dJopt_all = xr.concat([dJopt_all, dJ_all_gma], "mask")
dJ_all_gma = (dJful_all.sel(mask="global") - dJful_all.sel(mask="arct")).assign_coords(
    {"mask": "global-arct"}
)
dJful_all = xr.concat([dJful_all, dJ_all_gma], "mask")

fc = ut.get_soln(FCNAME, DATA_DIR)
fc_climanom, fc_mth = ut.soln_anoms(fc)
XX = fc_climanom
XX_smooth = ut.butter_ufunc(XX, 13, "time")
XX_high = XX - XX_smooth

YY_opt = dJopt_all.sel(mask=masks_plot).drop_vars("time")
YY_full = dJful_all.sel(mask=masks_plot).drop_vars("time")
YY_opt_lo = ut.butter_ufunc(YY_opt, 13, "dates")
YY_opt_hi = YY_opt - YY_opt_lo

YYopt_dt = (
    YY_opt.to_array()
    .copy(data=signal.detrend(YY_opt.to_array()))
    .to_dataset("variable")
    .rename({"dates": "time"})
)
YYopthi_dt = (
    YY_opt_hi.to_array()
    .copy(data=signal.detrend(YY_opt_hi.to_array()))
    .to_dataset("variable")
    .rename({"dates": "time"})
)
YYoptlo_dt = (
    YY_opt_lo.to_array()
    .copy(data=signal.detrend(YY_opt_lo.to_array()))
    .to_dataset("variable")
    .rename({"dates": "time"})
)
YYfull_dt = (
    YY_full.to_array()
    .copy(data=signal.detrend(YY_full.to_array()))
    .to_dataset("variable")
    .rename({"dates": "time"})
)

ev_opt = (1 - (XX - YYopt_dt).var("time") / XX.var("time")).load()
ev_full = (1 - (XX - YYfull_dt).var("time") / XX.var("time")).load()
ev_opt_hi = (1 - (XX_high - YYopthi_dt).var("time") / XX_high.var("time")).load()
ev_opt_lo = (1 - (XX_smooth - YYoptlo_dt).var("time") / XX_smooth.var("time")).load()

evopt_mean = (1 - (XX - YYopt_dt.mean("eyear")).var("time") / XX.var("time")).load()
evfull_mean = (1 - (XX - YYfull_dt.mean("eyear")).var("time") / XX.var("time")).load()
evopthi_mean = (
    1 - (XX_high - YYopthi_dt.mean("eyear")).var("time") / XX_high.var("time")
).load()
evoptlo_mean = (
    1 - (XX_smooth - YYoptlo_dt.mean("eyear")).var("time") / XX_smooth.var("time")
).load()

ev_opt = xr.concat([ev_opt, evopt_mean.assign_coords({"eyear": "ens_mean"})], "eyear")
ev_full = xr.concat(
    [ev_full, evfull_mean.assign_coords({"eyear": "ens_mean"})], "eyear"
)
ev_opt_hi = xr.concat(
    [ev_opt_hi, evopthi_mean.assign_coords({"eyear": "ens_mean"})], "eyear"
)
ev_opt_lo = xr.concat(
    [ev_opt_lo, evoptlo_mean.assign_coords({"eyear": "ens_mean"})], "eyear"
)


YY_comb = xr.concat(
    [YYfull_dt.wind_OCE, YYopt_dt.wind_OCE, YYoptlo_dt.wind_OCE, YYopthi_dt.wind_OCE],
    "recon_type",
).assign_coords({"recon_type": rtypes})
XX_comb = xr.concat([XX, XX, XX_smooth, XX_high], "recon_type").assign_coords(
    {"recon_type": rtypes}
)

ev_p_comb = []

for rtype in rtypes:
    print(rtype)
    ev_aa = xr.DataArray(
        dims=["mask", "eyear"],
        coords={
            "mask": (("mask", masks_plot)),
            "eyear": ("eyear", eyears + ["ens_mean"]),
        },
    )
    ev_aa_up = xr.zeros_like(ev_aa)
    ev_aa_lo = xr.zeros_like(ev_aa)
    pp = xr.zeros_like(ev_aa)

    for im, mask in enumerate(masks_plot):
        print(mask)
        for iy, eyear in enumerate(eyears):
            series1 = XX_comb.sel(recon_type=rtype)
            series2 = YY_comb.sel(mask=mask, eyear=eyear, recon_type=rtype)
            a = series2.var() / series1.var()
            pr = stats.pearsonr(series1, series2)
            pr_ci = pr.confidence_interval()
            ev_aa[im, iy] = -a + 2 * np.sqrt(a) * pr.statistic
            ev_aa_lo[im, iy] = -a + 2 * np.sqrt(a) * pr_ci.low
            ev_aa_up[im, iy] = -a + 2 * np.sqrt(a) * pr_ci.high
            pp[im, iy] = mont.MonteCarlo_EV(series1.data, series2.data, ev_aa[im,iy].data)

        iy += 1
        series2 = YY_comb.sel(mask=mask, recon_type=rtype).mean("eyear")
        a = series2.var() / series1.var()
        pr = stats.pearsonr(series1, series2)
        pr_ci = pr.confidence_interval()
        ev_aa[im, iy] = -a + 2 * np.sqrt(a) * pr.statistic
        ev_aa_lo[im, iy] = -a + 2 * np.sqrt(a) * pr_ci.low
        ev_aa_up[im, iy] = -a + 2 * np.sqrt(a) * pr_ci.high
        pp[im, iy] = mont.MonteCarlo_EV(series1.data, series2.data, ev_aa[im,iy].data)
    ev_aa.name = "EV"
    ev_aa_lo.name = "EV_lower"
    ev_aa_up.name = "EV_upper"
    ev_p_comb.append(xr.merge([ev_aa, ev_aa_lo, ev_aa_up]))
ev_p_comb = xr.concat(ev_p_comb, "recon_type").assign_coords({"recon_type": rtypes})
ev_p_comb.to_netcdf(f"{DATA_DIR}/recon_evstats.nc")

stats_comb = []

for rtype in rtypes:
    print(rtype)
    pp = xr.DataArray(
        dims=["mask", "eyear"],
        coords={
            "mask": (("mask", masks_plot)),
            "eyear": ("eyear", eyears + ["combined", "ens_mean"]),
        },
    )
    pp_alt = xr.zeros_like(pp)
    rr = xr.zeros_like(pp)
    rr_upper = xr.zeros_like(pp)
    rr_lower = xr.zeros_like(pp)

    for im, mask in enumerate(masks_plot):
        print(mask)
        for iy, eyear in enumerate(eyears):
            result = stats.pearsonr(
                XX_comb.sel(recon_type=rtype),
                YY_comb.sel(mask=mask, eyear=eyear, recon_type=rtype),
            )
            rr[im, iy] = result.statistic
            pp[im, iy] = result.pvalue
            rr_upper[im, iy] = result.confidence_interval().high
            rr_lower[im, iy] = result.confidence_interval().low
            pp_mc = mont.MonteCarlo_function(
                XX_comb.sel(recon_type=rtype).data,
                YY_comb.sel(mask=mask, eyear=eyear, recon_type=rtype).data,
                rr[im, iy].data,
            )
            pp_alt[im, iy] = pp_mc
        iy += 1
        z_mean = np.arctanh(rr[im, :iy]).mean()
        mean_corr = np.tanh(z_mean)
        corr_ci = [
            np.tanh(z_mean - 1.96 * 1 / np.sqrt(3 * len(XX) - 9)),
            np.tanh(z_mean + 1.96 * 1 / np.sqrt(3 * len(XX) - 9)),
        ]
        # print(mean_corr,corr_ci)
        rr[im, iy] = mean_corr
        rr_upper[im, iy] = corr_ci[1]
        rr_lower[im, iy] = corr_ci[0]
        chi = -2 * np.sum(np.log(pp[im, :iy].data))
        pp[im, iy] = 1 - chi2.cdf(chi, 6)
        chi = -2 * np.sum(np.log(pp_alt[im, :iy].data))
        pp_alt[im, iy] = 1 - chi2.cdf(chi, 6)
        iy += 1
        result = stats.pearsonr(
            XX_comb.sel(recon_type=rtype),
            YY_comb.sel(mask=mask, recon_type=rtype).mean("eyear"),
        )
        rr[im, iy] = result.statistic
        pp[im, iy] = result.pvalue
        rr_upper[im, iy] = result.confidence_interval().high
        rr_lower[im, iy] = result.confidence_interval().low
        pp_mc = mont.MonteCarlo_function(
            XX_comb.sel(recon_type=rtype).data,
            YY_comb.sel(mask=mask, recon_type=rtype).mean("eyear").data,
            rr[im, iy].data,
        )
        pp_alt[im, iy] = pp_mc
    rr.name = "R"
    rr_upper.name = "R_up"
    rr_lower.name = "R_lo"
    pp.name = "p"
    pp_alt.name = "p_alt"
    all_stats = xr.merge([rr, rr_upper, rr_lower, pp, pp_alt])
    stats_comb.append(all_stats)
stats_comb = xr.concat(stats_comb, "recon_type").assign_coords({"recon_type": rtypes})
stats_comb.to_netcdf(f"{DATA_DIR}/recon_correlationstats.nc")
