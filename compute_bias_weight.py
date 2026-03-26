import os, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import interpolate
from scipy.stats import binned_statistic

ARCSEC_TO_RAD = 1 / 206265.0

def load_obsdata(file_name):
    """load visibility data output from `ms_to_npz.py` """
    data = np.load(file_name)
    freq = data["freq_obs"] if "freq_obs" in data else None
    return data["u_obs"], data["v_obs"], data["vis_obs"], data["wgt_obs"], freq


def log_gridding_1d(xmin, xmax, n_d_log):
    """Log gridding"""
    dx = (np.log10(xmax) - np.log10(xmin)) / n_d_log
    xedges = np.linspace(np.log10(xmin), np.log10(xmax), n_d_log + 1)
    centers = 0.5 * dx + np.linspace(np.log10(xmin), np.log10(xmax), n_d_log + 1)
    coord = np.sort(np.append(10 ** xedges, 0))
    rep = np.sort(np.append(10 ** centers[:n_d_log], 0.5 * xmin))
    return coord, rep

def data_binning_1d(x, y, weights, grid):
    """1D binning (real & imaginary)"""
    w_sum, _, _ = binned_statistic(x, weights, "sum", bins=grid)
    r_sum, _, _ = binned_statistic(x, y.real * weights, "sum", bins=grid)
    i_sum, _, _ = binned_statistic(x, y.imag * weights, "sum", bins=grid)
    x_sum, _, _ = binned_statistic(x, x * weights, "sum", bins=grid)

    mask = w_sum > 0
    real = r_sum[mask] / w_sum[mask]
    imag = i_sum[mask] / w_sum[mask]
    x_m  = x_sum[mask] / w_sum[mask]

    vis = real + 1j * imag
    noise = (1 / np.sqrt(w_sum[mask])) * (1 + 1j)
    d_data = np.append(real, imag)
    sigma  = 1.0 / (np.append(noise.real, noise.imag) ** 2)
    return x_m, vis, noise, d_data, sigma

def deproject_uv_to_q(u, v, vis, cosi, pa, dx, dy):
    """Deproject (u,v) to q using disk deometry"""
    q = np.sqrt(
        (cosi**2) * (-u * np.cos(pa) + v * np.sin(pa))**2
        + (-u * np.sin(pa) - v * np.cos(pa))**2
    )
    phase = 2 * np.pi * (dx * u + dy * v)
    cosph, sinph = np.cos(phase), np.sin(phase)
    real = vis.real * cosph - vis.imag * sinph
    imag = vis.real * sinph + vis.imag * cosph
    return q, real + 1j * imag

def configure_plot():
    """ plot configure """
    mpl.rcParams['figure.autolayout'] = False
    mpl.rcParams['font.family'] = ['Times New Roman']
    mpl.rcParams['mathtext.fontset'] = 'cm'
    mpl.rcParams['mathtext.default'] = 'it'
    mpl.rcParams['font.size'] = 18
    mpl.rcParams['axes.grid'] = True
    mpl.rcParams['grid.linestyle'] = ''
    mpl.rcParams['axes.edgecolor'] = 'black'
    mpl.rcParams['xtick.color'] = 'black'
    mpl.rcParams['ytick.color'] = 'black'
    mpl.rcParams['figure.figsize'] = (6,5)

def plot_q_real(q, real, real_model):
    """ plot 1 """
    plt.plot(q, real, label="data")
    plt.plot(q, real_model, label="model")
    plt.xscale("log")
    plt.ylabel(r"Real Visibility [Jy]")
    plt.xlabel(r"Spatial frequency q [$\lambda$]")
    plt.legend()
    plt.show()

def plot_sigma_factor(q, real, imag, target_name =""):
    """ plot 2 """
    plt.plot(q, real, label="real")
    plt.xscale("log")
    plt.plot(q, imag, label="imag")
    plt.xscale("log")
    plt.ylabel(r"$\sigma_{\rm visibility}/\sigma_{\rm weight}$")
    plt.xlabel(r"Spatial frequency q [$\lambda$]")
    plt.title("%s" % target_name)
    plt.legend()
    plt.show()


def compute_sigma_bias(visfile, dx_arcsec, dy_arcsec, cosi, pa, q_min=100, n_grid=1000):
    """ Estimate biases in weights computing
    
            Ratios of standard deviation of visibilities of (data - model) and data error (1/sqrt(weight)).

        The ratios are gridded by deprojected spatial frequency (q).
        Model is constructed by assuming axisymmetric brightess distribution. 
        

    """
    #load data
    u, v, vis, wgt, _ = load_obsdata(visfile)
    dx, dy = dx_arcsec * ARCSEC_TO_RAD , dy_arcsec * ARCSEC_TO_RAD 

    #deproject q and phase shift
    q, vis_mod = deproject_uv_to_q(u, v, vis, cosi, pa, dx, dy)
    q_max = np.max(q)

    # bin q and visibility
    grid, _ = log_gridding_1d(q_min, q_max, n_grid)
    q_bin, vis_bin, noise, d_data, sigma = data_binning_1d(q, vis_mod.real, wgt, grid)

    # Make rough model using interpolation  & binning, and compute residual bewteen model and data. 
    f_interp = interpolate.interp1d(q_bin, vis_bin, fill_value="extrapolate")
    vis_model = f_interp(q)
    vis_model_bin= f_interp(q_bin)
    sigma_factor_r = (vis_mod.real - vis_model) * np.sqrt(wgt) ### (d_data-d_model)/sigma for real part
    sigma_factor_i = (vis_mod.imag) * np.sqrt(wgt) ### for imaginary part

    # Bin residual, q, and weights
    sigma_factor_grid, _, _ = binned_statistic(q, sigma_factor_r, "std", bins=grid) ## compute std of (d_data-d_model)/sigma in q-grid for real part
    sigma_factor_grid_i, _, _ = binned_statistic(q, sigma_factor_i, "std", bins=grid) ## for imaginary part
    q_w_grid, _, _ = binned_statistic(q, q * wgt, "sum", bins=grid)
    w_grid, _, _ = binned_statistic(q, wgt, "sum", bins=grid)
    mask = w_grid > 0 
    q_sigma_grid =q_w_grid[mask] / w_grid[mask] ## 
    
    # Biases in weights `\sigma_{visibility}/\sigma_{weight}` 
    # (Ratios of standard deviation of visibilities of (data - model) and data error (1/sqrt(weight))0

    sigma_factor_grid = sigma_factor_grid[mask]
    sigma_factor_grid_i = sigma_factor_grid_i[mask]

    return q_bin, vis_bin,  vis_model_bin, q_sigma_grid, sigma_factor_grid, sigma_factor_grid_i



if __name__ == '__main__':
	plot = False

    # Your visibility file (output from ms_to_npz.py)
	visfile ="./averaged_npz/AS209_continuum_averaged_corrected.vis.npz" 

    # Put your geometry:
    # dx_arcsec, dy_arcsec = disk centers (RA, Dec) in arcsec
    # cosi = cosine of inclination
    # pa (in unit of radian)
	dx_arcsec, dy_arcsec, cosi, pa = 1.9 * 0.001,-2.5* 0.001, np.cos(34.97 * np.pi/180.0),85.76 * np.pi/180.0

    # Main
	q_bin, vis_bin, vis_model_bin, q_sigma_grid, sigma_factor_grid, sigma_factor_grid_i = compute_sigma_bias(visfile, dx_arcsec, dy_arcsec, cosi, pa, q_min=100, n_grid=1000)

	# Weight bias
	weight_factor = np.median(sigma_factor_grid)**2  ## only real part but may use imaginary part

	print(" Weights in data (measurement set) should be divided by: %.4f" % weight_factor)

    # Can also check biases
	#print(np.median(sigma_factor_grid), np.median( sigma_factor_grid_i))
	#print(np.median(sigma_factor_grid)**2, np.median( sigma_factor_grid_i)**2)



	# plot 
	if plot:
		configure_plot()
		plot_q_real(q_bin, vis_bin, vis_model_bin)
		plot_sigma_factor(q_sigma_grid, sigma_factor_grid, sigma_factor_grid_i)
