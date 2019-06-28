import numpy as np
import copy
from scipy import interpolate

def cl_interval(lnL, sigma=None, CL=0.68):

    """ Calculate a confidence level interval from a log-likelihood image
    Simple area under the curve with the image collapsed along each
    dimension

    Parameters:
      lnL: np.array
        log-Likelihood image
      CL: float, optional
      sigma: float, optional
        Use to calculate confindence interval

    Returns:
      best_idx, all_error: Lists
        [best] [-, +] indices for each dimension
    """
    # Confidence limits
    if sigma is None:
        c0 = (1. - CL)/2.
        c1 = 1.-c0
    # Image dimensions
    shape = lnL.shape
    ndim = len(shape)
    slc = [slice(None)]*ndim
    # Find max
    norm_L = np.exp(np.maximum(lnL - np.max(lnL),-15.))
    # Find best indices
    indices = np.where(lnL == np.max(lnL))
    best_idx = [bi[0] for bi in indices]

    # Error intervals
    all_error = []
    for kk in range(ndim):
        # Collapse on this dimension
        slc = copy.deepcopy(best_idx)
        slc[kk] = slice(None)
        Lslice = norm_L[slc].flatten()
        # Interpolate and go
        cumul_area = np.cumsum(Lslice)
        f_area = interpolate.interp1d(cumul_area/cumul_area[-1], np.arange(len(Lslice)))
        # Here we go
        idx0 = int(np.round(f_area(c0)))
        idx1 = int(np.round(f_area(c1)))
        all_error.append([idx0,idx1])

    # Return
    return best_idx, all_error


def cl_indices(lnL, cl, sigma=False):
    """ Find the indices of a log-Likelihood grid encompassing a
    given confidence interval

    Parameters:
      lnL: np.array
        log-Likelihood image
      sigma: bool, optional
        Return as sigma values [not implemented]

    Returns:
      indices: Tuple of np.where output
    """
    # Max
    mxL = np.max(lnL)

    # Noramlize and flatten
    norm_img = lnL-mxL
    flat_img = norm_img.flatten()

    # Sort
    srt = np.argsort(flat_img)
    norm_lnL = flat_img[srt]

    # Sum
    cumulsum = np.cumsum(np.exp(np.maximum(norm_lnL,-15.)))
    cumul = cumulsum/cumulsum[-1]

    # Interpolation (smoothing a bit)
    fsum = interpolate.interp1d(norm_lnL, cumul)

    # Finish
    indices = np.where(fsum(norm_img) > (1-cl))

    # Return
    return indices

def fig_Ob_vs_H0(outfile='fig_Ob_vs_H0.pdf'):
    """
    Ob vs H0 contour plot

    Args:
        outfile:

    Returns:

    """

    like_file = '../Analysis/Output/Ob_H0_grid_small.fits'
    lnlhood, Obval, H0val = read_lhood(like_file)

    # Max
    iObmx, iH0mx = np.unravel_index(lnlhood.argmax(), lnlhood.shape)
    Ob_ML2 = Obval[iObmx]
    H0_ML2 = H0val[iH0mx]
    print("Ob ML2 = {}".format(Ob_ML2))
    print("H0 ML2 = {}".format(H0_ML2))

    # Error
    #best_idx, all_error = cas_cluster.cl_interval(rlnlhood)
    #rerr = rval[all_error[0]]
    #gerr = gval[all_error[1]]

    indices = frb_cosmology.cl_indices(lnlhood, 0.68)
    Oberr = Obval[[np.min(indices[0]), np.max(indices[0])]]
    H0err = H0val[[np.min(indices[1]), np.max(indices[1])]]

    # Get the contour levels
    mxL = np.max(lnlhood)
    norm_img = lnlhood-mxL
    flat_img = norm_img.flatten()
    # Sort
    srt = np.argsort(flat_img)
    norm_lnL = flat_img[srt]
    # Sum
    cumulsum = np.cumsum(np.exp(np.maximum(norm_lnL,-15.)))
    cumul = cumulsum/cumulsum[-1]
    # Interpolation (smoothing a bit)
    fsum = interp1d(cumul, norm_lnL)
    contours = [0.99,0.95,0.68]
    clevels=[]
    clrs = ('r','g','k')
    for contour in contours:
        clevels.append(fsum(1.-contour))
    # Add the zero contour
    clevels.append(0.)

    # Plot
    plt.figure(figsize=(6, 5))
    gs = gridspec.GridSpec(1,1)
    ax = plt.subplot(gs[0])

    X, Y = np.meshgrid(H0val, Obval)
    CS = ax.contourf(X, Y, norm_img, clevels,
                       colors=clrs,
                       origin='lower',
                       extend='neither')
    #CS.cmap.set_under('white')
    #CS.cmap.set_over('white')

    ax.scatter(H0_ML2, Ob_ML2, color='white', marker='*')

    # Values
    '''
    ax.text(0.55, 0.80, r'$r_0 = {:0.2f}'.format(r0_ML2)+'_{-'+'{:0.2f}'.format(
        r0_ML2-rerr[0])+'}^{+'+'{:0.2f}'.format(rerr[1]-r0_ML2)+'} \; (h_{100}^{-1} \, $ Mpc)',
            transform=ax.transAxes, size='large', ha='left')
    ax.text(0.55, 0.70, r'$\gamma = {:0.2f}'.format(gamma_ML2)+'_{-'+'{:0.2f}'.format(
        gamma_ML2-gerr[0])+'}^{+'+'{:0.2f}'.format(gerr[1]-gamma_ML2)+'}$',
            transform=ax.transAxes, size='large', ha='left')
    '''

    # Label me
    ax.set_xlabel(r'$H_0 \; \rm (km/s/Mpc)$')
    ax.set_ylabel(r'$\Omega_b$')
    ax.set_ylim(0., 0.1)

    # End
    plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 700
    plt.savefig(outfile, **kwargs)
    plt.close()