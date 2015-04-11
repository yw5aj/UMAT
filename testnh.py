# %% Definition section
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from constants import EPS

NPTS = 20 # number of points in plotting
NPOW = 17 # number of raised power + 1

def get_stress_numerical(F, params, eps_s=EPS, output='Cauchy'):
    def get_psi(C, params):
        J = np.sqrt(np.linalg.det(C))
        Cbar = J**(-2./3.) * C
        Ibar1 = np.trace(Cbar)
        if params['model'] == 'Neo-Hookean':
            psi = params['G'] * (Ibar1-3.) + 1./params['D']*(J-1.)**2
        elif params['model'] == 'Holzapfel':
            Ibar411 = params['a1'].dot(Cbar.dot(params['a1']))
            Ibar422 = params['a2'].dot(Cbar.dot(params['a2']))
            Ebar1 = params['Kappa']*(Ibar1-3)+(1-3*params['Kappa'])*(Ibar411-1)
            Ebar2 = params['Kappa']*(Ibar1-3)+(1-3*params['Kappa'])*(Ibar422-1)
            psi = params['C10']*(Ibar1-3.) + 1./params['D']*((J**2-1)\
                /2-np.log(J)) + params['K1']/2/params['K2']*(
                np.exp(params['K2']*(Ebar1+np.abs(Ebar1))**2/4) - 1 +
                np.exp(params['K2']*(Ebar2+np.abs(Ebar2))**2/4) - 1)
        return psi
    J = np.linalg.det(F)
    C = np.dot(F.T, F)
    S = np.empty((3, 3))
    psi = get_psi(C, params)
    for i in range(3):
        for j in range(3):
            e_i, e_j = np.eye(3)[[i, j]]
            psiptb = get_psi(C + eps_s*(np.tensordot(e_i, e_j, 0) +
                np.tensordot(e_j, e_i, 0)), params)
            S[i, j] = (psiptb - psi) / eps_s
    tau = F.dot(S).dot(F.T)
    sigma = tau / J
    if output == 'Cauchy':
        stress = sigma
    elif output == 'Kirchoff':
        stress = tau
    elif output == 'PK2':
        stress = S
    return stress


def get_stress_theoretical(F, params, output='Cauchy', **keywords):
    """
    Calculates Cauchy stress based on model type.
    """
    J = np.linalg.det(F)
    F_bar = J**(-1./3.) * F
    B_bar = F_bar.dot(F_bar.T)
    if params['model'] == 'Neo-Hookean':
        sigma = 2./J*params['G']*(B_bar - 1./3.*np.trace(B_bar)*np.eye(3)) + \
            2./params['D']*(J-1.)*np.eye(3)
    else:
        raise(NotImplementedError(params['model']))
    # Choose output type
    if output == 'Cauchy':
        stress = sigma
    elif output == 'Kirchoff':
        stress = J * sigma
    return stress


def get_C_CJ_ij(F, params, i, j, eps_c=EPS, eps_s=EPS,
                get_stress=get_stress_numerical):
    def perturb(F, i, j, eps_c=EPS):
        """
        To purturb deformation gradient on its (i, j) component
        """
        e_i, e_j = np.eye(3)[[i, j]]
        del_F = eps_c/2. * (np.outer(e_i, e_j).dot(F)
            + np.outer(e_j, e_i).dot(F))
        F_hat = F + del_F
        return F_hat
    J = np.linalg.det(F)
    tau = get_stress(F, params, output='Kirchoff', eps_s=eps_s)
    tau_perturb = get_stress(perturb(F, i, j, eps_c=eps_c), params,
                             output='Kirchoff', eps_s=eps_s)
    C_CJ_ij = 1. / J / eps_c * (tau_perturb - tau)
    return C_CJ_ij


def get_C_CJ_numerical(F, params, get_stress=get_stress_numerical,
                       eps_c=EPS, eps_s=EPS):
    C_CJ_numerical = np.empty([3, 3, 3, 3])
    for i in range(3):
        for j in range(3):
            C_CJ_numerical[:, :, i, j] = get_C_CJ_ij(F, params, i, j,
                eps_c=eps_c, eps_s=eps_s, get_stress=get_stress)
    return C_CJ_numerical


def ccj2ccc(ccj, sigma):
    delta = np.eye(3)
    ccc = ccj - 0.5 * (np.einsum('ik, jl', delta, sigma) +
                       np.einsum('jk, il', delta, sigma) +
                       np.einsum('il, jk', delta, sigma) +
                       np.einsum('jl, ik', delta, sigma))
    return ccc


def get_C_CJ_theoretical(F, params):
    J = np.linalg.det(F)
    F_bar = J**(-1./3.) * F
    B_bar = F_bar.dot(F_bar.T)
    C_CJ_theoretical = np.empty([3, 3, 3, 3])
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    if params['model'] == 'Neo-Hookean':
                        C_CJ_theoretical[i, j, k, l] = 2. / J * params['G'] *\
                            (0.5 * (np.eye(3)[i, k] * B_bar[j, l] +
                            B_bar[i, k] * np.eye(3)[j, l] +
                            np.eye(3)[i, l] * B_bar[j, k] +
                            B_bar[i, l] * np.eye(3)[j, k]) -
                            2./3. * (
                            np.eye(3)[i, j] * B_bar[k, l] +
                            B_bar[i, j] * np.eye(3)[k, l]) +
                            2./9. * np.eye(3)[i, j] * np.eye(3)[k, l] *
                            np.trace(B_bar)) + (2. / params['D'] * (2.*J - 1.)
                            * np.eye(3)[i, j] * np.eye(3)[k, l])
    return C_CJ_theoretical


def plot_error(dfgrd_list, maxindex, params, axs=None, add_color_map=False):
    # Check whether axs were passed in
    if axs is None:
        fig, axs = plt.subplots(3, 1, figsize=(3.27, 6.83))
    # Vary epsilon
    err_c = np.empty((NPOW, NPOW))
    err_s = np.empty(NPOW)
    for i in range(NPOW):
        eps_s = np.power(10., -1*i)
        # Select trace to plot
        do_plot = True if i in [6, 8] else False
        kwargs = {'ls': '-', 'color': 'k'} if i == 6 else {
            'ls': '-.', 'color':'k'}
        ssres, sstot = plot_stress(dfgrd_list, maxindex, params, eps_s, axs[0],
                    do_plot, **kwargs)
        err_s[i] = ssres / sstot
        for j in range(NPOW):
            eps_c = np.power(10., -1*j)
            ccj_theoretical = get_C_CJ_theoretical(dfgrd_list[-1], params)
            ccj_numerical  = get_C_CJ_numerical(dfgrd_list[-1], params,
                get_stress=get_stress_numerical, eps_c=eps_c, eps_s=eps_s)
            # Get ssres/sstot
            ccj_theoretical = np.nan_to_num(ccj_theoretical)
            ccj_numerical = np.nan_to_num(ccj_numerical)
            ssres = ((ccj_theoretical - ccj_numerical)**2).sum()
            sstot = ccj_theoretical.var() * ccj_theoretical.size
            err_c[i, j] = ssres / sstot
    # Plot the stress error
    axs[1].plot(range(NPOW), err_s, '-k')
    axs[1].set_yscale('log')
    xticks = axs[1].get_xticks()
    xticks = xticks[::2]
    axs[1].set_xticks(xticks)
    axs[1].set_xticklabels([r'$10^{-%d}$'%tick for tick in xticks])
    axs[1].set_xlabel(r'$\varepsilon_S$')
    axs[1].set_ylabel(r'$SS_{res}/SS_{tot}$ for $\sigma_{%d%d}$'%
        (maxindex[0]+1, maxindex[1]+1))
    # Plot the ccj contour plot
    from matplotlib.colors import LogNorm
    im = axs[2].imshow(err_c, norm=LogNorm(vmin=1e-9, vmax=1),
        interpolation='none')
    # Format the contour plot
    axs[2].set_xlabel(r'$\varepsilon_C$')
    axs[2].set_ylabel(r'$\varepsilon_S$')
    axs[2].set_xlim(1, NPOW-1)
    axs[2].set_ylim(1, NPOW-1)
    xticks = axs[2].get_xticks()
    yticks = axs[2].get_yticks()
    axs[2].set_xticklabels([r'$10^{-%d}$'%tick for tick in xticks],
        rotation='vertical')
    axs[2].set_yticklabels([r'$10^{-%d}$'%tick for tick in yticks])
    # Find the best solution
#    optimal = np.unravel_index(err_c.ravel().argmin(), err_c.shape)
#    axs[2].plot(optimal[1], optimal[0], '*', ms=10, mfc='k')
    # Add color map
    if add_color_map:
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        axin = inset_axes(axs[2], width='100%', height='100%',
                              bbox_to_anchor=(1.01, 0., .04, 1.),
                              bbox_transform=axs[2].transAxes, borderpad=0)
        plt.colorbar(im, cax=axin)
    return axs


def plot_stress(dfgrd_list, maxindex, params, eps_s, axs, do_plot=True,
                **kwargs):
    # Get the first plot
    stretch = []
    stress_theoretical = []
    stress_numerical = []
    for dfgrd in dfgrd_list:
        # Get stretch
        u = linalg.sqrtm(dfgrd.T.dot(dfgrd))
        stretch.append(u[maxindex])
        # Get stresses
        stress_theoretical.append(get_stress_theoretical(dfgrd, params)
            [maxindex])
        stress_numerical.append(get_stress_numerical(dfgrd, params, eps_s=eps_s
            )[maxindex])
    stress_theoretical = 1e-6*np.array(stress_theoretical)
    stress_numerical = 1e-6*np.array(stress_numerical)
    if do_plot:
        axs.plot(stretch, stress_theoretical, 'o', label='Analytic',
            color='k')
        axs.plot(stretch, stress_numerical, '-',
                 label=r'$\varepsilon_S=10^{%d}$'%np.log10(eps_s).astype(int),
                 **kwargs)
    axs.set_xlabel(r'$\lambda_{%d%d}$' % tuple(np.array(maxindex)+1))
    axs.set_ylabel(r'$\sigma_{%d%d}$ (MPa)' % tuple(np.array(maxindex)+1))
    ssres = ((stress_numerical - stress_theoretical)**2).sum()
    sstot = stress_theoretical.var()*stress_theoretical.shape[0]
    return ssres, sstot


# %% Main code
if __name__ == '__main__':
#    from constants import f as F
    # %% Get related quantities
    params = dict(G=80e3, D=2e-1, model='Neo-Hookean')
    # %% Generate the dfgrd_list for uniaxial compression-tension
    f11 = np.linspace(.25, 4., NPTS)
    f22 = np.sqrt(1./f11)
    f33 = f22
    dfgrd_list_uniaxial = []
    for i, dfgrd11 in enumerate(f11):
        dfgrd22 = f22[i]
        dfgrd33 = f33[i]
        dfgrd_list_uniaxial.append(np.diagflat([dfgrd11, dfgrd22, dfgrd33]))
    # Generate for biaxial stretch
    f11 = np.linspace(1., 4., NPTS)
    f22 = f11
    f33 = 1. / (f11 * f22)
    dfgrd_list_biaxial = []
    for i, dfgrd11 in enumerate(f11):
        dfgrd22 = f22[i]
        dfgrd33 = f33[i]
        dfgrd_list_biaxial.append(np.diagflat([dfgrd11, dfgrd22, dfgrd33]))
    # Generate for shear
    dfgrd_list_shear = [np.eye(3) for i in range(NPTS)]
    f12 = np.linspace(0., 1., NPTS)
    for i, dfgrd in enumerate(dfgrd_list_shear):
        dfgrd[0, 1] = f12[i]
    # Start plotting
    fig, axs = plt.subplots(3, 3, figsize=(6.83, 6.83))
    plot_error(dfgrd_list_uniaxial, (0, 0), params, axs[:, 0])
    plot_error(dfgrd_list_biaxial, (0, 0), params, axs[:, 1])
    plot_error(dfgrd_list_shear, (0, 1), params, axs[:, 2],
               add_color_map=True)
    axs[0, 0].set_ylim(-.5, 2.5)
    axs[0, 1].set_xlim(1, 4)
    axs[0, 1].set_ylim(-.0, 1)
    axs[0, 2].set_ylim(-.0, .2)
    axs[0, 0].set_xticks(np.arange(5))
    axs[0, 1].set_xticks(np.arange(1, 5))
    axs[0, 2].set_xticks(np.arange(0., .5, .1))
    # Add legend
    handles, labels = axs[0, 0].get_legend_handles_labels()
#    axs[0, 0].legend([handles[0]]+handles[1::2][::8],
#                     [labels[0]]+labels[1::2][::8], loc=2,
#                     fontsize=8, labelspacing=0.)
    axs[0, 0].legend([handles[0]]+handles[1::2][::1],
                     [labels[0]]+labels[1::2][::1], loc=2,
                     fontsize=8, labelspacing=0.)
    # Add titles
    axs[0, 0].set_title('Uniaxial compression/tension')
    axs[0, 1].set_title('Biaxial tension')
    axs[0, 2].set_title('Simple shear')
    axs[2, 0].set_title(r'$SS_{res}/SS_{tot}$ for $\mathbb{C}^{\sigma J}$')
    axs[2, 1].set_title(r'$SS_{res}/SS_{tot}$ for $\mathbb{C}^{\sigma J}$')
    axs[2, 2].set_title(r'$SS_{res}/SS_{tot}$ for $\mathbb{C}^{\sigma J}$')
    # Tighten figure and save
    fig.tight_layout()
    for axes_id, axes in enumerate(axs.ravel()):
        axes.text(-.3, 1.05, chr(65+axes_id), transform=axes.transAxes,
            fontsize=12, fontweight='bold', va='top')
    fig.subplots_adjust(right=.95)
    fig.savefig('./plots/singleelem.png')

