"""
Dadi demographic models for WY and CO populations of E. gillettii.
"""
import numpy
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum

#one ancestral population splits at time T ago with no migration
def bottleneck_split(params, (n1,n2), pts):
    nuW, nuC, T = params
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nuW, nuC)

    model_sfs = Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return model_sfs


#one dimensional demographic inference of bottlneck size and timing
def bottleneck_1d(params, n1, pts):
    nuC, T = params
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = Integration.one_pop(phi, xx, T, nuC)

    model_sfs = Spectrum.from_phi(phi, n1, (xx,))
    return model_sfs


#one ancestral population splits at time T ago, then exchanges m12 and m21 proportions of migrants
def split_migration(params, (n1,n2), pts):
    nuW, nuC, T, m12, m21 = params
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nuW, nuC, m12=m12, m21=m21)

    model_sfs = Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return model_sfs    


#one ancestral population splits at time T ago, then exchanges m12 proportions of migrants
def split_unidirectional_migration(params, (n1,n2), pts):
    nuW, nuC, T, m12  = params
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nuW, nuC, m12=m12)

    model_sfs = Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return model_sfs
                            


