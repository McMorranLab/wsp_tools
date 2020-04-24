import scipy.constants as sc

__all__ = ['s','m','kg','A','K','mol','cd','F','J','C','W','eV','c',
            'mu_0','eps_0','h','hbar','h_eV','hbar_eV','G','e','R','alpha',
            'N_A','k_B','k_B_eV','sigma','Wien','Rydberg','m_e','m_p',
            'm_n','pi']

def setUnits(second=1, meter=1, kilogram=1, Amp=1, Kelvin=1, mole=1, candela=1):
    global s,m,kg,A,K,mol,cd
    s,m,kg,A,K,mol,cd = second, meter, kilogram, Amp, Kelvin, mole, candela
    setConsts(s, m, kg, A, K, mol, cd)

def setConsts(s, m, kg, A, K, mol, cd):
    global F,J,C,W,eV,c,mu_0,eps_0,h,hbar,h_eV,hbar_eV,G,e,R,alpha,N_A,k_B
    global k_B_eV,sigma,Wien,Rydberg,m_e,m_p,m_n,pi
    F = s**4 * A**2 / m**2 / kg
    J = kg * m**2 / s**2
    C = A * s
    W = kg * m**2 / s**3
    eV = sc.e * J

    c       = sc.c          * m / s
    mu_0    = sc.mu_0       * m * kg / s**2 /A**2
    eps_0   = sc.epsilon_0  * F / m
    h       = sc.h          * J * s
    hbar    = sc.hbar       * J * s
    h_eV    = sc.h          * s * J / eV
    hbar_eV = sc.hbar       * s * J / eV
    G       = sc.G          * m**3 / kg / s**2
    g       = sc.g          * m / s**2
    e       = sc.e          * C
    R       = sc.R          * kg / s**2 / K / mol
    alpha   = sc.alpha
    N_A     = sc.N_A        / mol
    k_B     = sc.k          * J / K
    k_B_eV  = sc.k          * J / eV / K
    sigma   = sc.sigma      * W / m**2 / K**4
    Wien    = sc.Wien       * m * K
    Rydberg = sc.Rydberg    / m
    m_e     = sc.m_e        * kg
    m_p     = sc.m_p        * kg
    m_n     = sc.m_n        * kg
    pi      = sc.pi

setUnits()
