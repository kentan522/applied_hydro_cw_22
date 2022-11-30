from math import cos, exp, tanh, sinh, cosh, pi, sin
import cmath
from numpy import linspace
from scipy import special

class WaveCase():

    def __init__(self, h_max, T, D=40, d=100):
        self.h_max = h_max
        self.T = T
        
        # Fixed parameters
        self.a = self.h_max / 2 # Wave amplitude [m]
        self.D = D # Structure diameter [m]
        self.r = self.D / 2 # Structure radius [m]
        self.d = d # Water depth [m]
        self.rho = 1000 # Density of water [kg m^-3]
        self.g = 9.81 # Gravitational acceleration [ms^-2]

        # Calculated parameters
        self.omega = 2*pi / (self.T) # Angular velocity [rad/s]
        self.k = self.__calc_wave_num()
        self.umax = self.__calc_umax()
        self.kc = self.__calc_kc()
        self.wavelength = 2*pi / self.k
        self.d_over_lambda = self.D / self.wavelength

        # Forces and regimes
        self.force = self.__force_significance()
        self.regime = self.__regime()

        # Design calculations
        self.max_base_shear = self.__calc_base_shear()
        self.max_overturning_moment = self.__calc_overturning_moment()
        self.max_wave_surface = self.__calc_wave_surface_elevation()
        self.wave_run_up = self.__calc_wave_run_up()
        self.ka_val = self.__check_wave_breaking()
    
    def __calc_wave_num(self):

        k_old = (self.omega ** 2) / self.g
        diff = 1 # Initial difference

        while diff > 1e-5:
            k_new = (self.omega ** 2) / self.g*tanh(k_old * self.d)
            diff = abs(k_new - k_old)
            k_old = k_new
        
        return k_new
    
    def __calc_umax(self, z=0):

        nume = self.a * self.omega * cosh(self.k * (z + self.d))
        deno = sinh(self.k * self.d)

        return nume / deno

    def __calc_kc(self):

        return (self.umax*self.T) / self.D

    def __force_significance(self):

        if self.kc < 5:
            forces = 'Inertia Forces'
        elif self.kc > 20:
            forces = 'Drag Forces'
        else:
            forces = 'Inertia and Drag Forces'
        
        print(forces)
        return forces

    def __regime(self):
        
        if self.d_over_lambda >= 0.2:
            regime = 'Diffraction regime'
        else:
            regime = 'Long wave regime'
        
        print(regime)
        return regime
    
    def __calc_base_shear(self):
        
        # Long wave regime
        if self.regime == 'Long wave regime':
            t = linspace(0, self.T, 1000)
            F_tot = [0] * len(t)

            for i, time in enumerate(t):
                F1_1 = (2 * self.rho * self.g * self.a) *\
                     ((pi * self.D**2) / 4) *\
                     (cos(self.omega * time))
                F2_1 = 5/4 * (self.k * self.a * self.rho * self.g * self.a) *\
                     ((pi * self.D**2) / 4) *\
                     (sin(2 * self.omega * time))
                F3_1 = -1 * (self.k**2 * self.a**2 * self.rho * self.g * self.a) *\
                     ((pi * self.D**2) / 4) *\
                     (cos(3 * self.omega * time))
                F3_2 = (self.k**2 * self.a**2 * self.rho * self.g * self.a) *\
                     ((pi * self.D**2) / 4) *\
                     (cos(self.omega * time) - cos(3 * self.omega * time))

                F_tot[i] = F1_1 + F2_1 + F3_1 + F3_2
            
            F_max_shear = max(F_tot) / 1e6 # [MN]
        
        # Diffraction regime
        elif self.regime == 'Diffraction regime':
            Cm_star = 1.7 # Approx for D/lambda = 0.256 (Case 3)
            F_max_shear = Cm_star * ( (self.rho * pi**2 * self.D**2 * self.h_max * self.wavelength) / (4 * self.T ** 2) )
            F_max_shear /= 1e6

        print('Max base shear ', F_max_shear)

        return F_max_shear

    def __calc_overturning_moment(self):

        # Long wave regime
        if self.regime == 'Long wave regime':

            t = linspace(0, self.T, 1000)
            M_tot = [0] * len(t)

            for i, time in enumerate(t):
                M1 = 0.5 * pi * self.rho * self.g * self.k * self.a * (self.D**2) *\
                    cos(self.omega * time) *\
                    ((self.d / self.k) - (1 / (self.k ** 2)) *\
                    (1 - exp(-self.k * self.d)))
                M2 = 0.125 * pi * self.rho * self.g * (self.k**2) * (self.a**2) * (self.D**2) *\
                    sin(2 * self.omega * time) *\
                    ((self.d / (2 * self.k)) - (1 / (4 * (self.k ** 2))) + ((1 / (4 * (self.k ** 2))) *\
                    exp(-2 * self.k * self.d)))
                M3 = (self.d + (self.a * sin(self.omega * time))) *\
                    (0.25 * (self.k ** 2) * (self.a ** 3) * self.rho * self.g * pi * (self.D **2) *\
                    cos(self.omega * time) -\
                    0.5 * (self.k ** 2) * (self.a ** 3) * self.rho * self.g * pi * (self.D **2) *\
                    cos(3 * self.omega * time))
                
                M_tot[i] = M1 + M2 + M3
            
            M_max = max(M_tot) / 1e6 # [MNm]
        
        # Diffraction regime
        elif self.regime == 'Diffraction regime':
             Cm_star = 1.7 # Approx for D/lambda = 0.256 (Case 3)
             M_max = Cm_star * self.rho * self.g *\
                    ((self.h_max * self.wavelength * (self.D ** 2)) / 16) *\
                    (self.k * self.d * tanh(self.k * self.d) +\
                    (1 / cosh(self.k * self.d) - 1))
             M_max /= 1e6
        
        print('Max overturning moment ', M_max)

        return M_max
    
    def __bessel(self, m, ka):
        return special.jv(m, ka)
    
    def __bessel_prime(self, m, ka):
        return special.jv(m - 1, ka) - ((m / ka) * special.jv(m, ka))

    def __hankell(self, m, ka):
        return special.hankel1(m, ka)
    
    def __hankell_prime(self, m, ka):
        return ((m / ka) * special.hankel1(m, ka) - special.hankel1(m + 1, ka))

    def __calc_wave_surface_elevation(self):

        # Long wave regime
        if self.regime == 'Long wave regime':
            eta_max = self.a
        
        elif self.regime == 'Diffraction regime':
            theta = linspace(0, 2*pi, 100)
            t = linspace(0, self.T, 100)
            eta = len(theta) * [[0] * len(t)]

            for i, deg in enumerate(theta):
                for j, time in enumerate(t):
                    sum_m = 0
                    ka = (self.k * self.r)

                    for m in range(101):
                        if not m:
                            eps_m = 1
                        else:
                            eps_m = 2

                        sum_m += eps_m * (complex(0, 1)**(m + 1)) *\
                            (self.__bessel(m, ka) - (self.__bessel_prime(m, ka) / self.__hankell_prime(m, ka))*\
                            self.__hankell(m, ka)) *\
                            cos(m * deg) * cmath.exp(-1 * self.omega * time)

                    comp_eta = sum_m * self.h_max * 0.5

                    # Convert to real values
                    eta[i][j] = (comp_eta.real**2 + comp_eta.imag**2) ** 0.5
            
            eta_max = max(max(eta))

        print('Max water surface elevation ', eta_max)
        return eta_max

    def __calc_wave_run_up(self):

        if self.regime == 'Long wave regime':

            z = linspace(-self.d, self.a, 100)
            t = linspace(0, self.T, 100)

            u = len(t) * [[0] * len(z)]
            w = len(t) * [[0] * len(z)]
            u2w2 = len(t) * [[0] * len(z)]

            # Calculation of u
            for i, time in enumerate(t):
                for j, depth in enumerate(z):

                    u[i][j] = (self.a * self.omega * cosh(self.k * (depth + self.d))) * sin(self.omega * time)/\
                            sinh(self.k * self.d)
                    w[i][j] = (self.a * self.omega * sinh(self.k * (depth + self.d))) * cos(self.omega * time)/\
                            sinh(self.k * self.d)
                    u2w2[i][j] = u[i][j] ** 2 + w[i][j] ** 2

            R = max(max(u2w2))/(2*self.g) + self.max_wave_surface
        
        elif self.regime == 'Diffraction regime':

            R = (self.h_max / 2) * ((1 + 4 * ((self.r * self.k)**2))**0.5)

        print('Wave run-up on the face of the column ', R)
        return R

    def __check_wave_breaking(self):

        ka_val = self.k * self.a

        print('Value of ka is ', ka_val)
        if ka_val > 0.443: 
            print('Wave breaking occurs')
        else:
            print('Wave breaking does not occur')
        
        return
        
        






            

            











