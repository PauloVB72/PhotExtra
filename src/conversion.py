#code_working
from astropy import coordinates, units as u, wcs
from astropy.units import cds
from astropy.table.table import Table
import pandas as pd
import numpy as np


class Units:

    def __init__(self,unit = 'mJy'):

        self.units = unit

    def eval_unit(self):

        if self.units is 'mJy':
                cte_value = 1.
                print("units output: mJy")

        elif isinstance(self.units, str):
            units = {'F_w' : u.erg/(u.cm**2 * u.s *u.Angstrom) ,'F_v' : u.erg/(u.cm**2 * u.s *u.Hertz),
                            'f_w' : u.photon/(u.cm**2 * u.s *u.Angstrom)}
            try:
                assert self.units == 'F_w' or self.units == 'F_v' or self.units == 'f_w'

                if self.units == 'F_w':
                        cte_Fw = 3. * 10 ** -5 * 10 **-3
                        cte_value = cte_Fw
                        print("units output: " + str(units['F_w']))
                elif self.units == 'F_v':
                        cte_fw = 1.51 * 10 ** 3 * 10 **-3
                        cte_value = cte_fw
                        print("units output: " + str(units['F_v']))
                elif self.units == 'f_w':
                        cte_Fv = 10**-23 * 10 **-3
                        cte_value = cte_Fv
                        print("units output: " + str(units['f_w']))

            except AssertionError:
                    print("Error: incorrect name, input required a valid name.")
                    print("Avalibles: 'F_v: ergs/cm^2 s Hz', 'F_w: ergs/cm^2 s A ', 'f_w: photon/cm^2 s A '")

        return cte_value


class Surveys:

    def __init__(self,text,error = True):

        self.text = text
        self.name_srv = text.columns
        self.error = error

    def SDSS_u(self):

            fsdss_u = self.text['SDSS_u']
            magsdss_u = 22.5 - 2.5 * np.log10(fsdss_u) - 0.04
            fluxsdss_u = 10 ** ((22.5-magsdss_u)/2.5)
            sdss_u = fluxsdss_u * 3.631 * 10 ** -6 * 1000
            self.tbsdss_u = pd.DataFrame(sdss_u, columns= ['SDSS_u'])

            return self.tbsdss_u

    def SDSS_u_err(self):

            fsdss_u_e, x = self.text['SDSS_u_err'], self.tbsdss_u['SDSS_u']
            sdss_u_e = fsdss_u_e * 3.631 * 10 ** -6 * 1000 + x * 0.02
            self.tbsdss_u_e = pd.DataFrame(sdss_u_e, columns= ['SDSS_u_err'])

            return self.tbsdss_u_e

    def SDSS_g(self):

            fsdss_g  = self.text['SDSS_g']
            magsdss_g = 22.5 - 2.5 * np.log10(fsdss_g)
            fluxsdss_g = 10 ** ((22.5-magsdss_g)/2.5)
            sdss_g = fluxsdss_g* 3.631 * 10 ** -6 * 1000
            self.tbsdss_g = pd.DataFrame(sdss_g, columns= ['SDSS_g'])

            return self.tbsdss_g

    def SDSS_g_err(self):

            fsdss_g_e, x =  self.text['SDSS_g_err'], self.tbsdss_g['SDSS_g']
            sdss_g_e = fsdss_g_e * 3.631 * 10 ** -6 * 1000 + x * 0.01
            self.tbsdss_g_e = pd.DataFrame(sdss_g_e, columns= ['SDSS_g_err'])
            return  self.tbsdss_g_e

    def SDSS_r(self):

            fsdss_r = self.text['SDSS_r']
            magsdss_r = 22.5 - 2.5 * np.log10(fsdss_r)
            fluxsdss_r = 10 ** ((22.5-magsdss_r)/2.5)
            sdss_r = fluxsdss_r * 3.631 * 10 ** -6 * 1000
            self.tbsdss_r = pd.DataFrame(sdss_r, columns= ['SDSS_r'])

            return self.tbsdss_r

    def SDSS_r_err(self):

            fsdss_r_e, x = self.text['SDSS_r_err'],self.tbsdss_r['SDSS_r']
            sdss_r_e = fsdss_r_e * 3.631 * 10 ** -6 * 1000 + x * 0.01
            self.tbsdss_r_e = pd.DataFrame(sdss_r_e, columns= ['SDSS_r_err'])
            return self.tbsdss_r_e

    def SDSS_i(self):

            fsdss_i = self.text['SDSS_i']
            magsdss_i = 22.5 - 2.5 * np.log10(fsdss_i)
            fluxsdss_i = 10 ** ((22.5-magsdss_i)/2.5)
            sdss_i = fluxsdss_i * 3.631 * 10 ** -6 * 1000
            self.tbsdss_i = pd.DataFrame(sdss_i, columns= ['SDSS_i'])

            return self.tbsdss_i

    def SDSS_i_err(self):

            fsdss_i_e, x = self.text['SDSS_i_err'],self.tbsdss_i['SDSS_i']
            sdss_i_e = fsdss_i_e * 3.631 * 10 ** -6 * 1000 + x * 0.01
            self.tbsdss_i_e = pd.DataFrame(sdss_i_e, columns= ['SDSS_i_err'])
            return self.tbsdss_i_e

    def SDSS_z(self):

            fsdss_z = self.text['SDSS_z']
            magsdss_z = 22.5 - 2.5*np.log10(fsdss_z) + 0.02
            fluxsdss_z = 10 ** ((22.5-magsdss_z)/2.5)
            sdss_z = fluxsdss_z * 3.631 * 10 ** -6 * 1000
            self.tbsdss_z = pd.DataFrame(sdss_z, columns= ['SDSS_z'])

            return self.tbsdss_z

    def SDSS_z_err(self):

          fsdss_z_e, x =  self.text['SDSS_z_err'],self.tbsdss_z['SDSS_z']
          sdss_z_e = fsdss_z_e * 3.631 * 10 ** -6 * 1000 + x * 0.01
          self.tbsdss_z_e = pd.DataFrame(sdss_z_e, columns= ['SDSS_z_err'])
          return self.tbsdss_z_e

    def GALEX_FUV(self):

            lambda_fuv = 1538.6
            self.lambda_fuv= lambda_fuv
            fconv_f = 1.4 * 10 ** -15
            self.fconv_f = fconv_f

            fgalex_f = self.text['GALEX_FUV']
            fluxgalex_f = fconv_f * fgalex_f
            galex_f = fluxgalex_f * 3.34 * 10 ** 4 * lambda_fuv ** 2 * 1000
            galex_f = galex_f - 0.001 * galex_f
            self.tbgal_f = pd.DataFrame(galex_f, columns= ['GALEX_FUV'])

            return  self.tbgal_f

    def GALEX_FUV_err(self):

            fgalex_f_e, x = self.text['GALEX_FUV_err'], self.tbgal_f['GALEX_FUV']
            fluxgalex_f_e = self.fconv_f * fgalex_f_e
            galex_f_e = fluxgalex_f_e * 3.34 * 10 ** 4 * self.lambda_fuv ** 2 * 1000
            galex_f_e = galex_f_e + 0.01 * x
            self.tbgal_f_e = pd.DataFrame(galex_f_e, columns= ['GALEX_FUV_err'])
            return self.tbgal_f_e

    def GALEX_NUV(self):

            lambda_nuv = 2315.7
            self.lambda_nuv = lambda_nuv
            fconv_n = 2.06*10**-16 # conversion factor cps to erg/s/cm2/A
            self.fconv_n = fconv_n

            fgalex_n = self.text['GALEX_NUV']
            fluxgalex_n = fconv_n * fgalex_n
            galex_n = fluxgalex_n * 3.34 * 10 ** 4 * lambda_nuv ** 2 * 1000
            galex_n = galex_n -0.01 * galex_n
            self.tbgal_n = pd.DataFrame(galex_n, columns= ['GALEX_NUV'])
            return self.tbgal_n

    def GALEX_NUV_err(self):

            fgalex_n_e, x = self.text['GALEX_NUV_err'],self.tbgal_n['GALEX_NUV']
            fluxgalex_n_e = self.fconv_n * fgalex_n_e
            galex_n_e = fluxgalex_n_e * 3.34 * 10 ** 4 * self.lambda_nuv ** 2 * 1000
            galex_n_e = galex_n_e + 0.01 * x
            self.tbgal_n_e = pd.DataFrame(galex_n_e, columns= ['GALEX_NUV_err'])
            return self.tbgal_n_e

    def TWOMASS_J(self):

            MAGZP_J = 20.8749
            self.MAGZP_J = MAGZP_J
            FnuZP_J = 1594
            self.FnuZP_J = FnuZP_J
            lambda_J = 1.235
            self.lambda_J = lambda_J

            f2mass_j = self.text['TWOMASS_J']
            mag_j = MAGZP_J - 2.5*np.log10(f2mass_j)
            twomass_j = 10 ** (-mag_j/2.5) * FnuZP_J * 1000
            self.tbtwom_j = pd.DataFrame(twomass_j, columns= ['TWOMASS_J'])

            return self.tbtwom_j

    def TWOMASS_J_err(self):

            f2mass_j_e, x = self.text['TWOMASS_J_err'], self.tbtwom_j['TWOMASS_J']
            mag_j_e = self.MAGZP_J - 2.5*np.log10(f2mass_j_e)
            twomass_j_e = 10 ** (-mag_j_e/2.5) * self.FnuZP_J * 1000
            twomass_j_e = twomass_j_e + x * 0.01
            self.tbtwom_j_e = pd.DataFrame(twomass_j_e, columns= ['TWOMASS_J_err'])
            return self.tbtwom_j_e

    def TWOMASS_H(self):

            MAGZP_H = 20.4031
            self.MAGZP_H = MAGZP_H
            FnuZP_H = 1024
            self.FnuZP_H = FnuZP_H
            lambda_H = 1.662
            self.lambda_H = lambda_H

            f2mass_h = self.text['TWOMASS_H']
            mag_h = MAGZP_H - 2.5*np.log10(f2mass_h)
            twomass_h = 10 ** (-mag_h/2.5) * FnuZP_H * 1000
            self.tbtwom_h = pd.DataFrame(twomass_h, columns= ['TWOMASS_H'])

            return self.tbtwom_h

    def TWOMASS_H_err(self):

            f2mass_h_e, x = self.text['TWOMASS_H_err'],self.tbtwom_h['TWOMASS_H']
            mag_h_e = self.MAGZP_H - 2.5*np.log10(f2mass_h_e)
            twomass_h_e = 10 ** (-mag_h_e/2.5) * self.FnuZP_H * 1000
            twomass_h_e = twomass_h_e + x * 0.01
            self.tbtwom_h_e = pd.DataFrame(twomass_h_e, columns= ['TWOMASS_H_err'])
            return self.tbtwom_h_e

    def TWOMASS_K(self):

            MAGZP_K = 19.8759
            self.MAGZP_K = MAGZP_K
            FnuZP_K = 666.8
            self.FnuZP_K = FnuZP_K
            lambda_K = 2.159
            self.lambda_K = lambda_K

            f2mass_k = self.text['TWOMASS_K']
            mag_k = MAGZP_K - 2.5*np.log10(f2mass_k)
            twomass_k = 10 ** (-mag_k/2.5) * FnuZP_K * 1000
            self.tbtwom_k = pd.DataFrame(twomass_k, columns= ['TWOMASS_K'])

            return  self.tbtwom_k

    def TWOMASS_K_err(self):

            f2mass_k_e, x  = self.text['TWOMASS_K_err'],self.tbtwom_k['TWOMASS_K']
            mag_k_e = self.MAGZP_K - 2.5*np.log10(f2mass_k_e)
            twomass_k_e = 10 ** (-mag_k_e/2.5) * self.FnuZP_K * 1000
            twomass_k_e = twomass_k_e + x * 0.01
            self.tbtwom_k_e = pd.DataFrame(twomass_k_e, columns= ['TWOMASS_K_err'])

            return self.tbtwom_k_e

    def WISE_W1(self):

            M0_inst_W1 = 20.5
            self.M0_inst_W1 = M0_inst_W1
            AC_W1 = 0.222
            self.AC_W1 = AC_W1
            F0_W1 = 309.540
            self.F0_W1 = F0_W1
            Dm_W1 = 2.699
            self.Dm_W1 = Dm_W1

            fw1 = self.text['WISE_W1']
            Mcal_W1 = M0_inst_W1 - 2.5 * np.log10(fw1) - AC_W1
            F_W1 = F0_W1 * 10 ** (-1 * Mcal_W1/2.5)
            w1 = F_W1 * 10 ** (-1 * Dm_W1/2.5) * 1000
            self.tbw1 = pd.DataFrame(w1, columns= ['WISE_W1'])

            return self.tbw1

    def WISE_W1_err(self):

            fw1_e, x = self.text['WISE_W1_err'],self.tbw1['WISE_W1']
            Mcal_W1_e = self.M0_inst_W1 - 2.5 * np.log10(fw1_e) - self.AC_W1
            F_W1_e = self.F0_W1 * 10 ** (-1 * Mcal_W1_e/2.5)
            w1_e = F_W1_e * 10 ** (-1 * self.Dm_W1/2.5) * 1000
            w1_e = w1_e + x * 0.024
            self.tbw1_e = pd.DataFrame(w1_e, columns= ['WISE_W1_err'])
            return self.tbw1_e

    def WISE_W2(self):

            M0_inst_W2 = 19.5
            self.M0_inst_W2 = M0_inst_W2
            AC_W2 = 0.280
            self.AC_W2 = AC_W2
            F0_W2 = 171.787
            self.F0_W2 = F0_W2
            Dm_W2 = 3.339
            self.Dm_W2 = Dm_W2

            fw2 = self.text['WISE_W2']
            Mcal_W2 = M0_inst_W2 - 2.5 * np.log10(fw2) - AC_W2
            F_W2 = F0_W2 * 10 ** (-1 * Mcal_W2/2.5) * 1000
            w2 = F_W2 * 10 ** (-1 * Dm_W2/2.5)
            self.tbw2 = pd.DataFrame(w2, columns= ['WISE_W2'])

            return self.tbw2

    def WISE_W2_err(self):

            fw2_e, x = self.text['WISE_W2_err'],self.tbw2['WISE_W2']
            Mcal_W2_e = self.M0_inst_W2 - 2.5 * np.log10(fw2_e) - self.AC_W2
            F_W2_e = self.F0_W2 * 10 ** (-1 * Mcal_W2_e/2.5)
            w2_e = F_W2_e * 10 ** (-1 * self.Dm_W2/2.5) * 1000
            w2_e = w2_e + 0.028 * x
            self.tbw2_e = pd.DataFrame(w2_e, columns= ['WISE_W2_err'])

            return self.tbw2_e

    def WISE_W3(self):

            M0_inst_W3 = 18.0
            self.M0_inst_W3 = M0_inst_W3
            AC_W3 =  0.665
            self.AC_W3 = AC_W3
            F0_W3 = 31.674
            self.F0_W3 = F0_W3
            Dm_W3 = 5.174
            self.Dm_W3 = Dm_W3

            fw3 = self.text['WISE_W3']
            Mcal_W3 = M0_inst_W3 - 2.5 * np.log10(fw3) - AC_W3
            F_W3 = F0_W3 * 10 ** (-1 * Mcal_W3/2.5)
            w3 = F_W3 * 10 ** (-1 * Dm_W3/2.5) * 1000
            self.tbw3 = pd.DataFrame(w3, columns= ['WISE_W3'])

            return self.tbw3

    def WISE_W3_err(self):

            fw3_e, x = self.text['WISE_W3_err'], self.tbw3['WISE_W3']
            Mcal_W3_e = self.M0_inst_W3 - 2.5 * np.log10(fw3_e) - self.AC_W3
            F_W3_e = self.F0_W3 * 10 ** (-1 * Mcal_W3_e/2.5)
            w3_e = F_W3_e * 10 ** (-1 * self.Dm_W3/2.5) * 1000
            w3_e = w3_e + 0.045 * x
            self.tbw3_e = pd.DataFrame(w3_e, columns= ['WISE_W3_err'])

            return self.tbw3_e

    def selected_data(self):

        self.Galaxy = pd.DataFrame(np.asarray(self.text['object']), columns= ['Object'])
        TABLE_0 = list()
        TABLE_0.append(self.Galaxy)

        list_filter = {'SDSS_u' : self.SDSS_u,'SDSS_u_err' : self.SDSS_u_err,'SDSS_g' : self.SDSS_g,'SDSS_g_err' : self.SDSS_g_err,
                       'SDSS_r' : self.SDSS_r,'SDSS_r_err' : self.SDSS_r_err,
                       'SDSS_i' : self.SDSS_i,'SDSS_i_err' : self.SDSS_i_err,'SDSS_z' : self.SDSS_z,'SDSS_z_err' : self.SDSS_z_err,
                       'GALEX_FUV' : self.GALEX_FUV,'GALEX_FUV_err' : self.GALEX_FUV_err,
                       'GALEX_NUV' : self.GALEX_NUV,'GALEX_NUV_err' : self.GALEX_NUV_err,'TWOMASS_J' : self.TWOMASS_J,'TWOMASS_J_err' : self.TWOMASS_J_err,
                       'TWOMASS_H' : self.TWOMASS_H,'TWOMASS_H_err' : self.TWOMASS_H_err,'TWOMASS_K' : self.TWOMASS_K,
                       'TWOMASS_K_err' : self.TWOMASS_K_err,'WISE_W1' : self.WISE_W1,'WISE_W1_err' : self.WISE_W1_err,
                       'WISE_W2' : self.WISE_W2,'WISE_W2_err' : self.WISE_W2_err,'WISE_W3' : self.WISE_W3,'WISE_W3_err' : self.WISE_W3_err}
        for filter in list_filter.keys():
            for txt_filter in self.name_srv:
                if filter == txt_filter:
                    TABLE_0.append(list_filter[filter]())

        TB = pd.concat(TABLE_0, axis = 1)
        return TB


class conversion(Surveys):

    def __init__(self,txt,unit = 'mJy'):

        super().__init__(txt)
        self.TB = super().selected_data()
        self.units = unit
        self.units_cte = Units(unit).eval_unit()

    def data_conv(self):

        list_filter = {'SDSS_u' : 3543,'SDSS_g' : 4770, 'SDSS_r' : 6231,'SDSS_i' : 7625, 'SDSS_z' : 9134,

                'GALEX_FUV' : 1538.6,'GALEX_NUV' : 2315.7,'TWOMASS_J' : 1235,'TWOMASS_H' : 1662,'TWOMASS_K' : 2159,

                 'WISE_W1' : 3368,'WISE_W2' : 4618,'WISE_W3' : 12082,'UNWISE_W1': 3368,'UNWISE_W2' : 4618,

                 'UNWISE_W3' :  12082}

        for col in self.TB.columns:

            for col2  in list_filter.keys():

                if col == col2:
                    if self.units == 'mJy':
                        self.TB[col] =  self.TB[col] *1
                        if col+'_err' in self.TB.columns:
                            self.TB[col+'_err'] = self.TB[col+'_err'] *1

                    elif self.units == 'F_w':
                        self.TB[col] =  self.units_cte *self.TB[col] / list_filter[col] ** 2
                        if col+'_err' in self.TB.columns:
                            self.TB[col+'_err'] =  self.units_cte *self.TB[col+'_err'] / list_filter[col] ** 2

                    elif self.units == 'f_w':
                        self.TB[col] = self.units_cte * self.TB[col]  / list_filter[col]
                        if col+'_err' in self.TB.columns:
                            self.TB[col+'_err'] = self.units_cte * self.TB[col+'_err']  / list_filter[col]

                    elif self.units == 'F_v':
                        self.TB[col] = self.units_cte * self.TB[col]
                        if col+'_err' in self.TB.columns:
                            self.TB[col+'_err'] = self.units_cte * self.TB[col+'_err']
        return self.TB

class PYFlux(conversion):

    def __init__(self,txt,unit):

        super().__init__(txt,unit)
        super().data_conv()
        self.data = super().data_conv()
        print(self.data)

    def save_tb(self,path):

        self.data.to_csv(path,index=False)
        return "save file"

def randomgalname(number, nmb_name):

  count = 0
  lst = list()

  while count < number:

    count = count + 1
    c = string.ascii_letters
    val = [''.join(random.choice(c) for x in range(nmb_name))]
    kk = np.random.rand(1,26)[0].tolist()

    val.extend(kk)
    lst.append(val)

  df = pd.DataFrame(lst,columns=['object','SDSS_u','SDSS_u_err','SDSS_g','SDSS_g_err','SDSS_r','SDSS_r_err','SDSS_i','SDSS_i_err',
                                 'SDSS_z', 'SDSS_z_err','GALEX_FUV','GALEX_FUV_err','GALEX_NUV','GALEX_NUV_err','TWOMASS_J','TWOMASS_J_err',
                                 'TWOMASS_H','TWOMASS_H_err','TWOMASS_K','TWOMASS_K_err','WISE_W1','WISE_W1_err','WISE_W2','WISE_W2_err','WISE_W3','WISE_W3_err'])

  return df


txt_new = randomgalname(10,13)

path = 'hola.txt'

PYFlux(txt_new,'F_w')

