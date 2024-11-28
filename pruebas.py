import numpy as np
import pandas as pd


surveys = ['PS1','SDSS','GALEX','2MASS','WISE','unWISE','LegacySurvey','VISTA','SPLUS','UKIDSS']
filters = [['g','r','i','z','y'],['u','g','r','i','z'],['FUV','NUV'] , ['J','H','Ks'], ['W1','W2','W3','W4'],['W1','W2','W3','W4'],['g','r','z'],
           ['Z','Y','J','H','Ks'],['u','F378','F395','F410','F430','g','F515','r','F660','i','F861','z' ],['Z','Y','J','H','K']]
m_zp = [25,22.5,[18.82,20.08],'header','header',22.5,22.5,'header','header','header']
pixel_scale = [0.25,0.396,1.5,1.0,1.375,2.75,0.262,0.339,0.55,0.4]
resolution = [0.,0.,1.,1.0,1.,2.,0.,0.,0.,0.]


df = pd.DataFrame({'survey':surveys,'filters':filters,'m_zp':m_zp,'pizxel_scale':pixel_scale,'resolution':resolution})

print('i' in ''.join(df.loc[df['survey'] == 'SDSS', 'filters'].astype(str)))
print(type(np.array(df['m_zp'][df['survey'] == 'GALEX'])[0][0]))
print('SDSS' in np.array(df['survey']))
#indice = df.loc[df['survey'] == 'SDSS', 'filters'].apply(lambda x: 'i' in x).idxmax()
#filtered = df.loc[df['survey'] == 'SDSS']
#indices = filtered.index[filtered['filters'].str.contains('i')]

VALUE = {'SDSS':['g','r','i'],'WISE':['W1','W2']}
print(VALUE['SDSS'])




from hostphot.cutouts import download_images

name = 'SN2004eo'
host_ra, host_dec = 308.2092, 9.92755  # coords of host galaxy of SN2004eo
survey = 'PS1'
download_images(name, host_ra, host_dec, survey=survey,filters='g')

#df.to_csv('prm_config.csv')