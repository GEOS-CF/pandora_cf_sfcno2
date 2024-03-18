#!/usr/local/other/python/GEOSpyD/2019.03_py3.7/2019-04-22/bin/python
'''
Python script to read Pandora L2_rnvh3p1-8 surface NO2 observations and calculate
corresponding quantities derived from GEOS-CF. These quantities are:
- cf_no2_sfcmr: GEOS-CF surface mixing ratio of NO2, in ppbv
- cf_no2_sfcconc: GEOS-CF surface concentration of NO2, in mol/m3
- cf_no2_l1col: GEOS-CF partial column amount in layer 1 [moles per square meter]
- cf_no2_pbl: GEOS-CF boundary layer height, [meters]

The corresponding Pandora quantities are labeled pandora_*. Conversion between
surface mixing ratios and concentrations is done using environmental quantities 
obtained from GEOS-CF.

EXAMPLES:
wget https://data.pandonia-global-network.org/WashingtonDC/Pandora140s1/L2/Pandora140s1_WashingtonDC_L2_rnvh3p1-8.txt
python pandora_cf_sfcno2.py -i Pandora140s1_WashingtonDC_L2_rnvh3p1-8.txt

HISTORY: 
20240205 - christoph.a.keller@nasa.gov - initial version
'''
import numpy as np
import os
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import xarray as xr
import argparse
import json

VV_TO_MOLEC = 1000./(9.80665*28.9644)
MINDATE=dt.datetime(2020,1,1)

# read pandora 
def main(args):
    '''
    Read Pandora observations and calculate corresponding GEOS-CF quantities
    '''
    locations = json.load(open(args.locations))
    assert args.nsite<len(locations),"index out of range: {} vs {}".format(args.nsite,len(locations)) 
    iloc = locations[args.nsite]
    basename = iloc.get('pandora_url').split('/')[-1]
    ifile = "obs/"+basename
    if not os.path.isfile(ifile):
        print("obs file does not exist, skip: {}".format(ifile))
        return

    ofile = "merged_csv/"+basename.replace(".txt","+GEOSCF.csv")
    if os.path.isfile(ofile):
        print("file exists, don't do anything: {}".format(ofile))
        return

    pand,lat,lon = _read_pandora(ifile)
    # qc flags?
    ##pand = pand.loc[pand['qval']==10,].copy()
    # limit data to minimum date and three days from now (due to CF latency)
    maxdate = dt.datetime.today() - dt.timedelta(days=3) 
    maxdate = dt.datetime(maxdate.year,maxdate.month,maxdate.day)
    pand = pand.loc[pand['date']>=MINDATE,].copy()
    pand = pand.loc[pand['date']<=maxdate,].copy()
    pand = pand.loc[(pand['pandora_no2_l1hgt']>0.0)&(pand['pandora_no2_l1hgt']<15.),].copy()
 
    # create empty entries for CF fields 
    pand['pandora_no2_sfcmr'] = [np.nan for i in range(pand.shape[0])]
    pand['cf_no2_sfcmr'] = [np.nan for i in range(pand.shape[0])]
    pand['cf_no2_sfcconc'] = [np.nan for i in range(pand.shape[0])]
    pand['cf_no2_l1col'] = [np.nan for i in range(pand.shape[0])]
    pand['cf_no2_pbl'] = [np.nan for i in range(pand.shape[0])]
    
    # populate CF fields
    fnl={}; dsl={}
    for irow, row in enumerate(pand.itertuples()):
        # get cf values for this entry
        sfcmr,sfcconc,l1col,pbl,sfcmr_pandora,fnl,dsl = _match_cf(args,row,lat,lon,fnl,dsl)
        # add to pandas array
        pand['pandora_no2_sfcmr'].values[irow] = sfcmr_pandora
        pand['cf_no2_sfcmr'].values[irow] = sfcmr
        pand['cf_no2_sfcconc'].values[irow] = sfcconc
        pand['cf_no2_l1col'].values[irow] = l1col
        pand['cf_no2_pbl'].values[irow] = pbl

    # add latitude and longitude to file
    pand['lat'] = [lat for i in range(pand.shape[0])]
    pand['lon'] = [lon for i in range(pand.shape[0])]
    
    # write to file
    pand.to_csv(ofile,index=False,date_format="%Y-%m-%d %H:%M")
    print("data written to {}".format(ofile))
    
    # plot l1 columns
    if False:
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
        plt.plot(pand['date'].values,pand['l1col'].values*1.0e6,label="Pandora",color="gray",alpha=0.5)
        plt.plot(pand['date'].values,pand['cf_l1col'].values*1.0e6,label="GEOS-CF",color="red",alpha=0.5)
        plt.gcf().autofmt_xdate()
        plt.ylabel('NO2 column in PANDORA layer 1 [1.0e-6 mol/m3]')
        plt.ylim(0,200)
        plt.legend()
        plt.gcf().set_size_inches(10,5)
        plt.savefig('pandora_l1_cols_DC.png',dpi=100)
        print('figure written to pandora_l1_cols_DC.png')
        plt.close()

    # plot surface mixing ratio 
    if False:
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
        plt.plot(pand['date'].values,pand['sfcmr'].values,label="Pandora",color="gray",alpha=0.5)
        plt.plot(pand['date'].values,pand['cf_sfcmr'].values,label="GEOS-CF",color="red",alpha=0.5)
        plt.gcf().autofmt_xdate()
        plt.ylabel('NO2 surface mixing ratio [ppbv]')
        plt.ylim(0,40.)
        plt.legend()
        plt.gcf().set_size_inches(10,5)
        plt.savefig('pandora_sfcmr_DC.png',dpi=100)
        print('figure written to pandora_sfcmr_DC.png')
        plt.close()

    return    
   
 
def _match_cf(args,row,lat,lon,fnl,dsl):
    '''
    Read GEOS-CF fields and calculate matching quantities at provided location (lat,lon)
    '''
    # initialize values
    sfcmr=np.nan; sfcconc=np.nan; l1col=np.nan; pbl=np.nan; sfcmr_pandora=np.nan
    # read new file if needed. Round to hour to match up with GEOS-CF time stamps
    idate = row.date.round('H')
    templ = idate.strftime(args.cf_template)
    # chm collection
    ifilec = templ.replace('<col>','chm')
    readchm = True
    if 'chm' in fnl:
        if fnl['chm']==ifilec:
            dsc = dsl['chm']
            readchm = False
    if readchm:
        if not os.path.isfile(ifilec):
           print('Warning - file not found: {}'.format(ifilec))
           return sfcmr, sfcconc, l1col, pbl, sfcmr_pandora, fnl, dsl
        print('reading {}'.format(ifilec))
        dsc = xr.open_dataset(ifilec)
        fnl['chm']=ifilec
        dsl['chm']=dsc.copy()
    # met collection
    ifilem = templ.replace('<col>','met')
    readmet = True
    if 'met' in fnl:
        if fnl['met']==ifilem:
            dsm = dsl['met']
            readmet = False
    if readmet:
        print('reading {}'.format(ifilem))
        dsm = xr.open_dataset(ifilem)
        fnl['met']=ifilem
        dsl['met']=dsm.copy()
    # for pbl collection, use tavg collection. Hence no rounding of input date
    ifilep = row.date.strftime(args.pbl_template)
    readpbl = True
    if 'pbl' in fnl:
        if fnl['pbl']==ifilep:
            dsp = dsl['pbl']
            readpbl = False
    if readpbl:
        print('reading {}'.format(ifilep))
        dsp = xr.open_dataset(ifilep)
        fnl['pbl'] = ifilep
        dsl['pbl'] = dsp.copy()
    # get pbl
    pbl = dsp.sel(lon=lon,lat=lat,method='nearest')['ZPBL'].values[0] 
    # surface NO2 mixing ratios in ppb and mol/m3
    iconc = dsc.sel(lon=lon,lat=lat,method='nearest')['NO2'].values[0,-1]
    prs = dsm.sel(lon=lon,lat=lat,method='nearest')['PS'].values[0]
    temp = dsm.sel(lon=lon,lat=lat,method='nearest')['T'].values[0,-1]
    sfcmr = iconc * 1.0e9
    sfcconc = iconc * prs / ( 8.314 * temp )
    # partial column in moles/m2
    ihgt = row.pandora_no2_l1hgt * 1000.
    delp = dsm.sel(lon=lon,lat=lat,method='nearest')['DELP'].values[0,::-1]
    zl   = dsm.sel(lon=lon,lat=lat,method='nearest')['ZL'].values[0,::-1]
    q    = dsm.sel(lon=lon,lat=lat,method='nearest')['Q'].values[0,::-1]
    no2  = dsc.sel(lon=lon,lat=lat,method='nearest')['NO2'].values[0,::-1]
    inML = True
    iL = 0
    l1col = 0.0
    while inML:
        # layer top height
        toph = (zl[iL]+zl[iL+1])/2. 
        if (toph>ihgt):
            both = (zl[iL]+zl[iL-1])/2. if iL>0 else 0.0
            frac = (ihgt-both) / (toph-both)
        else:
            frac = 1.0
        l1col += no2[iL] * delp[iL] * (1.-q[iL]) * VV_TO_MOLEC * frac
        iL += 1
        if frac < 1.0: 
            inML = False
        if iL>=len(zl)+1:
            print('Warning: Pandora layer height greater than entire atmosphere: {}'.format(ihgt))
            inML: False
    # convert reported pandora surface concentration from mol/m3 to ppbv. Inverse of
    # the GEOS-CF conversion above. Also dry out.
    sfcmr_pandora = row.pandora_no2_sfcconc / prs * ( 8.314 * temp ) * 1.0e9 / (1.-q[0])
    return sfcmr, sfcconc, l1col, pbl, sfcmr_pandora, fnl, dsl


def _read_pandora(ifile):
    '''
    Read pandora observations (L2_rnvh3p1-8)
    '''
    print('Reading {}'.format(ifile))
    dat = pd.read_csv(ifile,sep=',',skiprows=93,encoding="utf-8")
    tmp = []
    for l in range(dat.shape[0]):
        iline = dat.iloc[l]
        vals = iline[0].split(" ")
        idate = dt.datetime.strptime(vals[0],"%Y%m%dT%H%M%S.%fz")
        qval = int(vals[52])
        sfcmr = float(vals[55])
        l1hgt = float(vals[67])
        l1col = float(vals[68])
        irow = pd.DataFrame({"date":[idate],"pandora_no2_qval":[qval],"pandora_no2_sfcconc":[sfcmr],"pandora_no2_l1hgt":[l1hgt],"pandora_no2_l1col":[l1col]})
        tmp.append(irow)
    alldat = pd.concat(tmp)
    # get latitude and longitude
    f = open(ifile,'rb')
    lat = None 
    lon = None 
    while lat is None or lon is None: 
        line = str(f.readline())
        if "Location latitude" in line:
            lat = float(str(line.split(':')[1]).replace(" ","").replace("\\n","").replace("'",""))
            print("location latitude: {}".format(lat))
        if "Location longitude" in line:
            lon = float(str(line.split(':')[1]).replace(" ","").replace("\\n","").replace("'",""))
            print("location longitude: {}".format(lon))
    return alldat,lat,lon


def parse_args():
    p = argparse.ArgumentParser(description='Undef certain variables')
    p.add_argument('-l', '--locations',type=str,help='locations',default="PANDORA_Locations.json")
    p.add_argument('-n', '--nsite',type=int,help='location index in json list',default=0)
    #p.add_argument('-i', '--ifile',type=str,help='input pandora file',default=None)
    p.add_argument('-c', '--cf_template',type=str,help='GEOS-CF file template',default="/discover/nobackup/projects/gmao/geos_cf/pub/GEOS-CF_NRT/ana/Y%Y/M%m/D%d/GEOS-CF.v01.rpl.<col>_inst_1hr_g1440x721_v72.%Y%m%d_%H00z.nc4")
    p.add_argument('-p', '--pbl_template',type=str,help='GEOS-CF pbl file template',default="/discover/nobackup/projects/gmao/geos_cf/pub/GEOS-CF_NRT/ana/Y%Y/M%m/D%d/GEOS-CF.v01.rpl.met_tavg_1hr_g1440x721_x1.%Y%m%d_%H30z.nc4")
    return p.parse_args()

 
if __name__ == '__main__':
    main(parse_args())
#
