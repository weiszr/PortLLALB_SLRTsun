
    
def add_values(arr,mc_seal):
    from numpy import shape, zeros
    nrows, ncols = shape(arr)
    new_array = zeros([mc_seal, ncols])
    if mc_seal/nrows < 2:
        for i_col in range(ncols):
            new_array[nrows-1:mc_seal-1,i_col] = arr[0:(mc_seal-nrows)-1,i_col]
    else:
        for i_col in range(ncols):
            for i in range(mc_seal-nrows):
                index = i % nrows
                new_array[nrows-1+i,i_col] = arr[index,i_col]
    return new_array

def add_values_single(arr,mc_seal):
    from numpy import zeros
    nrows = len(arr)
    new_array = zeros(mc_seal)
    if mc_seal/nrows < 2:
        new_array[nrows-1:mc_seal-1] = arr[0:(mc_seal-nrows)-1]
    else:
        for i in range(mc_seal-nrows):
            index = i % nrows
            new_array[nrows-1+i] = arr[index]
    return new_array



def read_sea_level_pd(scen,years,mc_seal=1000):
    from numpy import shape,zeros,array,loadtxt,genfromtxt,where,in1d,insert
    import scipy
    import pandas as pd
    import itertools
    cols = []
    for i in range(len(years)):
        cols.append(str(int(years[i])))
    if scen == 'RCP85NA':
        dummy_array = loadtxt('loc_slr_mcsamps_losangeles_rcp85_kopp2014.txt',skiprows=4)
        dummy_nx,dummy_ny=shape(dummy_array)
        sl_data=zeros([dummy_nx,dummy_ny+1])
        sl_data[:,1::] = loadtxt('loc_slr_mcsamps_losangeles_rcp85_kopp2014.txt',skiprows=4)
        with open('loc_slr_mcsamps_losangeles_rcp85_kopp2014.txt') as f:
            year_file=f.readlines()[2].split()
        year_file = array(year_file,dtype=int)
        indices = where(in1d(year_file, years))[0]
        indices = insert(indices, 0, -1)
        indices[:] = indices[:]+1
        sl_data = sl_data * 0.001
        sl_data = array(sl_data)
        nmc,nyear = shape(sl_data)
        if nmc<mc_seal:
            sl_data = add_values(sl_data,mc_seal)
        newsl_data = zeros([mc_seal, len(indices)])
        for i in range(len(indices)):
            i_year = indices[i]
            dd=sl_data[:,i_year]
            if min(dd) == 0.0 and max(dd) == 0.0:
                newsl_data[:,i] = 0.0
            else:
               # print(len(dd))
                sample_pdf = scipy.stats.gaussian_kde(sl_data[:,i_year])
                dummy = sample_pdf.resample(mc_seal).T[:,0]
                newsl_data[:,i] = dummy
        df = pd.DataFrame(newsl_data, columns=cols)
        return df
    if scen == 'RCP85WA':   
        dummy_array = loadtxt('LA_slr_mc_subset_rcp85_dpais.txt',skiprows=1)
        dummy_nx,dummy_ny=shape(dummy_array)
        sl_data=zeros([dummy_nx,dummy_ny+1])
        sl_data[:,1::] = loadtxt('LA_slr_mc_subset_rcp85_dpais.txt',skiprows=1)
        with open('LA_slr_mc_subset_rcp85_dpais.txt') as f:
            year_file=f.readlines()[0].split()
        year_file = array(year_file,dtype=int)
        indices = where(in1d(year_file, years))[0]
        indices = insert(indices, 0, -1)
        indices[:] = indices[:]+1
        sl_data = sl_data * 0.001
        nmc,nyear = shape(sl_data)
        newsl_data = zeros([mc_seal, len(indices)])
        if nmc<mc_seal:
            sl_data = add_values(sl_data,mc_seal)
        for i in range(len(indices)):
            i_year = indices[i]
            dd=sl_data[:,i_year]
           # print(dd)
            if min(dd) == 0.0 and max(dd) == 0.0:
                newsl_data[:,i] = 0.0
            else:
                #print(len(dd))
                sample_pdf = scipy.stats.gaussian_kde(sl_data[:,i_year])
                dummy = sample_pdf.resample(mc_seal).T[:,0]
                newsl_data[:,i] = dummy
        df = pd.DataFrame(newsl_data, columns=cols)
        return df
    if scen == 'RCP45WA':    
        dummy_array = loadtxt('LA_slr_mc_subset_rcp45_dpais.txt',skiprows=1)
        dummy_nx,dummy_ny=shape(dummy_array)
        sl_data=zeros([dummy_nx,dummy_ny+1])
        sl_data[:,1::] = loadtxt('LA_slr_mc_subset_rcp45_dpais.txt',skiprows=1)
        with open('LA_slr_mc_subset_rcp45_dpais.txt') as f:
            year_file=f.readlines()[2].split()
        year_file = array(year_file,dtype=int)
        indices = where(in1d(year_file, years))[0]
        indices = insert(indices, 0, -1)
        indices[:] = indices[:]+1
        sl_data = sl_data * 0.001
        nmc,nyear = shape(sl_data)
        newsl_data = zeros([mc_seal, len(indices)])
        if nmc<mc_seal:
            sl_data = add_values(sl_data,mc_seal)
        for i in range(len(indices)):
            i_year = indices[i]
            dd=sl_data[:,i_year]
            #print(dd)
            if min(dd) == 0.0 and max(dd) == 0.0:
                newsl_data[:,i] = 0.0
            else:
              #  print(len(dd))
                sample_pdf = scipy.stats.gaussian_kde(dd[:].T)
                dummy = sample_pdf.resample(mc_seal).T[:,0]
                newsl_data[:,i] = dummy
        df = pd.DataFrame(newsl_data, columns=cols)
        return df
    if scen == 'RCP45NA':    
        dummy_array = loadtxt('LA_slr_mc_subset_rcp45_k2014.txt.txt',skiprows=1)
        dummy_nx,dummy_ny=shape(dummy_array)
        sl_data=zeros([dummy_nx,dummy_ny+1])
        sl_data[:,1::] = loadtxt('LA_slr_mc_subset_rcp45_k2014.txt.txt',skiprows=1)
        with open('LA_slr_mc_subset_rcp45_k2014.txt.txt') as f:
            year_file=f.readlines()[2].split()
        year_file = array(year_file,dtype=int)
        indices = where(in1d(year_file, years))[0]
        indices = insert(indices, 0, -1)
        indices[:] = indices[:]+1
        sl_data = sl_data * 0.001
        nmc,nyear = shape(sl_data)
        newsl_data = zeros([mc_seal, len(indices)])
        if nmc<mc_seal:
            sl_data = add_values(sl_data,mc_seal)
        for i in range(len(indices)):
            i_year = indices[i]
            dd=sl_data[:,i_year]
           # print(dd)
            if min(dd) == 0.0 and max(dd) == 0.0:
                newsl_data[:,i] = 0.0
            else:
                #print(len(dd))
                sample_pdf = scipy.stats.gaussian_kde(sl_data[:,i_year])
                dummy = sample_pdf.resample(mc_seal).T[:,0]
                newsl_data[:,i] = dummy
        df = pd.DataFrame(newsl_data, columns=cols)
        return df    
    if scen == 'RCP26NA':
        dummy_array = loadtxt('loc_slr_mcsamps_losangeles_rcp26_k2014.txt',skiprows=4)
        dummy_nx,dummy_ny=shape(dummy_array)
        sl_data=zeros([dummy_nx,dummy_ny+1])
        sl_data[:,1::] = loadtxt('loc_slr_mcsamps_losangeles_rcp26_k2014.txt',skiprows=4)
        with open('loc_slr_mcsamps_losangeles_rcp26_k2014.txt') as f:
            year_file=f.readlines()[2].split()
        year_file = array(year_file,dtype=int)
        indices = where(in1d(year_file, years))[0]
        indices = insert(indices, 0, -1)
        indices[:] = indices[:]+1
        sl_data = sl_data * 0.001
        nmc,nyear = shape(sl_data)
        newsl_data = zeros([mc_seal, len(indices)])
        if nmc<mc_seal:
            sl_data = add_values(sl_data,mc_seal)
        for i in range(len(indices)):
            i_year = indices[i]
            dd=sl_data[:,i_year]
           # print(dd)
            if min(dd) == 0.0 and max(dd) == 0.0:
                newsl_data[:,i] = 0.0
            else:
                #print(len(dd))
                sample_pdf = scipy.stats.gaussian_kde(sl_data[:,i_year])
                dummy = sample_pdf.resample(mc_seal).T[:,0]
                newsl_data[:,i] = dummy
        df = pd.DataFrame(newsl_data, columns=cols)
        return df

    if scen == 'RCP26WA':
        dummy_array = loadtxt('loc_slr_mcsamps_losangeles_rcp26_dpais.txt',skiprows=4)
        dummy_nx,dummy_ny=shape(dummy_array)
        sl_data=zeros([dummy_nx,dummy_ny+1])
        sl_data[:,1::] = loadtxt('loc_slr_mcsamps_losangeles_rcp26_dpais.txt',skiprows=4)
        with open('loc_slr_mcsamps_losangeles_rcp26_dpais.txt') as f:
            year_file=f.readlines()[2].split()
        year_file = array(year_file,dtype=int)
        indices = where(in1d(year_file, years))[0]
        indices = insert(indices, 0, -1)
        indices[:] = indices[:]+1
        sl_data = sl_data * 0.001
        nmc,nyear = shape(sl_data)
        newsl_data = zeros([mc_seal, len(indices)])
        if nmc<mc_seal:
            sl_data = add_values(sl_data,mc_seal)
        for i in range(len(indices)):
            i_year = indices[i]
            dd=sl_data[:,i_year]
           # print(dd)
            if min(dd) == 0.0 and max(dd) == 0.0:
                newsl_data[:,i] = 0.0
            else:
                #print(len(dd))
                sample_pdf = scipy.stats.gaussian_kde(sl_data[:,i_year])
                dummy = sample_pdf.resample(mc_seal).T[:,0]
                newsl_data[:,i] = dummy
        df = pd.DataFrame(newsl_data, columns=cols)
        return df

    
def adjust_tsunami_dicts(tsu_d,seal,year,n_sea):
    cols=seal.columns.values
    data_out = []
    new_col =[]
    n_se = len(seal)
    m=-1
    data_frame = {}
    se = seal[year].values
    col_names = array(list(tsu_d.columns))
    eq_dd = tsu_d[col_names[0]].values
    mc_eq=len(eq_dd)
    n_col= len(col_names)
    data = zeros([15, len(se)*mc_eq])
    m=-1
    for i_eq in trange(n_col,desc='Processing year {t1}'.format(t1=year),disable=mode):
        te = tsu_d[col_names[i_eq]].values
        dummy = []
        for i_eq1 in range(len(te)):
            for i_seal in range(len(se)):
                dummy.append((te[i_eq1]+se[i_seal]))
        dummy = array(dummy)
        data[i_eq,:] = dummy[:]
    data_frame= pd.DataFrame(transpose(data),columns=["8.0","8.1","8.2","8.3","8.4","8.5","8.6","8.7","8.8","8.9","9.0","9.1","9.2","9.3","9.4"]) 
    return data_frame

def adjust_tsunami_dicts1(tsu_d,seal,year,mode_l):
    from numpy import array,zeros,shape,transpose,mean
    from tqdm import trange
    import pandas as pd
    cols=seal.columns.values
    data_out = []
    new_col =[]
    n_se = len(seal)
    m=-1
    data_frame = {}
    se = seal[year].values
    col_names = array(list(tsu_d.columns))
    eq_dd = tsu_d[col_names[0]].values
    mc_eq=len(eq_dd)
    n_col= len(col_names)
    data = zeros([15, len(se)*mc_eq])
    m=-1
    for i_eq in trange(n_col,desc='Tsunami+Tide+SLR in {t1}'.format(t1=year),disable=mode_l):
        te = tsu_d[col_names[i_eq]].values
        dummy = []
        for i_eq1 in range(len(te)):
#            for i_seal in range(len(se)):
            dummy.append((te[i_eq1]+se[:]))
        dummy = array(dummy)
        #print(shape(dummy))
        d_nx,d_ny=shape(dummy)
        dummy=dummy.reshape(d_nx*d_ny)
        data[i_eq,:] = dummy[:]
    data_frame= pd.DataFrame(transpose(data),columns=["8.0","8.1","8.2","8.3","8.4","8.5","8.6","8.7","8.8","8.9","9.0","9.1","9.2","9.3","9.4"]) 
    return data_frame




def read_data(fname):
    from numpy import loadtxt,shape,linspace
    array1 = loadtxt(fname)
    return array1

def read_eq_data(water_depth,nrep):
    from numpy import shape,array,sum, linspace,asarray,argwhere,zeros,concatenate,where
    from numpy import random
    import pandas as pd
   # from scipy.interpolate import interp1d
    fname = "./RCP85WA_p50_2000_a.dat"
    max_data = read_data(fname)
    max_data = (max_data - water_depth)
    dummy_x, dummy_y = shape(max_data)
    n_max_data = zeros([dummy_x*nrep,dummy_y])
    for i in range(dummy_y):
        mm=-1
        for jj in range(nrep):
            for j in range(dummy_x):
                mm+=1
                n_max_data[mm,i] = (max_data[j,i] 
    eq = linspace(8.0,9.4,15)
    cols = []
    for i in range(len(eq)):
        cols.append(str(eq[i]))
    dummy_eq = pd.DataFrame(n_max_data, columns=cols)

    return dummy_eq


def smooth(x,window_len=10,window='hanning'):
    #from numpy import *
    #import numpy
    from numpy import r_, convolve,ones, hanning
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len<3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
    s=r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=ones(window_len,'d')
    else:
        w=eval(window+'(window_len)')
    y=convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]

def create_tsunami_mgsep_tide(tsu_df,newtide_data_l,mhw):
    from numpy import array, shape,zeros, append
    import pandas as pd
    col_names = array(list(tsu_df.columns))
    n_rows,n_cols = shape(tsu_df)
    n_tide = len(newtide_data_l)
    data_n = zeros([n_rows*n_tide,n_cols])
    for i_cols in range(n_cols):
        dummy=[]
        tsu_dl = tsu_df.loc[:,col_names[i_cols]]
        for i_tide in range(n_tide):#,desc='Combining EQ {t1}'.format(t1=col_names[i_cols])):
            d1=newtide_data_l[i_tide]+tsu_dl[:]-mhw
            dummy=append(dummy,d1)
        dummy=array(dummy)
        data_n[:,i_cols]=dummy
    tsu_tide = pd.DataFrame(data_n, columns=col_names)
    return tsu_tide

def calculate_flooding(tsu_dd,seal_scen,mode_l):
    import pandas as pd
    if mode_l == False:
        print()
        print('Create linear combination')
    t2000 =adjust_tsunami_dicts1(tsu_dd,seal_scen,'2000',mode_l)
    t2010 =adjust_tsunami_dicts1(tsu_dd,seal_scen,'2010',mode_l)
    t2020 =adjust_tsunami_dicts1(tsu_dd,seal_scen,'2020',mode_l)
    t2030 =adjust_tsunami_dicts1(tsu_dd,seal_scen,'2030',mode_l)
    t2040 =adjust_tsunami_dicts1(tsu_dd,seal_scen,'2040',mode_l)
    t2050 =adjust_tsunami_dicts1(tsu_dd,seal_scen,'2050',mode_l)
    t2060 =adjust_tsunami_dicts1(tsu_dd,seal_scen,'2060',mode_l)
    t2070 =adjust_tsunami_dicts1(tsu_dd,seal_scen,'2070',mode_l)
    t2080 =adjust_tsunami_dicts1(tsu_dd,seal_scen,'2080',mode_l)
    t2090 =adjust_tsunami_dicts1(tsu_dd,seal_scen,'2090',mode_l)
    t2100 =adjust_tsunami_dicts1(tsu_dd,seal_scen,'2100',mode_l)
    year_array=['2000','2010','2020','2030','2040','2050','2060','2070','2080','2090','2100']
    colss = t2000.columns
    d_f={}
    d_f[year_array[0]] = pd.DataFrame(t2000.values,columns=colss)
    d_f[year_array[1]] = pd.DataFrame(t2010.values,columns=colss)
    d_f[year_array[2]] = pd.DataFrame(t2020.values,columns=colss)
    d_f[year_array[3]] = pd.DataFrame(t2030.values,columns=colss)
    d_f[year_array[4]] = pd.DataFrame(t2040.values,columns=colss)
    d_f[year_array[5]] = pd.DataFrame(t2050.values,columns=colss)
    d_f[year_array[6]] = pd.DataFrame(t2060.values,columns=colss)
    d_f[year_array[7]] = pd.DataFrame(t2070.values,columns=colss)
    d_f[year_array[8]] = pd.DataFrame(t2080.values,columns=colss)
    d_f[year_array[9]] = pd.DataFrame(t2090.values,columns=colss)
    d_f[year_array[10]] = pd.DataFrame(t2100.values,columns=colss)
    return d_f
    
def calc_floodheigth_exceedance(d_f,flood_height_l,mode_l):
    import pandas as pd
    from tqdm import trange
    from numpy import linspace, argwhere,interp
    cols=d_f['2000'].columns
    year_array=['2000','2010','2020','2030','2040','2050','2060','2070','2080','2090','2100']
    c_len=len(cols)
    value_501 = []
    for i in trange(len(year_array),desc='Floodheight: {t1:3.2f} m'.format(t1=flood_height_l), disable=mode_l):
        year=year_array[i]
        dd = d_f[year]
        value_dd = []
        for i_eq in range(c_len):
            dummy = dd[cols[i_eq]].values
            value_dd.append(float(len(argwhere(dummy>=flood_height_l)))/float(len(dummy)))
        for i in range(1,len(value_dd)):
            if value_dd[i-1] > value_dd[i]:
                value_dd[i-1] = value_dd[i]
        eq = linspace(8.0,9.4,15)
        eq1=linspace(8.0,9.4,1000)
        y1=smooth(interp(eq1, eq, value_dd),200) 
        index_c = -1
        for i_y1 in range(len(y1)):
            if y1[i_y1]>0.5:
                index_c = i_y1
                break
        if index_c>-1:
            ff = eq1[index_c]
        else:
            ff = float(9.4)
        value_501.append(ff)
    return year_array, value_501


