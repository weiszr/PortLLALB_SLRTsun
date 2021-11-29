

def main_fig3(mc_tide,mcs,mode):
    import pandas as pd
    import scipy.stats
    from tqdm import trange
    from random import random
    import math
    from numpy import linspace, zeros, histogram, round, shape, array
    import pickle
    from tools_main import read_eq_data, create_tsunami_mgsep_tide
    from tools_main import read_sea_level_pd, calculate_flooding
    from tools_main import adjust_tsunami_dicts1
    
    from numpy import ones_like, savetxt, append
    rcp_scenario = ['RCP85NA', 'RCP85WA','RCP26NA', 'RCP26WA'] 
    if mode == False:
        print('\t Flood-height comparison for 2000, 2050,2070, and 2100 for all SLR scenarios')
        print()
        print('Sea-level subsample size:',mcs) 
        print('Tide subsample size:',mc_tide)
        print()
    tsunami_mode = 'tsunami_tide'
    if tsunami_mode == 'tsunami_tide':
        data_df=pd.read_csv('LA_tide_MSL.dat', delimiter = ' ')
        data_arr = []
        data_arr_d = data_df['value'].values
        for ii in range(len(data_arr_d)):
            if math.isnan(data_arr_d[ii]) == False:
                data_arr.append(data_arr_d[ii])
        sample_pdf = scipy.stats.gaussian_kde(data_arr)
        newtide_data = sample_pdf.resample(mc_tide).T[:,0]
        tsu_data = read_eq_data(-0.2,1)
        tsu_data1=create_tsunami_mgsep_tide(tsu_data,newtide_data,0.0)

 
    years = [2000,2050,2070,2100]
    cols = []
    for i in range(len(years)):
        cols.append(str(years[i]))
    d_f={}
    
    seal85NA=read_sea_level_pd(rcp_scenario[0],years,mcs)
    nx1,nt = shape(seal85NA)
    nx,ny = shape(tsu_data1)
    bins = linspace(-1.6,7.11,101)
    for i in range(len(rcp_scenario)):
        if mode == False:
            print('Sea-level Scenario: {t1}'.format(t1=rcp_scenario[i]))
        seal85NA=read_sea_level_pd(rcp_scenario[i],years,mcs)
        df_years = {}
        for j in range(len(years)):
            tt = adjust_tsunami_dicts1(tsu_data1,seal85NA,str(years[j]),mode)
            dummy = tt.values
            nx, ny = shape(dummy)
            dummy = dummy.reshape(nx*ny)
            weights = ones_like(dummy)/float(len(dummy))
            x1_1,y1_1 = histogram(dummy,bins,weights=weights)
            fname1_1 = 'file1_{t1}_{t2}.dat'.format(t1=rcp_scenario[i],t2=years[j])
            savetxt(fname1_1,list(zip(x1_1,y1_1)))
# TODO: 
#       - multithreading vs serial computing
#       - How should data and file for tides should be handled?
#       - How should the data and files for sea-level rise should be handled?
#       - Add file that contain the table with the statistics

def main_floodheight_t(rcp_scenario,tsunami_mode,mc_tide,mcs,flood_height,mode):
    import pandas as pd
    import scipy.stats
    from tqdm import trange, tqdm
    from random import random
    import math
    from numpy import linspace,zeros,histogram,round
    from numpy import shape,argwhere,interp,savetxt,round
    import pickle
    import pandas as pd
    from tools_main import read_eq_data,create_tsunami_mgsep_tide
    from tools_main import read_sea_level_pd, calculate_flooding
    from tools_main import calc_floodheigth_exceedance,smooth
    from tools_main import adjust_tsunami_dicts1
    #rcp_scenario = ['RCP85NA', 'RCP85WA','RCP26NA', 'RCP26WA'] 
    if mode == False:
#         print('\t Flood-height comparison for 2000, 2050,2070, and 2100 for all SLR scenarios')
#         print()
        print('Sea-level subsample size:',mcs) 
        print('Tide subsample size:',mc_tide)
    tsunami_mode = 'tsunami_tide'
    flood_height = round(flood_height,2)
    if mode ==False:
        print('Flood Heights:', flood_height)
        print()
    if tsunami_mode == 'tsunami_tide':
        data_df=pd.read_csv('LA_tide_MSL.dat', delimiter = ' ')
        data_arr = []
        data_arr_d = data_df['value'].values
        for ii in range(len(data_arr_d)):
            if math.isnan(data_arr_d[ii]) == False:
                data_arr.append(data_arr_d[ii])
        sample_pdf = scipy.stats.gaussian_kde(data_arr)
        newtide_data = sample_pdf.resample(mc_tide).T[:,0]
        tsu_data = read_eq_data(-0.2,1)
        tsu_data1=create_tsunami_mgsep_tide(tsu_data,newtide_data,0.0)
    years = linspace(2000,2100,11,dtype=int)
    cols = []
    for i in range(len(years)):
        cols.append(str(years[i]))
    seal85NA=read_sea_level_pd(rcp_scenario,years,mcs)
    nx1,nt = shape(seal85NA)
    nx,ny = shape(tsu_data1)
    if mode == False:
        print('Sea-level Scenario: {t1}'.format(t1=rcp_scenario))
    seal85NA=read_sea_level_pd(rcp_scenario,years,mcs)
    df_years = {}
    value_501 = []
    fl_data = zeros([len(years),len(flood_height)+1])
    fl_data[:,0]=years[:]
    for j in range(len(years)):
        tt = adjust_tsunami_dicts1(tsu_data1,seal85NA,str(years[j]),mode)
        dummy = tt.values
        nx, ny = shape(dummy)
        eq = linspace(8.0,9.4,15)
        eq1=linspace(8.0,9.4,1000)
        jjj=0
        value_dd = []
        for jj in trange(ny*len(flood_height),disable=mode):
            i = jj % ny
#             print(j,jj,i,jjj,flood_height[jjj])
            if jj%ny==0 and jj>0:
                for i_d in range(1,len(value_dd)):
                    if value_dd[i_d-1] > value_dd[i_d]:
                        value_dd[i_d-1] = value_dd[i_d]
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
                fl_data[j,jjj+1]=ff
        #        print(fl_data[j,jjj])
                value_dd = []
                jjj = jjj+1
        #    print(jj,len(value_dd),i,jjj,jjj+1)
            value_dd.append(float(len(argwhere(dummy[:,i]>=flood_height[jjj])))/float(len(dummy[:,i])))
#        for i_d in range(1,len(value_dd)):
        if value_dd[i_d-1] > value_dd[i_d]:
                value_dd[i_d-1] = value_dd[i_d]
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
        fl_data[j,jjj+1]=ff
 
        if mode==False:
            print()
    fname1_1 = 'exe1_{t1}.dat'.format(t1=rcp_scenario)
    print(fname1_1)
    savetxt(fname1_1,fl_data,fmt='%3.2f')
# TODO: 
#       - multithreading vs serial computing
#       - How should data and file for tides should be handled?
#       - How should the data and files for sea-level rise should be handled?




if __name__ == '__main__':
    import argparse
    from numpy import linspace
    parser = argparse.ArgumentParser(description='Sea-level rise, Tsunami and Tides')
    parser.add_argument('-run', '--runmode', help='flood_height,or distri', required=True)
    parser.add_argument('-s','--scenario',help='RCP Scenario',required=False)
    parser.add_argument('-m','--mode',help='Mode (tsunami, tsunami_tide)',required=False)
    parser.add_argument('-sti','--subs_tide',help='Subsample size of tide',required=False)
    parser.add_argument('-sse','--subs_seal',help='Subsample size of sea level',required=False)
    parser.add_argument('-fh','--flood_h',help='Flood Heights',required=False)
    parser.add_argument('-p','--production', action='store_true')
    parser.set_defaults(production=False)
    
    args = parser.parse_args()
    main_mcs=50
    main_mc_tide = 50
    flood_height_main = linspace(0.5,1.5,3)
    if str(args.runmode) != 'None':
        run_mode = str(args.runmode)
    if str(args.scenario) != 'None':
        main_rcp_scenario = str(args.scenario)
    if str(args.mode) != 'None':
        main_mode = str(args.mode)
    if str(args.subs_tide) != 'None':
       main_mc_tide = int(args.subs_tide)
    if str(args.subs_seal) != 'None':
        main_mcs = int(args.subs_seal)
    if str(args.flood_h) != 'None':
        my_list = [float(item) for item in args.flood_h.split(',')]
        flood_height_main = linspace(my_list[0],my_list[1],int(my_list[2]))
    main_prod_mode =args.production
    if main_prod_mode == False:
        print("\t\t \033[1m Sea-level rise, Tsunami and Tides\033[0m")
        print()
    if run_mode == 'flood_height':
        if run_mode != 'None' and str(args.scenario) != 'None':
            if main_prod_mode == False:
                print("\t\t \033[1m Flood-Height Calculation\033[0m")
                print()
            main_floodheight_t(main_rcp_scenario,main_mode,main_mc_tide,main_mcs,flood_height_main,main_prod_mode)
        else:
            print("\t\t \033[1m Flood-Height Exceedance Calculation\033[0m")
            print()
            print('Please choose -s option (rcp scenario) and -m option (tsunami, tsunami_tide)')
            exit()
    if run_mode == 'distribution':
        if main_prod_mode==False:
            print("\t \033[1m Flood Height Distributions\033[0m")
            print()
        main_fig3(main_mc_tide,main_mcs,main_prod_mode)
    if run_mode != 'flood_height' and run_mode != 'distribution':
        print('Not a valid option!')
