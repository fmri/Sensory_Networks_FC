######################################################################################
'''
#Create Design Matrix
'''
######################################################################################
'''
print("Now creating the design matrix")

for subjid in subjects:
    for parcel_num in sig_parcels:
        # Remove newline character
        subjid = subjid.strip()

        # read in the SEED txt file in array
        seed_run_1 = np.loadtxt("/home/aamoore/scc/projectnb/perlab/projects/LLFC/parcels/EPI_parcel_Timeseries/{}/run_1/froi-parcel{}_avg_thr90_100_ts.txt".format(subjid,parcel_num))
        seed_run_2 = np.loadtxt("/home/aamoore/scc/projectnb/perlab/projects/LLFC/parcels/EPI_parcel_Timeseries/{}/run_2/froi-parcel{}_avg_thr90_100_ts.txt".format(subjid,parcel_num))
        seed = np.append(seed_run_1, seed_run_2)
        
        # Load design.mat in a dataframe
        design_file_run_1 = ("/home/aamoore/scc/projectnb/perpwm/fmri_analysis/experiments/LangLoc/{}/model/smoothed/run_1/design.mat").format(subjid)
        design_file_run_2 = ("/home/aamoore/scc/projectnb/perpwm/fmri_analysis/experiments/LangLoc/{}/model/smoothed/run_2/design.mat").format(subjid)
        
        design_df_run_1 = pd.read_csv(design_file_run_1, sep='\s+', skiprows=5, header=None, usecols=[0,1])
        design_df_run_2 = pd.read_csv(design_file_run_2, sep='\s+', skiprows=5, header=None, usecols=[0,1])

        # Scale, Center and Demean arrays
        def scaleCenterDemean(array):
            scaled = array * 1 / (np.max(array) - np.min(array))
            centered = scaled - np.min(scaled) - 0.5
            demeaned = centered - np.mean(centered)
            return demeaned        
        def rmsToReference(reference,target):
            # Calculate RMS of the reference array
            rms_reference_value = np.sqrt(np.mean(reference**2))
            rms_target_intermediate = np.sqrt(np.mean(target**2))
            # Normalize the target array based on the RMS of the reference array
            rms_target = target * (rms_reference_value / rms_target_intermediate) # flipping based on tyler's change
            #rms_target = target * (rms_target_intermediate / rms_reference_value)
            return rms_target   
            
        # Scale, Center and Demean arrays 
        intact_1 = scaleCenterDemean(design_df_run_1.iloc[:,1])
        intact_2 = scaleCenterDemean(design_df_run_2.iloc[:,1])
        
        degr_1 = scaleCenterDemean(design_df_run_1.iloc[:,0])
        degr_2 = scaleCenterDemean(design_df_run_2.iloc[:,0])
        
        seed_1_sCD = scaleCenterDemean(seed_run_1)
        seed_2_sCD = scaleCenterDemean(seed_run_2)
        
        
        # RMS seed to int and degr regressors
        seed_1_rms_int = rmsToReference(intact_1,seed_1_sCD)
        seed_2_rms_int = rmsToReference(intact_2,seed_2_sCD)
        
        seed_1_rms_degr = rmsToReference(degr_1,seed_1_sCD)
        seed_2_rms_degr = rmsToReference(degr_2,seed_2_sCD)

        # Shift the Seed 
        seed_shifted_1_int = seed_1_rms_int - np.min(seed_1_rms_int)
        seed_shifted_2_int = seed_2_rms_int - np.min(seed_2_rms_int) 
        seed_shifted_1_degr = seed_1_rms_degr - np.min(seed_1_rms_degr)
        seed_shifted_2_degr = seed_2_rms_degr -  np.min(seed_2_rms_degr)
        
        
        # Do the math
        # make a 2x1 array for each type to run through np.prod()
        intxseed_run_1 = np.array([seed_shifted_1_int, intact_1])
        intxseed_run_2 = np.array([seed_shifted_2_int, intact_2])
        degrxseed_run_1 = np.array([seed_shifted_1_degr, degr_1])
        degrxseed_run_2 = np.array([seed_shifted_2_degr, degr_2])
        
        
        # Take the product for each time point across condition and seed and concatenate
        # Intact
        intact_seed_run_1 = np.prod(intxseed_run_1, axis=0) 
        intact_seed_run_2 = np.prod(intxseed_run_2, axis=0)
        intact_seed = np.append(intact_seed_run_1, intact_seed_run_2)

        
        # Degraded
        degraded_seed_run_1 = np.prod(degrxseed_run_1, axis=0)
        degraded_seed_run_2 = np.prod(degrxseed_run_2, axis=0)
        degraded_seed = np.append(degraded_seed_run_1,degraded_seed_run_2)

        # Save seed only
        seed = np.append(seed_1_rms_int, seed_2_rms_int) 
        
        #Final scale product arrays
        intact_seed = rmsToReference(intact_1, intact_seed)
        degraded_seed = rmsToReference(degr_1, degraded_seed)
        
        #Add Intact and Degraded conditions without the seeds
        Intact = np.append(intact_1, intact_2)
        Degraded = np.append(degr_1, degr_2)
        # Runs
        # make a  column with the correct run number and length
        run_1 = ['1']*len(seed_run_1)
        run_2 = ['2']*len(seed_run_2)
        run = np.append(run_1,run_2)
        
        
    
    
        # Create dataframe for regressor file
        regressor_df = pd.DataFrame({'seed': seed, 'int_seed':intact_seed, 'degr_seed':degraded_seed, 'Intact': Intact, 'Degraded': Degraded,  'run': run}).reindex(columns=['degr_seed','int_seed','seed','run'])

        # Make output path I can change this later to move them into the correct folders usually "brains/{subjid}/design"
        #output_dir = "/home/aamoore/scc/projectnb/perpwm/brains/{subjid}/design/".format(subjid=subjid)
        output_dir = "/home/aamoore/scc/projectnb/perlab/projects/LLFC/brains/{}/design/".format(subjid)
        #if not os.path.exists(output_dir):
        #    os.makedirs(output_dir)
        print("Saving parcel {} regressor for subject {}".format(parcel_num,subjid))
        regressor_df.to_csv(os.path.join(output_dir, "gPPI_seed_parcel{}_regressor.csv".format(parcel_num)), index=False)


############################################
Next Function
############################################

def createBetaContrastAvgMatrixBothRuns(contrast,i):
    altmodel = "parcel" #"LLCB"
    pivot_tables = []
    runs = ["1","2"]
    for run in runs:
        appended_dict_df = pd.DataFrame()
        all_index = pd.DataFrame()
        for subject in subjects:
            # Remove newline character
            subject = subject.strip()
            print(subject)
            subject_beta_df = pd.DataFrame()
            for index, s in sig_parcels_and_names.iterrows():
                seed = s['Parcel_Number']
                all_dict = {"Seed_X":[],"Target_Y":[],f"Beta_{subject}":[]}
                design = pd.read_csv(f"/projectnb/perlab/projects/LLFC/analysis/fROI-{altmodel}{seed}/{subject}/model/smoothed/run_{run}/design.mat",header=None, skiprows=5,delimiter='\t')
                # column 4 is the seed column i (2) is the int interaction
                X = design
                for index, t in sig_parcels_and_names.iterrows():
                    target = t['Parcel_Number']
                    y = pd.read_csv(f"/projectnb/perlab/projects/LLFC/analysis/fROI-{altmodel}{target}/{subject}/model/smoothed/run_{run}/design.mat",header=None, skiprows=5,delimiter='\t')
                    Y = y.iloc[:,4]
                    beta = np.linalg.lstsq(X, Y, rcond=None)[0]
                    all_dict["Seed_X"].append(f"{s['Parcel_Number']}")
                    all_dict["Target_Y"].append(f"{t['Parcel_Number']}")
                    all_dict[f"Beta_{subject}"].append(beta[i])
                    current_dict = pd.DataFrame(all_dict[f"Beta_{subject}"], columns=[f"Beta_{subject}"])
                    current_index = pd.DataFrame({"Seed_X": all_dict["Seed_X"], "Target_Y": all_dict["Target_Y"]})
                all_index = pd.concat([all_index,current_index], axis=0,ignore_index=True)
                subject_beta_df = pd.concat([subject_beta_df,current_dict], axis=0,ignore_index=True)
            appended_dict_df = pd.concat([appended_dict_df,subject_beta_df], axis=1)
            appended_mean = pd.Series(appended_dict_df.mean(axis=1), name=f"Beta_{contrast}")
        df = pd.concat([all_index[:len(appended_dict_df)],appended_mean],axis=1)
        pivot = pd.pivot_table(df, index='Seed_X', columns='Target_Y', values=f"Beta_{contrast}")
        pivot_tables.append(pivot)
    print(f"There are {len(pivot_tables)} runs being averaged.")
    averaged_pivot = (sum(pivot_tables) / len(pivot_tables))
    #averaged_pivot_names_nums = pd.concat(averaged_pivot
    averaged_pivot.to_csv(f"/projectnb/perlab/projects/LLFC/scripts/notebooks/outputs/{contrast}_beta_table_avg_across_runs_nums.csv", sep=",",header=True, index=True)
    #df.to_csv(f"/projectnb/perlab/projects/LLFC/scripts/notebooks/outputs/int_seed_r{run}_df_num.csv", sep=",",header=True, index=True)
    return averaged_pivot
int_only = createBetaContrastAvgMatrixBothRuns("int_seed", 2)
