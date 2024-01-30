# this script is a python function to clean and parcellate fmriprep-processed functional connectomes across 3 sessions
# (ses-BAS1, ses-FLU1, and ses-FLU2) from the nki-rockland sample longitudinal discovery of brain development
# trajectories sub-study. the function takes 1 argument:participant id. note that loading the python packages requires
# activating the neuroconda v2 environment.
def clean_and_parcellate_fmri_nki(participant_id):
    # step 1 - set up the work space
    from nilearn import image as nimg, input_data
    from nilearn.connectome import ConnectivityMeasure
    import numpy as np
    import h5py
    import pandas as pd
    from glob import glob
    # step 2 - find how many sessions this participant has
    func_path = '/imaging/projects/external/nkir/analyses/enhanced_nkir/fmriprep/outdir/final_output/'
    sub_session_paths = glob(func_path + participant_id + '/ses-*')
    # extract the sessions
    session_list = []
    for session_idx in range(0, len(sub_session_paths)):
        session_list.append(sub_session_paths[session_idx][98:])
    # step 3 - for each session, extract the bold and confound files!
    for session in session_list:
        confound_file = pd.read_csv(func_path + participant_id + '/' + session + '/' + participant_id + '_' + session
                                    + '_task-rest_acq-1400_desc-confounds_timeseries.tsv', delimiter='\t')
        bold_img = nimg.load_img(func_path + participant_id + '/' + session + '/' + participant_id + '_' + session +
                                 '_task-rest_acq-1400_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii.gz')
        # extract the confounds which we will control for in the bold signal. based on evidence from prior work (parkes
        # et al., 2018), we regress out the first derivatives of cerebrospinal fluid (csf) and white matter (wm) from
        # the aroma-denoised signal, as this combination was shown to perfom well across several fc quality-control (qc)
        # measures,such as qc-fc distance-dependence and loss of temporal degrees of freedom.
        extracted_confounds = confound_file[['csf_derivative1', 'white_matter_derivative1']].values
        # replace nan in extracted_confounds with 0. this corresponds to the first row only
        extracted_confounds[np.isnan(extracted_confounds)] = 0
        # set the measure of connectivity
        correlation_measure = ConnectivityMeasure(kind='correlation')
        # set the list of parcellations we'll use
        parcellation_list = ['aal116', 'schaefer100x7', 'schaefer200x7', 'brainnetome246', 'schaefer400x7']
        # create a list to hold the cleaned and parcellated connectomes
        cleaned_parcellated_correlated_connectome_list = []
        # step 4 - now loop over the parcellations, cleaning and constructing functional connectomes!
        for parcellation in parcellation_list:
            # load the parcellation nifti file
            if parcellation == 'aal116':
                parcellation_img = nimg.load_img('/imaging/astle/am10/atlases/aal_SPM12/aal/atlas/AAL.nii')
            elif parcellation == 'schaefer100x7':
                parcellation_img = nimg.load_img('/imaging/astle/am10/atlases/schaefer_2018/'
                                                 'Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.nii.gz')
            elif parcellation == 'schaefer200x7':
                parcellation_img = nimg.load_img('/imaging/astle/am10/atlases/schaefer_2018/'
                                                 'Schaefer2018_200Parcels_7Networks_order_FSLMNI152_1mm.nii.gz')
            elif parcellation == 'brainnetome246':
                parcellation_img = nimg.load_img('/imaging/astle/am10/atlases/brainnetome/BN_Atlas_246_1mm.nii.gz')
            else:
                parcellation_img = nimg.load_img('/imaging/astle/am10/atlases/schaefer_2018/'
                                                 'Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm.nii.gz')
            # create a masker object which will parcellation, clean, and average the functional image to create a
            # functional connectivity matrix! we use a low-pass filter of .1Hz, in line with prior calm (Jones et al.,
            # 2021) and NKI (Tobe et al., 2022) studies, and set t_r (repetition time) as 1.4 seconds.
            masker = input_data.NiftiLabelsMasker(labels_img=parcellation_img, standardize=True, memory='nlearn_cache',
                                                  verbose=0, detrend=True, low_pass=.1, high_pass=.01, t_r=1.4)
            # apply masker to clean and average the BOLD data
            cleaned_and_averaged_time_series = masker.fit_transform(bold_img, extracted_confounds)
            # find pearson correlations between pairs of regions
            correlation_matrix = np.squeeze(correlation_measure.fit_transform([cleaned_and_averaged_time_series]))
            # set the diagonal to 0
            np.fill_diagonal(correlation_matrix, 0)
            # extract the lower triangle of the array with an offset of 1, apply fisher's z-transformation to edge
            # weights, and assign to output array
            cleaned_parcellated_correlated_connectome = np.arctanh(np.tril(correlation_matrix, k=1))
            # and append this output to the output list
            cleaned_parcellated_correlated_connectome_list.append(cleaned_parcellated_correlated_connectome)
            # clear up the work space
            del cleaned_parcellated_correlated_connectome
        # step 4 - save output and clean work space
        save_file_name = func_path + participant_id + '_' + session + '_cleaned_functional_connectomes.h5'
        hf = h5py.File(save_file_name, 'w')
        for parcellation_idx in range(0, len(parcellation_list)):
            hf.create_dataset(parcellation_list[parcellation_idx],
                              data=cleaned_parcellated_correlated_connectome_list[parcellation_idx])
        hf.close()
        # and update user
        print("processed functional connectomes for {} in {}. \n".format(participant_id, session))
