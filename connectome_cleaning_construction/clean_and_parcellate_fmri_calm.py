# this script is a function to clean and parcellate fmriprep-processed functional connectomes from the baseline and
# longitudinal time points in the calm (centre for attention, learning, and memory) data set. the function takes 2
# arguments: data set, and participant id. note that loading the python packages requires activating the neuroconda v2
# environment.
def clean_and_parcellate_fmri_calm(dataset, participant_id):
    # step 1 - set up the work space
    from nilearn import image as nimg, input_data
    from nilearn.connectome import ConnectivityMeasure
    import numpy as np
    import h5py
    import pandas as pd
    # step 2 - load participant confound and bold files
    # set the path to the files depending on input data set
    if dataset == 'baseline_calm':
        func_path = '/imaging/projects/external/nkir/analyses/calm-I/fmriprep/outdir/final_output/'
    else:
        func_path = '/imaging/projects/external/nkir/analyses/calm-II/fmriprep/outdir/final_output/'
    # load the confound and bold files for this participant
    confound_file = pd.read_csv(
        func_path + participant_id + '/' + participant_id + '_task-rest_desc-confounds_timeseries.tsv', delimiter='\t')
    bold_img = \
        nimg.load_img(func_path + participant_id + '/' +
                      participant_id + '_task-rest_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii.gz')
    # extract the confounds which we will control for in the bold signal. based on evidence from prior work (parkes et
    # al., 2018), we regress out the first derivatives of cerebrospinal fluid (csf) and white matter (wm) from the
    # aroma-denoised signal, as this combination was shown to perfom well across several fc quality-control (qc)
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
    # step 3 - for each parcellation, clean functional connectomes, correlate, and apply a fisher transformation
    for parcellation in parcellation_list:
        # load the parcellation nifti file
        if parcellation == 'aal116':
            parcellation_img = nimg.load_img('/imaging/astle/am10/atlases/aal_SPM12/aal/atlas/AAL.nii')
            nroi = 116
        elif parcellation == 'schaefer100x7':
            parcellation_img = nimg.load_img(
                '/imaging/astle/am10/atlases/schaefer_2018/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.nii.gz')
            nroi = 100
        elif parcellation == 'schaefer200x7':
            parcellation_img = nimg.load_img(
                '/imaging/astle/am10/atlases/schaefer_2018/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_1mm.nii.gz')
            nroi = 200
        elif parcellation == 'brainnetome246':
            parcellation_img = nimg.load_img('/imaging/astle/am10/atlases/brainnetome/BN_Atlas_246_1mm.nii.gz')
            nroi = 246
        else:
            parcellation_img = nimg.load_img(
                '/imaging/astle/am10/atlases/schaefer_2018/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm.nii.gz')
            nroi = 400
        # create a masker object which will parcellation, clean, and average the functional image to create a functional
        # connectivity matrix! we use a low-pass filter of .1Hz, in line with prior calm (Jones et al., 2021) and NKI
        # (Tobe et al., 2022) studies, and set t_r (repetition time) as 2.25 seconds.
        masker = input_data.NiftiLabelsMasker(labels_img=parcellation_img, standardize=True, memory='nlearn_cache',
                                              verbose=0, detrend=True, low_pass=.1, high_pass=.01, t_r=2.25)
        # create an output array for the correlation matrix
        cleaned_parcellated_correlated_connectome = np.zeros([nroi, nroi])
        # apply masker to clean and average the BOLD data
        cleaned_and_averaged_time_series = masker.fit_transform(bold_img, extracted_confounds)
        # find pearson correlations between pairs of regions
        correlation_matrix = np.squeeze(correlation_measure.fit_transform([cleaned_and_averaged_time_series]))
        # set the diagonal to 0
        np.fill_diagonal(correlation_matrix, 0)
        # apply fisher's z-transformation to edge weights and assign to output array
        cleaned_parcellated_correlated_connectome = np.arctanh(correlation_matrix)
        # and append this output to the output list
        cleaned_parcellated_correlated_connectome_list.append(cleaned_parcellated_correlated_connectome)
        # clear up the work space
        del cleaned_parcellated_correlated_connectome
    # step 4 - save output and clean work space
    save_file_name = func_path + '/cleaned_connectomes/' + participant_id + '_cleaned_functional_connectomes.h5'
    hf = h5py.File(save_file_name, 'w')
    for parcellation_idx in range(0, len(parcellation_list)):
        hf.create_dataset(parcellation_list[parcellation_idx],
                          data=cleaned_parcellated_correlated_connectome_list[parcellation_idx])
    hf.close()
    # and update user
    print("processed functional connectomes for {} in the {} data set\n".format(participant_id, dataset))
