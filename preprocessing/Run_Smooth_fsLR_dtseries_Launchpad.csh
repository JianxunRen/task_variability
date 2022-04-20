#! /bin/csh

set InPath =  #Change this
set OutPath = /autofs/space/bidlin9_001/users/HCP_FIX/DataProcessed_FIX_s4/Q3
set Downsamplefolder = /autofs/space/bidlin9_001/users/HCP_FIX/DataProcessed_FIX_s4/Q3
set TemplatePath = # midthickness
set count = 1
set stop = 100 #change

set FWHM = 4
set sigma = 1.47 #FWHM = 4, sigma = FHWM/2.355
set TemporalFilter = 200
set att_file = list_100Unrelated.txt
set runs = (LR RL)
set taskn = (LANGUAGE)

while($count <= $stop)
    set sub = `head -n $count $att_file | tail -n 1 | awk '{print $1}'`
    echo "${count}:${sub}"
	foreach run($runs)
	set filename = tfMRI_${taskn}_$run_Atlas
	set cmd = "wb_command -cifti-smoothing ${InPath}/${sub}/${filename}.dtseries.nii ${sigma} ${sigma} COLUMN ${InPath}/${sub}/${filename}_s${FWHM}.dtseries.nii -left-surface $Downsamplefolder/${sub}/MNINonLinear/fsaverage_LR32k/${sub}.L.midthickness.32k_fs_LR.surf.gii -right-surface $Downsamplefolder/${sub}/MNINonLinear/fsaverage_LR32k/${sub}.R.midthickness.32k_fs_LR.surf.gii"
	eval $cmd
	
	wb_command -cifti-convert -to-nifti ${InPath}/${sub}/${filename}_s${FWHM}.dtseries.nii ${InPath}/${sub}/${filename}_s${FWHM}_FAKENIFTI.nii.gz
	fslmaths ${InPath}/${sub}/${filename}_s${FWHM}_FAKENIFTI.nii.gz -Tmean ${InPath}/${sub}/${filename}_s${FWHM}_FAKENIFTI_mean.nii.gz
	fslmaths ${InPath}/${sub}/${filename}_s${FWHM}_FAKENIFTI.nii.gz -bptf `echo "0.5 * $TemporalFilter / 0.72" | bc -l` 0 -add ${InPath}/${sub}/${filename}_s${FWHM}_FAKENIFTI_mean.nii.gz ${InPath}/${sub}/${filename}_s${FWHM}_FAKENIFTI.nii.gz
	wb_command -cifti-convert -from-nifti ${InPath}/${sub}/${filename}_s${FWHM}_FAKENIFTI.nii.gz ${InPath}/${sub}/${filename}_hp200_s${FWHM}.dtseries.nii
	rm ${InPath}/${sub}/${filename}_s${FWHM}_FAKENIFTI.nii.gz ${InPath}/${sub}/${filename}_s${FWHM}_FAKENIFTI_mean.nii.gz
	end

  @ count = $count + 1
end



