
Yeo BT, Krienen FM, Sepulcre J, Sabuncu MR, Lashkari D, Hollinshead M, Roffman JL, Smoller JW, Zollei L., Polimeni JR, Fischl B, Liu H, Buckner RL. The organization of the human cerebral cortex estimated by intrinsic functional connectivity. J Neurophysiol 106(3):1125-65, 2011.

The Yeo split label cortical atlas in MNI152 space can be downloaded at the following link:

https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/Yeo_JNeurophysiol11_SplitLabels

Directory structure should correspond to the following:

Yeo_JNeurophysiol11_SplitLabel/:

    ./Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask+orig.BRIK
    ./Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask+orig.HEAD
    ./Yeo_JNeurophysiol11_SplitLabels_README
    ./Yeo_JNeurophysiol11_SplitLabels_readme.txt

./MNI152/:
    17Networks_ColorLUT_freeview.txt
    17Networks_ColorLUT_fslview.lut
    7Networks_ColorLUT_freeview.txt
    7Networks_ColorLUT_fslview.lut
    FSL_MNI152_FreeSurferConformed_1mm.nii.gz
    MNI152_T1_2mm_brain.nii.gz
    Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask+orig.BRIK
    Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask+orig.HEAD
    Yeo2011_17Networks_N1000.split_components.FSL_MNI152_2mm.nii.gz
    Yeo2011_17Networks_N1000.split_components.FSL_MNI152_FreeSurferConformed_1mm.nii
    Yeo2011_17Networks_N1000.split_components.FSL_MNI152_FreeSurferConformed_1mm.nii.gz
    Yeo2011_7Networks_N1000.split_components.FSL_MNI152_2mm.nii.gz
    Yeo2011_7Networks_N1000.split_components.FSL_MNI152_FreeSurferConformed_1mm.nii.gz

./fsaverage5/:
    label/
    mri/
    scripts/
    surf/
    lh.reg.template.tif
    rh.reg.template.tif
    stats/
    tmp/

./scripts/:
    DownsampleMNI1mmParcellationTo2mm.csh
    ProjectSplitLabels2MNI1mm.m
    lh.Yeo2011_17Networks_N1000.split_components.txt
    lh.Yeo2011_7Networks_N1000.split_components.txt
    rh.Yeo2011_17Networks_N1000.split_components.txt
    rh.Yeo2011_7Networks_N1000.split_components.txt
