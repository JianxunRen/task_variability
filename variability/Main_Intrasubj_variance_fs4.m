clear,clc;

DataPath = '/autofs/space/bidlin9_001/users/CoRR/DataProcessed';
OutPath = '/autofs/space/bidlin9_001/users/CoRR/Variability_4subs/Intra';
mkdir(OutPath)
subs = textread(['/autofs/space/bidlin9_001/users/CoRR/Lists/list_name.txt'], '%s');
Netpath = '.';
Net1_lh = squeeze(load_mgh([Netpath '/lh_network_1.mgh']));
Net1_rh = squeeze(load_mgh([Netpath '/rh_network_1.mgh']));
Net1 = [Net1_lh;Net1_rh];
IndNet = find(Net1 == 0);

num_fs4=2562;

nsubs = 4;
nsess = 8;
IntraVariance = zeros(nsubs, num_fs4);
for s = 1:nsubs %length(subs)
    sub = subs{s}
    Rmat = zeros(num_fs4, length(IndNet), nsess);
    for i = 1:nsess
	    i
        sub_name = [sub '_' num2str(i)];	
        % Load Surface Signal
        FileName=[DataPath '/' sub_name '/surf/lh.' sub_name '_bld004_rest_reorient_skip_faln_mc_g1000000000_bpss_resid_fsaverage6_sm6_fsaverage4.nii.gz'];
        hdr  = MRIread(FileName);
        data_lh = reshape(hdr.vol, [hdr.nvoxels, hdr.nframes]);
        data_fs4_lh= data_lh';
        
        FileName=[DataPath '/' sub_name '/surf/rh.' sub_name '_bld004_rest_reorient_skip_faln_mc_g1000000000_bpss_resid_fsaverage6_sm6_fsaverage4.nii.gz'];
        hdr  = MRIread(FileName);
        data_rh = reshape(hdr.vol, [hdr.nvoxels, hdr.nframes]);
        data_fs4_rh= data_rh';

        FileName=[DataPath '/' sub_name '/surf/lh.' sub_name '_bld004_rest_reorient_skip_faln_mc_g1000000000_bpss_resid_fsaverage6_sm6_fsaverage3.nii.gz'];
        hdr  = MRIread(FileName);
        data_lh = reshape(hdr.vol, [hdr.nvoxels, hdr.nframes]);
        data_lh=squeeze(data_lh);

        FileName=[DataPath '/' sub_name '/surf/rh.' sub_name '_bld004_rest_reorient_skip_faln_mc_g1000000000_bpss_resid_fsaverage6_sm6_fsaverage3.nii.gz'];
        hdr  = MRIread(FileName);
        data_rh = reshape(hdr.vol, [hdr.nvoxels, hdr.nframes]);
        data_rh=squeeze(data_rh);
        data_tmp=[data_lh;data_rh]';
        
        data_fs4_all = single([data_fs4_lh, data_fs4_rh]);
        Data_All = single(data_tmp(:, IndNet));
        Rmat(:,:,i) = my_corr(data_fs4_all, Data_All);
    end 
    Rmat(isnan(Rmat)) = 0;
    count = 0;
    AveRmat = zeros(num_fs4, 1);
    for m = 1:nsess
        for n = m+1:nsess
            count = count + 1;
            tmp = my_corr(squeeze(Rmat(:,:,m))', squeeze(Rmat(:,:,n))');
            tmp(isnan(tmp)) = 0;
            AveRmat = AveRmat + diag(tmp);
        end
    end
    count
    IntraVariance(s, :) = 1 - AveRmat/count;
    %save_mgh(IntraVariance(s,:), [OutPath '/lh.' sub '_intravariance_fs4.mgh'],eye(4))
end
meanIntraVariance = mean(IntraVariance);
save_mgh(meanIntraVariance, [OutPath '/lh.meanIntravariance_acrossCoRR10_fs4.mgh'],eye(4))
