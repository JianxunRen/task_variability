clear,clc;

DataPath = '/autofs/space/bidlin9_001/users/CoRR/DataProcessed';
Inpath = '/autofs/space/bidlin9_001/users/CoRR/Variability_4subs/Intra'; % Intrasubject results
OutPath = '/autofs/space/bidlin9_001/users/CoRR/Variability_4subs/Inter';
mkdir(OutPath)
subs = textread(['/autofs/space/bidlin9_001/users/CoRR/Lists/list_name.txt'], '%s');
Netpath = '/autofs/space/bidlin4_001/users/Share/NetworkLabels/Surf_7NetworkLabels/fsaverage3';
Net1_lh = squeeze(load_mgh([Netpath '/lh_network_1.mgh']));
Net1_rh = squeeze(load_mgh([Netpath '/rh_network_1.mgh']));
Net1 = [Net1_lh;Net1_rh];
IndNet = find(Net1 == 0);

num_fs4=2562;
nsess = 8;
nsubs = 4;
InterSimilarity = zeros(nsess, num_fs4);
for i = 1:nsess
    i
    Rmat = zeros(num_fs4, length(IndNet), nsubs);
    for s = 1:nsubs %length(subs)
        sub = subs{s}
        sub_name = [sub '_' num2str(i)];	
        % Load Surface Signal
        FileName=[DataPath '/' sub_name '/surf/lh.' sub_name '_bld004_rest_reorient_skip_faln_mc_g1000000000_bpss_resid_fsaverage6_sm6_fsaverage4.nii.gz'];
        hdr  = MRIread(FileName);
        data_lh = reshape(hdr.vol, [hdr.nvoxels, hdr.nframes]);
        data_fs4_lh= data_lh';

        FileName=[DataPath '/' sub_name '/surf/lh.' sub_name '_bld004_rest_reorient_skip_faln_mc_g1000000000_bpss_resid_fsaverage6_sm6_fsaverage3.nii.gz'];
        hdr  = MRIread(FileName);
        data_lh = reshape(hdr.vol, [hdr.nvoxels, hdr.nframes]);
        data_lh=squeeze(data_lh);

        FileName=[DataPath '/' sub_name '/surf/rh.' sub_name '_bld004_rest_reorient_skip_faln_mc_g1000000000_bpss_resid_fsaverage6_sm6_fsaverage3.nii.gz'];
        hdr  = MRIread(FileName);
        data_rh = reshape(hdr.vol, [hdr.nvoxels, hdr.nframes]);
        data_rh=squeeze(data_rh);
        data_tmp=[data_lh;data_rh]';
        
        data_fs4_all = single(data_fs4_lh);
        Data_All = single(data_tmp(:, IndNet));
        Rmat(:,:,s) = my_corr(data_fs4_all, Data_All);
    end 
    Rmat(isnan(Rmat)) = 0;
    count = 0;
    AveRmat = zeros(num_fs4, 1)';
    for m = 1:nsubs %length(subs)
        for n = m+1:nsubs %length(subs)
            count = count + 1;
            tmp = my_matcorr(squeeze(Rmat(:,:,m))', squeeze(Rmat(:,:,n))');
            tmp(isnan(tmp)) = 0;
            AveRmat = AveRmat + tmp;
        end
    end
    count
    InterSimilarity(i, :) = AveRmat/count;
    save_mgh(InterSimilarity(i,:), [OutPath '/lh.session' num2str(i) '_similairity_fs4.mgh'],eye(4))
end
meanIntraVariance = mean(InterSimilarity);
save_mgh(meanIntraVariance, [OutPath '/lh.meanSimilarity_across10session_fs4.mgh'],eye(4))

% Regress Out
Variability_norm = zeros(nsess,num_fs4);
Variability = zeros(nsess,num_fs4);
for i = 1:nsess
    Intra = load_mgh([Inpath '/lh.meanIntravariance_acrossCoRR10_fs4.mgh']);
    Intra_value = Intra;
    Inter = (1 - InterSimilarity(i,:))';
    X = [Intra_value', ones(num_fs4,1)];
    beta = pinv(X)*Inter;
    tmp = Inter - X*beta;
    Variability_norm(i,:) = tmp;
    save_mgh(Variability_norm(i,:), [OutPath '/lh.session' num2str(i) '_intervariability_norm_fs4.mgh'], eye(4))
    tmp = Inter - X(:,1:end-1)*beta(1:end-1);
    Variability(i,:) = tmp;
    save_mgh(Variability(i,:), [OutPath '/lh.session' num2str(i) '_intervariability_fs4.mgh'], eye(4))
end
meanVariability = mean(Variability_norm);
save_mgh(meanVariability, [OutPath '/lh.meanInterVariability_across10session_norm_fs4.mgh'],eye(4))
meanVariability = mean(Variability);
save_mgh(meanVariability, [OutPath '/lh.meanInterVariability_across10session_fs4.mgh'],eye(4))
