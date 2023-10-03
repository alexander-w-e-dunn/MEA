
edge_thresh = 0.6;
node_thresh = edge_thresh;
scale = 'rate';
cblim = [0 30];

filename = '200114_FTDOrg_GrpB_2A_Slice8_mSpikes_3_adjM_0.175.mat'
directory = 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output';
    plot_FR_edgenetwork_fcn(directory, filename, edge_thresh, node_thresh, scale, cblim);

    filename = '200114_FTDOrg_GrpC_8B_Slice5_mSpikes_3_adjM_0.175.mat'
        plot_FR_edgenetwork_fcn(directory, filename, edge_thresh, node_thresh, scale, cblim);

    filename = '200114_FTDOrg_GrpA_1A_Slice9_mSpikes_3_adjM_0.175.mat'
        plot_FR_edgenetwork_fcn(directory, filename, edge_thresh, node_thresh, scale, cblim);

            filename = '200114_FTDOrg_GrpD_2A_Slice7_mSpikes_3_adjM_0.175.mat'
        plot_FR_edgenetwork_fcn(directory, filename, edge_thresh, node_thresh, scale, cblim);
