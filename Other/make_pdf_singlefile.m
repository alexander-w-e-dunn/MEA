% append a specific recording's pdf figures
% first move the pdfs to chosen directory

% cd to chosen directory
cd 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'

filename = '191210_slice4_DIV_g04_2018_secondarea_mSpikes_3.mat'
method      = 'mSpikes_'; %use mSpikes_, cSpikes_L or aSpikes (no _ if aSpikes)
parameter   = '3'; % put value in inv. commas ('')

fig_PDFnames = dir(strcat(filename(1:end-4),'*.pdf'));
        for figname = 1:length(fig_PDFnames)
            fig_PDFnames_list{figname,1} = fig_PDFnames(figname).name;
        end
        
        try
            append_pdfs(strcat(filename(1:end-4),'_',method,parameter,'_all_figs.pdf'),fig_PDFnames_list{:});
        catch
            disp('already done')
        end

        % now delete the old pdfs from this directory if copied/duplicated

