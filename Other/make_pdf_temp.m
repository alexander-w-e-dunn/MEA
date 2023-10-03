% combine figures (PNGs) produced by main script to a pdf
% create a folder and save the pdf and PNGs in there for each recording
%   may require moving the PNGs already created

% this script requires installing ghostscript:
% https://www.ghostscript.com/download.html
% this script also require the export_fig package
% https://uk.mathworks.com/matlabcentral/fileexchange/23629-export_fig

% inputs: files; spike method and parameter
clear all
files = dir('*slice*.mat'); 
files = files(~contains({files.name}, 'TTX','IgnoreCase',true));
files = files(~~contains({files.name}, 'FTD'));
files = files(~contains({files.name}, 'spine'));
files = files(~contains({files.name}, '190830'));
files = files(~contains({files.name}, '_2min'));
files = files([16 22 28 30]);

files = files(~contains({files.name}, 'Spikes'));
files = files(~contains({files.name}, 'stim'));
files = files(~contains({files.name}, 'cleaned'));
files = files(~contains({files.name}, 'Filtd'));

method      = 'mSpikes_'; %use mSpikes_, cSpikes_ or aSpikes (no _ if aSpikes)
parameter   = '4'; % put value in inv. commas ('')

progressbar('files','figures')
pb_pos = gcf;
pb_pos.Position = [0.7200 0.4738 0.2400 0.0525];
warning('off','MATLAB:print:FigureTooLargeForPage');
for file = 1:length(files)
    filename = files(file).name;
    
    fig_names = dir(strcat(filename(1:end-4),'*.PNG'));
    fig_names = fig_names(~~contains({fig_names.name},strcat(method,parameter,'_')));
    fig_names = fig_names(~contains({fig_names.name},'based_on'));
    
    for fig = 1:length(fig_names)
        figname = fig_names(fig).name;
        if ~exist(strcat(figname(1:end-4),'.pdf'))
            figure
            imshow(imread(figname));
            h=gcf;
            set(h,'PaperPositionMode','auto');
            set(h,'PaperOrientation','landscape');
            set(h,'Position',[50 50 1200 800]);
            print(gcf, '-dpdf', strcat(figname(1:end-4),'.pdf'))
            
            close(h)
        end
        progressbar(file/length(files),fig/length(fig_names))
    end
    
    % append figures and save
    
    % create folder for the recording and move all files into that
    if ~exist(filename(1:end-4), 'dir')
        
        %put pdf figures in list of cells
        fig_PDFnames = dir(strcat(filename(1:end-4),'*.pdf'));
        
        for figname = 1:length(fig_PDFnames)
            fig_PDFnames_list{figname,1} = fig_PDFnames(figname).name;
        end
        
        mkdir(filename(1:end-4))
        dir1_made = 1;
        try
            append_pdfs(strcat(filename(1:end-4),'_',method,parameter,'_all_figs.pdf'),fig_PDFnames_list{:});
        catch
            disp('already done')
        end
        
        files2move = dir(strcat(filename(1:end-4),'*.p*'));
        files2move = files2move(~~contains({files2move.name},strcat(method,parameter,'_')));
        
        cd(filename(1:end-4))
        mkdir(strcat(filename(1:end-4),'_',method,parameter))
        dir2_made = 1;
        
        cd ..
        
        for file_m = 1:length(files2move)
            movefile(files2move(file_m).name,strcat(filename(1:end-4),'/',filename(1:end-4),'_',method,parameter));
        end
        
    else
        
        cd(filename(1:end-4))
        
        if ~exist(strcat(filename(1:end-4),'_',method,parameter), 'dir')
            mkdir(strcat(filename(1:end-4),'_',method,parameter))
            dir2_made = 1;
            
            cd ..
            
            if ~~exist(strcat(filename(1:end-4),'_',method,parameter,'_all_figs.pdf'));
                delete(strcat(filename(1:end-4),'_',method,parameter,'_all_figs.pdf'))
            end
            
            fig_PDFnames = dir(strcat(filename(1:end-4),'*.pdf'));
            
            for figname = 1:length(fig_PDFnames)
                fig_PDFnames_list{figname,1} = fig_PDFnames(figname).name;
            end
            
            try
                append_pdfs(strcat(filename(1:end-4),'_',method,parameter,'_all_figs.pdf'),fig_PDFnames_list{:});
            catch
                disp('already done')
            end
            
            files2move = dir(strcat(filename(1:end-4),'*.p*'));
            files2move = files2move(~~contains({files2move.name},strcat(method,parameter,'_')));
            
            for file_m = 1:length(files2move)
                movefile(files2move(file_m).name,strcat(filename(1:end-4),'/',filename(1:end-4),'_',method,parameter));
            end
            
        else
            
            cd ..
            
            files2move = dir(strcat(filename(1:end-4),'*.p*'));
            files2move = files2move(~~contains({files2move.name},strcat(method,parameter,'_')));
            
            for file_m = 1:length(files2move)
                movefile(files2move(file_m).name,strcat(filename(1:end-4),'/',filename(1:end-4),'_',method,parameter));
            end
            
            cd(filename(1:end-4))
            cd(strcat(filename(1:end-4),'_',method,parameter))
            
            to_delete = strcat(filename(1:end-4),'_',method,parameter,'_all_figs.pdf');
            delete(to_delete)
            
            fig_PDFnames = dir(strcat(filename(1:end-4),'*.pdf'));
            
            for figname = 1:length(fig_PDFnames)
                fig_PDFnames_list{figname,1} = fig_PDFnames(figname).name;
            end
            
            try
                append_pdfs(strcat(filename(1:end-4),'_',method,parameter,'_all_figs.pdf'),fig_PDFnames_list{:});
            catch
                disp('already done')
            end
            
            cd ..
            cd ..
            
        end
    end
    
    progressbar(file/length(files),fig/length(fig_names))
    clear fig_names fig_PDFnames fig_PDFnames_list to_delete ...
        dir1_made dir2made
end
warning('on','MATLAB:print:FigureTooLargeForPage');
