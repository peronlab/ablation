function [pre_path post_path abl_curation_file_path subvol_restrict] = get_data_path(anim, abl)
%
% "Database" of where to get particular data for a particular animal-ablation type combination
%
%  anim - id string of animal ('274424', for instance)
%  abl - ablation type ('touch','silent','whisking')
%
    settings = get_settings; 

    % this is the naming scheme for animals with a single ablation; if multiple, first will be preablate, then postablate_xxxxx where xxxxx is ablation type
    pre_default = 'session_neuropilone_preablate_merged';
    post_default = 'session_neuropilone_postablate_merged';

    % for animals where, e.g., you imaged L4 and want to exclude
    subvol_restrict = nan; 

    switch anim

        % ----------- NYU ------------
        case '274424' % verified 
            switch abl
                case 'touch'
                    pre_path = [settings.data_path 'an' anim '/' pre_default '/'];
                    post_path = [settings.data_path 'an' anim '/' post_default '/'];
                    abl_curation_file_path = [settings.data_path 'an' anim '/2019_04_07_ablation/manualCurationData.mat'];

                otherwise 
                    disp(sprintf('No %s ablation data for %s' , abl, anim));
            end

        case '272761' % verified
            switch abl
                case 'touch'
                    pre_path = [settings.data_path 'an' anim '/' pre_default '/'];
                    post_path = [settings.data_path 'an' anim '/' post_default '/'];
                    abl_curation_file_path = [settings.data_path 'an' anim '/2019_03_24_ablation/manualCurationData.mat'];

                otherwise 
                    disp(sprintf('No %s ablation data for %s' , abl, anim));
            end

        case '274577' % verified
            switch abl
                case 'touch'
                    pre_path = [settings.data_path 'an' anim '/' pre_default '/'];
                    post_path = [settings.data_path 'an' anim '/' post_default '/'];
                    abl_curation_file_path = [settings.data_path 'an' anim '/2019_04_04_ablation/manualCurationData.mat'];

                otherwise 
                    disp(sprintf('No %s ablation data for %s' , abl, anim));
            end

        case '275798' % verified
            switch abl
                case 'silent'
                    post_default = 'session_neuropilone_postablate_silent_merged';
                    pre_path = [settings.data_path 'an' anim '/' pre_default '/'];
                    post_path = [settings.data_path 'an' anim '/' post_default '/'];
                    abl_curation_file_path = [settings.data_path 'an' anim '/2019_05_07_ablation/manualCurationData.mat'];
                
                case 'whisking'
                    pre_default = 'session_neuropilone_postablate_silent_merged';
                    post_default = 'session_neuropilone_postablate_whisking_merged';
                    pre_path = [settings.data_path 'an' anim '/' pre_default '/'];
                    post_path = [settings.data_path 'an' anim '/' post_default '/'];
                    abl_curation_file_path = [settings.data_path 'an' anim '/2019_05_17_ablation/manualCurationData.mat'];

                otherwise 
                    disp(sprintf('No %s ablation data for %s' , abl, anim));
            end

        case '275801' % verified
            switch abl
                case 'silent'
                    post_default = 'session_neuropilone_postablate_silent_merged';
                    pre_path = [settings.data_path 'an' anim '/' pre_default '/'];
                    post_path = [settings.data_path 'an' anim '/' post_default '/'];
                    abl_curation_file_path = [settings.data_path 'an' anim '/2019_05_21_ablation/manualCurationData.mat'];

                 % a whisking ablation attempt was made prior to touch, after silent ; very few viable targets, only 4 were killed,
                 %   so this data was excluded, but this is why we have postablate_whisking as the pre-ablation directory
                 case 'touch' 
                    pre_default = 'session_neuropilone_postablate_whisking_merged';
                    post_default = 'session_neuropilone_postablate_touch_merged';
                    pre_path = [settings.data_path 'an' anim '/' pre_default '/'];
                    post_path = [settings.data_path 'an' anim '/' post_default '/'];
                    abl_curation_file_path = [settings.data_path 'an' anim '/2019_05_31_ablation/manualCurationData.mat'];

                otherwise 
                    disp(sprintf('No %s ablation data for %s' , abl, anim));
            end

        case '276013' % verified
            switch abl
                case 'silent'
                    post_default = 'session_neuropilone_postablate_silent_merged';
                    pre_path = [settings.data_path 'an' anim '/' pre_default '/'];
                    post_path = [settings.data_path 'an' anim '/' post_default '/'];
                    abl_curation_file_path = [settings.data_path 'an' anim '/2019_06_11_ablation_silent/manualCurationData.mat'];

                case 'whisking'
                    pre_default = 'session_neuropilone_postablate_silent_merged';
                    post_default = 'session_neuropilone_postablate_whisking_merged';
                    pre_path = [settings.data_path 'an' anim '/' pre_default '/'];
                    post_path = [settings.data_path 'an' anim '/' post_default '/'];
                    abl_curation_file_path = [settings.data_path 'an' anim '/2019_06_17_ablation_whisking/manualCurationData.mat'];

                otherwise 
                    disp(sprintf('No %s ablation data for %s' , abl, anim));
            end


        case '278288' % verified
            switch abl 
                case 'whisking'
                    pre_path = [settings.data_path 'an' anim '/' pre_default '/'];
                    post_path = [settings.data_path 'an' anim '/' post_default '/'];
                    abl_curation_file_path = [settings.data_path 'an' anim '/2019_07_08_ablation_whisking/manualCurationData.mat'];

                otherwise 
                    disp(sprintf('No %s ablation data for %s' , abl, anim));
            end

        case '278759' % verified
            switch abl
                case 'whisking'
                    pre_path = [settings.data_path 'an' anim '/' pre_default '/'];
                    post_path = [settings.data_path 'an' anim '/' post_default '/'];
                    abl_curation_file_path = [settings.data_path 'an' anim '/2019_07_17_whisking_ablation/manualCurationData.mat'];

                otherwise 
                    disp(sprintf('No %s ablation data for %s' , abl, anim));
            end

        % ----------- JRC ------------
        case '250220' % verified
            switch abl
                case 'touch'
                    pre_path = [settings.data_path 'an' anim '/' pre_default '/'];
                    post_path = [settings.data_path 'an' anim '/' post_default '/'];
                    abl_curation_file_path = [settings.data_path 'an' anim '/an' anim '_all.mat'];

                otherwise 
                    disp(sprintf('No %s ablation data for %s' , abl, anim));
            end

       case '257218' % verified
            switch abl
                case 'silent'
                    pre_path = [settings.data_path 'an' anim '/' pre_default '/'];
                    post_path = [settings.data_path 'an' anim '/' post_default '/'];
                    abl_curation_file_path = [settings.data_path 'an' anim '/an' anim '_2014_09_01.mat'];

                otherwise 
                    disp(sprintf('No %s ablation data for %s' , abl, anim));
            end

       case '257220' % verified
            switch abl
                case 'touch'
                    pre_path = [settings.data_path 'an' anim '/' pre_default '/'];
                    post_path = [settings.data_path 'an' anim '/' post_default '/'];
                    abl_curation_file_path = [settings.data_path 'an' anim '/an' anim '_2014_09_17.mat'];

                otherwise 
                    disp(sprintf('No %s ablation data for %s' , abl, anim));
            end

        case '258836' % verified
            switch abl
               case 'silent'
                    post_default = 'session_neuropilone_postablate_silent_merged';
                    pre_path = [settings.data_path 'an' anim '/' pre_default '/'];
                    post_path = [settings.data_path 'an' anim '/' post_default '/'];
                    abl_curation_file_path = [settings.data_path 'an' anim '/an' anim '_2014_10_11.mat'];

               case 'touch'
                    pre_default = 'session_neuropilone_postablate_silent_merged';
                    post_default = 'session_neuropilone_postablate_touch_merged';
                    pre_path = [settings.data_path 'an' anim '/' pre_default '/'];
                    post_path = [settings.data_path 'an' anim '/' post_default '/'];
                    abl_curation_file_path = {[settings.data_path 'an' anim '/an' anim '_2014_10_14.mat'], [settings.data_path 'an' anim '/an' anim '_2014_10_16.mat']};

                otherwise 
                    disp(sprintf('No %s ablation data for %s' , abl, anim));
            end

        case '271211' % verified
            switch abl
                case 'silent'
                    post_default = 'session_neuropilone_postablate_silent_merged';
                    pre_path = [settings.data_path 'an' anim '/' pre_default '/'];
                    post_path = [settings.data_path 'an' anim '/' post_default '/'];
                    abl_curation_file_path = [settings.data_path 'an' anim '/an' anim '_2014_11_12.mat'];

                case 'touch'
                    pre_default = 'session_neuropilone_postablate_silent_merged';
                    post_default = 'session_neuropilone_postablate_touch_merged';
                    pre_path = [settings.data_path 'an' anim '/' pre_default '/'];
                    post_path = [settings.data_path 'an' anim '/' post_default '/'];
                    abl_curation_file_path = [settings.data_path 'an' anim '/an' anim '_2014_11_16.mat'];

                otherwise 
                    disp(sprintf('No %s ablation data for %s' , abl, anim));
            end

        case '278937' % verified
            switch abl
                case 'silent'
                    post_default = 'session_neuropilone_postablate_silent_merged';
                    pre_path = [settings.data_path 'an' anim '/' pre_default '/'];
                    post_path = [settings.data_path 'an' anim '/' post_default '/'];
                    abl_curation_file_path = [settings.data_path 'an' anim '/an' anim '_2015_01_30.mat'];

                case 'whisking'
                    pre_default = 'session_neuropilone_postablate_silent_merged'; 
                    post_default = 'session_neuropilone_postablate_whisking_merged';
                    pre_path = [settings.data_path 'an' anim '/' pre_default '/'];
                    post_path = [settings.data_path 'an' anim '/' post_default '/'];
                    abl_curation_file_path = [settings.data_path 'an' anim '/an' anim '_2015_02_04.mat'];

                otherwise 
                    disp(sprintf('No %s ablation data for %s' , abl, anim));
            end

        case '278939' % verified
            switch abl
                case 'silent'
                    post_default = 'session_neuropilone_postablate_silent_merged';
                    pre_path = [settings.data_path 'an' anim '/' pre_default '/'];
                    post_path = [settings.data_path 'an' anim '/' post_default '/'];
                    abl_curation_file_path = [settings.data_path 'an' anim '/an' anim '_2015_01_30.mat'];

                case 'whisking'
                    pre_default = 'session_neuropilone_postablate_silent_merged';
                    post_default = 'session_neuropilone_postablate_whisking_merged';
                    pre_path = [settings.data_path 'an' anim '/' pre_default '/'];
                    post_path = [settings.data_path 'an' anim '/' post_default '/'];
                    abl_curation_file_path = [settings.data_path 'an' anim '/an' anim '_2015_02_04.mat'];

                otherwise 
                    disp(sprintf('No %s ablation data for %s' , abl, anim));
            end


        case '281915' % verified
            switch abl
                case 'whisking'
                    pre_default = 'session_neuropilone_preablate_merged';
                    post_default = 'session_neuropilone_postablate_whisking_merged';
                    pre_path = [settings.data_path 'an' anim '/' pre_default '/'];
                    post_path = [settings.data_path 'an' anim '/' post_default '/'];
                    abl_curation_file_path = [settings.data_path 'an' anim '/an' anim '_2015_02_23.mat'];

                 case 'touch'
                    pre_default = 'session_neuropilone_postablate_whisking_merged';
                    post_default = 'session_neuropilone_postablate_touch_merged';
                    pre_path = [settings.data_path 'an' anim '/' pre_default '/'];
                    post_path = [settings.data_path 'an' anim '/' post_default '/'];
                    abl_curation_file_path = [settings.data_path 'an' anim '/an' anim '_2015_03_02.mat'];

               otherwise 
                    disp(sprintf('No %s ablation data for %s' , abl, anim));
            end

        otherwise 
            disp(sprintf('No data for animal %s (did you accidentally include the an tag?)', anim));
    end
