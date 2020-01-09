function abl_print_fig (fh, output_filename)
%
% So you can do all your prettification of figures in one place
%
    settings = get_settings;
    settings.print_on = 0;

    if (~settings.print_on) ; return ; end
    
    if (~exist(settings.print_dir, 'dir'))
        mkdir(settings.print_dir);
    end
    output_full_path = [settings.print_dir filesep output_filename '.eps'];

    set(fh,'PaperPositionMode','auto');

    ppos = get(fh,'PaperPosition');
    su = get(fh,'Units');
    pu = get(fh,'PaperUnits');  
    set(fh,'Units',pu);
    spos = get(fh,'Position');
    set(fh,'Position',[spos(1) spos(2) ppos(3) ppos(4)])
    set(fh,'Units',su)
    set(fh,'renderer','Painters');

    print (fh, output_full_path ,'-dpsc2', '-painters', '-noui');
