clear all
close all
format long
 
%USER PARAMETER

%value of the maximal level for Gauss Lobatto Chebyshev quadratures
lvlmax = 9;
%min/max values for parameter 1 
xmin_param_1 = 0;
xmax_param_1 = 10;
%min/max values for parameter 2
xmin_param_2 = 10;
xmax_param_2 = 20;
%min/max values for parameter 3
xmin_param_3 = -10;
xmax_param_3 = 5;
% !!! TO CLOSE ONE DIMENSION/PARAMETER, JUST PUT THE SAME VALUE FOR XMIN AND
%XMAX !!!

%arrays containing all level of the quadratures for each parameter
all_Level_quadrature_param_1 = zeros(10, 10000);
all_Level_quadrature_param_2 = zeros(10, 10000);
all_Level_quadrature_param_3 = zeros(10, 10000);
%array containing the size of each level of the Gauss Lobatto Chebyshev quadratures
level_size_GLC_quadrature = zeros(50, 1);

%arrays containing all level of the quadratures values for each parameter INLINE
%with their associated level
all_Level_param_1_val_and_lvl_inline = zeros(10000, 2);
all_Level_param_2_val_and_lvl_inline = zeros(10000, 2);
all_Level_param_3_val_and_lvl_inline = zeros(10000, 2);

%GENERATING ALL LEVEL QUADRATURES FOR ALL PARAMETERS

%Gauss Lobatto Chebyshev quadrature for parameter 1

%two first values for level 0
all_Level_param_1_val_and_lvl_inline(1,1) = xmin_param_1;
all_Level_param_1_val_and_lvl_inline(1,2) = 0;
all_Level_param_1_val_and_lvl_inline(2,1) = xmax_param_1;
all_Level_param_1_val_and_lvl_inline(2,2) = 0;
compt = 3;

for lvl = 1 : 1 : lvlmax
    current_quad = GLC_quadrature( xmin_param_1, xmax_param_1, lvl );
    size_of_current_quad = size(current_quad,2);
    level_size_GLC_quadrature(lvl) = size_of_current_quad;
    all_Level_quadrature_param_1(lvl, 1:size_of_current_quad) = current_quad;
    
    %inline storage of the parameter values
    transp_current_quad = current_quad';
    all_Level_param_1_val_and_lvl_inline(compt:compt+size_of_current_quad-1,1) = transp_current_quad;
    %and their associated levels
    all_Level_param_1_val_and_lvl_inline(compt:compt+size_of_current_quad-1,2)=lvl;
    
    compt = compt + size_of_current_quad;
end
%resizing the table
final_size = compt-1;
all_Level_param_1_val_and_lvl_inline(final_size+1:10000,:) = [];


%Gauss Lobatto Chebyshev quadrature for parameter 2

%two first values for level 0
all_Level_param_2_val_and_lvl_inline(1,1) = xmin_param_2;
all_Level_param_2_val_and_lvl_inline(1,2) = 0;
all_Level_param_2_val_and_lvl_inline(2,1) = xmax_param_2;
all_Level_param_2_val_and_lvl_inline(2,2) = 0;
compt = 3;

for lvl = 1 : 1 : lvlmax
    current_quad = GLC_quadrature( xmin_param_2, xmax_param_2, lvl );
    size_of_current_quad = size(current_quad,2);
    all_Level_quadrature_param_2(lvl, 1:size_of_current_quad) = current_quad;
    
    %inline storage of the parameter values
    transp_current_quad = current_quad';
    all_Level_param_2_val_and_lvl_inline(compt:compt+size_of_current_quad-1,1) = transp_current_quad;
    %and their associated levels
    all_Level_param_2_val_and_lvl_inline(compt:compt+size_of_current_quad-1,2)=lvl;
    
    compt = compt + size_of_current_quad;
end
%resizing the table
final_size = compt-1;
all_Level_param_2_val_and_lvl_inline(final_size+1:10000,:) = [];


%Gauss Lobatto Chebyshev quadrature for parameter 3

%two first values for level 0
all_Level_param_3_val_and_lvl_inline(1,1) = xmin_param_3;
all_Level_param_3_val_and_lvl_inline(1,2) = 0;
all_Level_param_3_val_and_lvl_inline(2,1) = xmax_param_3;
all_Level_param_3_val_and_lvl_inline(2,2) = 0;
compt = 3;

for lvl = 1 : 1 : lvlmax
    current_quad = GLC_quadrature( xmin_param_3, xmax_param_3, lvl );
    size_of_current_quad = size(current_quad,2);
    all_Level_quadrature_param_3(lvl, 1:size_of_current_quad) = current_quad;
    
    %inline storage of the parameter values
    transp_current_quad = current_quad';
    all_Level_param_3_val_and_lvl_inline(compt:compt+size_of_current_quad-1,1) = transp_current_quad;
    %and their associated levels
    all_Level_param_3_val_and_lvl_inline(compt:compt+size_of_current_quad-1,2)=lvl;
    
    compt = compt + size_of_current_quad;
end
%resizing the table
final_size = compt-1;
all_Level_param_3_val_and_lvl_inline(final_size+1:10000,:) = [];



%RESEARCH OF ALL COMBO FOR THE WANTED LEVEL
% research of all combinations for the wanted level
% !!! 3 DIMENSION/PARAMETER
if(xmin_param_1 ~= xmax_param_1 && xmin_param_2 ~= xmax_param_2 && xmin_param_3 ~= xmax_param_3)
    combos = zeros(final_size^2,3);
    i_combo = 1;
    nb_combos_found = 0;

    for(p1 = 1 : 1 : final_size)
        for(p2 = 1 : 1 : final_size)
            for(p3 = 1 : 1 : final_size)
                %scanning all combo and compute their cumulated levels
                combo_lvl = all_Level_param_1_val_and_lvl_inline(p1,2)+ all_Level_param_2_val_and_lvl_inline(p2,2)+all_Level_param_3_val_and_lvl_inline(p3,2);
                %adding combo if it matches with the wanted level
                if(combo_lvl == lvlmax)
                    combos(i_combo,1) = all_Level_param_1_val_and_lvl_inline(p1,1);
                    combos(i_combo,2) = all_Level_param_2_val_and_lvl_inline(p2,1);
                    combos(i_combo,3) = all_Level_param_3_val_and_lvl_inline(p3,1);
                    i_combo = i_combo + 1;
                    nb_combos_found = nb_combos_found +1;
                end            
            end
        end
    end
    %resizing the table
    final_size_combos = size(combos,1);
    %all_Level_param_3_val_and_lvl_inline(nb_combos_found+1:final_size_combos,:) = [];
    combos(nb_combos_found+1:final_size_combos,:) = [];

    %PLOT THE COMBOS
    %creation of vectors for plotting and saving the DOE table
    param_1_table = zeros (nb_combos_found,1);
    param_2_table = zeros (nb_combos_found,1);
    param_3_table = zeros (nb_combos_found,1);
    DOE_to_be_written = zeros (nb_combos_found,3);


    for(i = 1 : 1 : nb_combos_found)
        param_1_table(i) = combos(i,1);
        param_2_table(i) = combos(i,2);
        param_3_table(i) = combos(i,3);
    end

    scatter3(param_1_table,param_2_table,param_3_table, 'filled');

    %WRITING LIST OF PARAMETERS COMBOS
    %for avoiding merging with longer previously generated file
    delete 'my_DOE.xlsx'

    %writing if three dimensions/parameters
    xlswrite('my_DOE.xlsx',combos);
end




%RESEARCH OF ALL COMBO FOR THE WANTED LEVEL
% research of all combinations for the wanted level
% !!! 2 DIMENSION/PARAMETER
if(xmin_param_1 ~= xmax_param_1 && xmin_param_2 ~= xmax_param_2 && xmin_param_3 == xmax_param_3)
    combos = zeros(final_size^2,3);
    i_combo = 1;
    nb_combos_found = 0;

    for(p1 = 1 : 1 : final_size)
        for(p2 = 1 : 1 : final_size) 
            %scanning all combo and compute their cumulated levels
            combo_lvl = all_Level_param_1_val_and_lvl_inline(p1,2)+ all_Level_param_2_val_and_lvl_inline(p2,2);
            %adding combo if it matches with the wanted level
            if(combo_lvl == lvlmax)
                combos(i_combo,1) = all_Level_param_1_val_and_lvl_inline(p1,1);
                combos(i_combo,2) = all_Level_param_2_val_and_lvl_inline(p2,1);
                i_combo = i_combo + 1;
                nb_combos_found = nb_combos_found +1;
            end            
        end
    end
    %resizing the table
    final_size_combos = size(combos,1);
    %all_Level_param_3_val_and_lvl_inline(nb_combos_found+1:final_size_combos,:) = [];
    combos(nb_combos_found+1:final_size_combos,:) = [];

    %PLOT THE COMBOS
    %creation of vectors for plotting and saving the DOE table
    param_1_table = zeros (nb_combos_found,1);
    param_2_table = zeros (nb_combos_found,1);
    DOE_to_be_written = zeros (nb_combos_found,2);


    for(i = 1 : 1 : nb_combos_found)
        param_1_table(i) = combos(i,1);
        param_2_table(i) = combos(i,2);
        param_3_table(i) = 0;
    end

    scatter3(param_1_table,param_2_table,param_3_table, 'filled');

    %WRITING LIST OF PARAMETERS COMBOS
    %for avoiding merging with longer previously generated file
    delete 'my_DOE.xlsx'
    %writing if two dimensions/parameters
    xlswrite('my_DOE.xlsx',combos(:,1:2));
end


%!!! 1 DIMENSION/PARAMETER : NO COMBO RESEARCH, JUST PLOT THE POINTS
%GIVEN BY THE QUADRATURES AT ALL LEVELS
    for(i = 1 : 1 : size(all_Level_param_3_val_and_lvl_inline,1))
        param_1_table(i) = all_Level_param_3_val_and_lvl_inline(i,1);
        param_2_table(i) = 0;
        param_3_table(i) = 0;
    end

    scatter3(param_1_table,param_2_table,param_3_table, 'filled');
%writing if one dimensions/parameters
if(xmin_param_1 ~= xmax_param_1 && xmin_param_2 == xmax_param_2 && xmin_param_3 == xmax_param_3)
    xlswrite('my_DOE.xlsx',all_Level_param_3_val_and_lvl_inline(:,1));
end







