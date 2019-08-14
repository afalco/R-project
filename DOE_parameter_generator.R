GLC_quadrature <- function( xmin, xmax, wanted_level)
{
# Generates an array of point corresponding to the Gauss-Lobatto-Tchebishev
# quadrature comprised between user parameters xmin and xmax. The wanted
# level for the quadrature is set by the wanted_level argument
    x <- rep(0,2^(wanted_level-1))
    Points <- 2^wanted_level
    index <- seq(from=2, to=Points,by=2)
    aux <- -cos(pi*(index-1)/Points)
    x<- xmin + (xmax-xmin)*(aux+1)/2
    return(x)
}
#we use the library('scatterplot3d') and library('WriteXLS')
library('scatterplot3d')
library('WriteXLS')

#value of the maximal level for Gauss Lobatto Chebyshev quadratures
lvlmax <- 9
#min/max values for parameter 1 
xmin_param_1 <- 0
xmax_param_1 <- 10
#min/max values for parameter 2
xmin_param_2 <- 10
xmax_param_2 <-20
#min/max values for parameter 3
xmin_param_3 <- -10
xmax_param_3 <- 5
# !!! TO CLOSE ONE DIMENSION/PARAMETER, JUST PUT THE SAME VALUE FOR XMIN AND
#XMAX !!!


#arrays containing all level of the quadratures for each parameter
all_Level_quadrature_param_1 <-  matrix(0,nrow=10,ncol=10000)
all_Level_quadrature_param_2 <-  matrix(0,nrow=10,ncol=10000)
all_Level_quadrature_param_3 <-  matrix(0,nrow=10,ncol=10000)
#array containing the size of each level of the Gauss Lobatto Chebyshev quadratures
level_size_GLC_quadrature <-  matrix(0,nrow=50,ncol=1)

#arrays containing all level of the quadratures values for each parameter INLINE
#with their associated level
all_Level_param_1_val_and_lvl_inline <- matrix(0,nrow=10000,ncol=2)
all_Level_param_2_val_and_lvl_inline <- matrix(0,nrow=10000,ncol=2)
all_Level_param_3_val_and_lvl_inline <- matrix(0,nrow=10000,ncol=2)

#GENERATING ALL LEVEL QUADRATURES FOR ALL PARAMETERS

#Gauss Lobatto Chebyshev quadrature for parameter 1

#two first values for level 0
all_Level_param_1_val_and_lvl_inline[1,1] <- xmin_param_1
all_Level_param_1_val_and_lvl_inline[1,2] <- 0
all_Level_param_1_val_and_lvl_inline[2,1] <- xmax_param_1
all_Level_param_1_val_and_lvl_inline[2,2] <- 0
compt <- 3

for (lvl in 1:lvlmax){
    current_quad <- GLC_quadrature(xmin_param_1, xmax_param_1,lvl)
    size_of_current_quad <- length(current_quad)
    level_size_GLC_quadrature[lvl] <- size_of_current_quad
    all_Level_quadrature_param_1[lvl,1:size_of_current_quad] <- current_quad   
    #inline storage of the parameter values
    transp_current_quad <- t(current_quad)
    all_Level_param_1_val_and_lvl_inline[compt:(compt+size_of_current_quad-1),1] <- transp_current_quad
    #and their associated levels
    all_Level_param_1_val_and_lvl_inline[compt:(compt+size_of_current_quad-1),2] <- lvl  
    compt <- compt + size_of_current_quad
}
#resizing the table
final_size <- compt-1
#all_Level_param_1_val_and_lvl_inline[final_size+1:10000,] <- []


#Gauss Lobatto Chebyshev quadrature for parameter 2

#two first values for level 0
all_Level_param_2_val_and_lvl_inline[1,1] <- xmin_param_2
all_Level_param_2_val_and_lvl_inline[1,2] <- 0
all_Level_param_2_val_and_lvl_inline[2,1] <- xmax_param_2
all_Level_param_2_val_and_lvl_inline[2,2] <- 0
compt <- 3

for (lvl in 1:lvlmax){
    current_quad <- GLC_quadrature( xmin_param_2, xmax_param_2, lvl )
    size_of_current_quad <- length(current_quad)
    all_Level_quadrature_param_2[lvl, 1:size_of_current_quad]<- current_quad
    #inline storage of the parameter values
    transp_current_quad <- t(current_quad)
    all_Level_param_2_val_and_lvl_inline[compt:(compt+size_of_current_quad-1),1] <- transp_current_quad
    #and their associated levels
    all_Level_param_2_val_and_lvl_inline[compt:(compt+size_of_current_quad-1),2]<-lvl
    compt <- compt + size_of_current_quad
}
#resizing the table
final_size <- compt-1
#all_Level_param_2_val_and_lvl_inline(final_size+1:10000,:) <- []


#Gauss Lobatto Chebyshev quadrature for parameter 3

#two first values for level 0
all_Level_param_3_val_and_lvl_inline[1,1] <- xmin_param_3
all_Level_param_3_val_and_lvl_inline[1,2] <- 0
all_Level_param_3_val_and_lvl_inline[2,1] <- xmax_param_3
all_Level_param_3_val_and_lvl_inline[2,2] <- 0
compt <- 3

for (lvl in 1:lvlmax){
    current_quad <- GLC_quadrature(xmin_param_3, xmax_param_3, lvl )
    size_of_current_quad <- length(current_quad)
    all_Level_quadrature_param_3[lvl, 1:size_of_current_quad] <- current_quad  
    #inline storage of the parameter values
    transp_current_quad <- t(current_quad)
    all_Level_param_3_val_and_lvl_inline[compt:(compt+size_of_current_quad-1),1] <- transp_current_quad
    #and their associated levels
    all_Level_param_3_val_and_lvl_inline[compt:(compt+size_of_current_quad-1),2]<-lvl
    compt <- compt + size_of_current_quad
}
#resizing the table
final_size <- compt-1
#all_Level_param_3_val_and_lvl_inline(final_size+1:10000,:) <- []



#RESEARCH OF ALL COMBO FOR THE WANTED LEVEL
# research of all combinations for the wanted level
# !!! 3 DIMENSION/PARAMETER
if (xmin_param_1 != xmax_param_1 && xmin_param_2 != xmax_param_2 && xmin_param_3 != xmax_param_3){
    combos <- matrix(0,final_size^2,3)
    i_combo <- 1
    nb_combos_found <- 0

    for (p1 in 1:final_size){
        for (p2 in 1:final_size){
            for(p3 in 1:final_size){
                #scanning all combo and compute their cumulated levels
                combo_lvl <- all_Level_param_1_val_and_lvl_inline[p1,2]+all_Level_param_2_val_and_lvl_inline[p2,2]+all_Level_param_3_val_and_lvl_inline[p3,2]
                #adding combo if it matches with the wanted level
                    if(combo_lvl ==lvlmax){
                    combos[i_combo,1] <- all_Level_param_1_val_and_lvl_inline[p1,1]
                    combos[i_combo,2] <- all_Level_param_2_val_and_lvl_inline[p2,1]
                    combos[i_combo,3] <- all_Level_param_3_val_and_lvl_inline[p3,1]
                    i_combo <- i_combo + 1
                    nb_combos_found <- nb_combos_found +1
                    }
                }            
            }
        }
    #resizing the table
    final_size_combos <- length(combos)[1]
    #all_Level_param_3_val_and_lvl_inline(nb_combos_found+1:final_size_combos,:) <- []
    #combos((nb_combos_found+1):final_size_combos,:) <- []

    #PLOT THE COMBOS
    #creation of vectors for plotting and saving the DOE table
    param_1_table <- rep(0,nb_combos_found)
    param_2_table <- rep(0,nb_combos_found)
    param_3_table <- rep(0,nb_combos_found)
    DOE_to_be_written <- matrix(0,nb_combos_found,3)


    for(i in 1:nb_combos_found){
        param_1_table[i] <- combos[i,1]
        param_2_table[i] <- combos[i,2]
        param_3_table[i] <- combos[i,3]
    }
    scatterplot3d(param_1_table,param_2_table,param_3_table)

    #WRITING LIST OF PARAMETERS COMBOS
    #for avoiding merging with longer previously generated file
    #delete 'my_DOE.xlsx'

    #writing if three dimensions/parameters
    #write.csv(combos,'./my_DOE.csv')
    WriteXLS(data.frame(combos),'./my_DOE.xlsx')
}


#RESEARCH OF ALL COMBO FOR THE WANTED LEVEL
# research of all combinations for the wanted level
# !!! 2 DIMENSION/PARAMETER
if(xmin_param_1 != xmax_param_1 && xmin_param_2 != xmax_param_2 && xmin_param_3 == xmax_param_3){
    combos <- rep(0,final_size^2,3)
    i_combo <- 1
    nb_combos_found <- 0

    for(p1 in 1:final_size){
        for(p2 in 1:final_size){
            #scanning all combo and compute their cumulated levels
            combo_lvl <- all_Level_param_1_val_and_lvl_inline[p1,2]+ all_Level_param_2_val_and_lvl_inline[p2,2]
            #adding combo if it matches with the wanted level
            if(combo_lvl==lvlmax){
                combos[i_combo,1] <- all_Level_param_1_val_and_lvl_inline[p1,1]
                combos[i_combo,2] <- all_Level_param_2_val_and_lvl_inline[p2,1]
                i_combo <- i_combo + 1
                nb_combos_found <- nb_combos_found +1
            }            
        }
    }
    #resizing the table
    final_size_combos <- length(combos)
    #all_Level_param_3_val_and_lvl_inline(nb_combos_found+1:final_size_combos,:) <- []
    #combos(nb_combos_found+1:final_size_combos,:) <- []

    #PLOT THE COMBOS
    #creation of vectors for plotting and saving the DOE table
    param_1_table <- rep(0,nb_combos_found)
    param_2_table <- rep(0,nb_combos_found)
    DOE_to_be_written <- rep(0,nb_combos_found,2)


    for(i in 1: nb_combos_found)
        param_1_table(i) <- combos(i,1)
        param_2_table(i) <- combos(i,2)
        param_3_table(i) <- 0
    end
    scatterplot3d(param_1_table,param_2_table,param_3_table)

    #WRITING LIST OF PARAMETERS COMBOS
    #for avoiding merging with longer previously generated file
    #delete 'my_DOE.xlsx'
    #writing if two dimensions/parameters
    #write.csv(combos[,1:2],'./my_DOE.csv')
    WriteXLS(data.frame(combos[,1:2]),'./my_DOE.xlsx')
}


#!!! 1 DIMENSION/PARAMETER : NO COMBO RESEARCH, JUST PLOT THE POINTS
#GIVEN BY THE QUADRATURES AT ALL LEVELS
    for(i in 1:length(all_Level_param_3_val_and_lvl_inline)){
        param_1_table[i] <- all_Level_param_3_val_and_lvl_inline[i]
        param_2_table[i] <- 0
        param_3_table[i] <- 0
    }

scatterplot3d(param_1_table,param_2_table,param_3_table)
#writing if one dimensions/parameters
if(xmin_param_1 != xmax_param_1 && xmin_param_2 == xmax_param_2 && xmin_param_3 == xmax_param_3){
    #write.csv(all_Level_param_3_val_and_lvl_inline,'./my_DOE.csv')
    WriteXLS(data.frame(all_Level_param_3_val_and_lvl_inline),'./my_DOE.xlsx')
}
