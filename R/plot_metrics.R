library(tidyverse)
library(ggplot2)
library(yaml)
library(cowplot)
library(magrittr)
library(optparse)
options(stringsAsFactors = F)
options(bitmaptype="cairo" )

scale_y_origin <- function(...) {
	scale_y_continuous(expand = expansion(mult = c(0, 0.05)), ...)
}

panpipes_theme <- function(font_size=10) {
	theme(
		text = element_text(size=font_size),
		
		panel.border = element_rect(colour = "black", fill=NA),
		# color background 2)
		panel.background = element_rect(fill = "white"),
		# modify grid 3)
		# panel.grid.major.x = element_line(colour = "black", linetype = 3, size = 0.5),
		# panel.grid.minor.x = element_blank(),
		panel.grid.major =  element_line(colour = "grey", linetype = 3, size = 0.5),
		# panel.grid.minor.y = element_blank(),
		# modify text, axis and colour 4) and 5)	
		axis.text = element_text(size=font_size),
		axis.title = element_text(size=font_size),
		# axis.ticks = element_line(colour = "steelblue"),
		plot.title = element_text(hjust = 0.5),
		# legend at the bottom 6)
		legend.position = "bottom",
		
	)
}

parse_cell_metadata <- function(cmtd, cat_vars, grp_vars){
	colnames(cmtd) <- gsub(":", "_", colnames(cmtd))
	if("rep_has_ir" %in% colnames(cmtd)){
		cmtd$rep_has_ir = cmtd$rep_has_ir %>% replace_na(FALSE)
	}
	cmtd %<>% mutate(across(cat_vars$variable, factor))
	cmtd %<>% mutate(across(grp_vars$variable, factor))
	return(cmtd)
}

parse_vars <- function(cat_vars){
    if (is.list(cat_vars)){
        uniq_cat_vars <- unique(unlist(cat_vars))
        uniq_cat_vars
        # cat vars
        cat_split = data.frame(str_split(uniq_cat_vars, ":", simplify = T)) %>%
            mutate(mod =ifelse(X1 %in% c('rna', 'prot', 'atac', 'rep'), X1, "multimodal"),
                        variable = gsub(":", "_", uniq_cat_vars)) %>% 
            select(mod, variable)
    }else{
        cat_split = data.frame(variable=gsub(":", "_", cat_vars))
    }
	return(cat_split)
}

# Get options ----------------------------------------------------------------
option_list <- list(
	make_option(c("--mtd_object"),help="tsv file containing per cell metrics for plotting"),

	make_option(c("--params_yaml"), default="pipeline.yml",help="yml file in the format expected by panpipes vis")
	
)

opt <- parse_args(OptionParser(option_list=option_list))
message("Running with options:")
print(opt)

# load_data ----------------------------------------------------------------

cmtd = data.table::fread(opt$mtd_object, sep='\t', na.strings = "")

PARAMS = yaml::read_yaml(opt$params_yaml)

grp_vars =  parse_vars(PARAMS$grouping_vars)

cont_vars = parse_vars(PARAMS$continuous_vars)
cat_vars = parse_vars(PARAMS$categorical_vars)

cmtd <- parse_cell_metadata(cmtd, cat_vars, grp_vars)



# Barplots ----------------------------------------------------------------
if( PARAMS$do_plots$categorical_barplots){
    for( mod in unique(cat_vars$mod)){
        if (!dir.exists(file.path(mod))) dir.create(file.path(mod))
        print(mod)
        uniq_cat_vars = cat_vars[cat_vars$mod == mod, 'variable']
        
        for (cat_v in uniq_cat_vars){
            print(cat_v)
            plt_title <- cat_v
            plot_dat = cmtd %>% 
                mutate(!! cat_v  := factor(!! rlang::sym(cat_v))) %>% 
                group_by(!!sym(cat_v)) %>% 
                summarise(n_cells=n()) %>% 
                drop_na()
            
            # coord flip if samples have long names
            flip_axes <- max(nchar(as.character(plot_dat[[cat_v]] ))) > 10

            p1 <- ggplot(plot_dat,  aes_string(x=cat_v, y="n_cells", fill=cat_v )) +
                geom_col() +
                ggtitle(cat_v) +
                scale_y_origin() + 
                panpipes_theme() + 
                theme(legend.position = "none")
            if (flip_axes){
                p1 <- p1 + coord_flip()
            }
        
        print(paste0("saving to: ", file.path(mod, paste0("bar_", cat_v, ".png"))))
        save_plot(filename=file.path(mod, paste0("bar_", cat_v, ".png")),
                            plot=p1, ncol=1, 
                            nrow=1,base_asp = 1.3, base_height = 5) 
    }
}
}
# Stacked barplots by grouping var ----------------------------------------------------------------
if( PARAMS$do_plots$categorical_stacked_barplots){

    for( mod in unique(cat_vars$mod)){
        print(mod)
        uniq_cat_vars = cat_vars[cat_vars$mod == mod, 'variable']
        
        for( gv in grp_vars$variable){
            print(gv)
            if (!dir.exists(file.path(mod, gv))) dir.create(file.path(mod, gv), recursive=T)
            stacked_barplots <- list()
        
            for (cat_v in uniq_cat_vars){
                
                if( gv == cat_v){
                    print("skipping")
                    next
                }
                
                plot_dat = cmtd %>% select(!!sym(gv), !!sym(cat_v)) %>% 
                    drop_na() %>% 
                    group_by( !!sym(gv), !!sym(cat_v)) %>% 
                    summarise(n_cells=n()) %>% ungroup() %>% group_by(!!sym(gv)) %>% 
                    mutate(proportion= n_cells / sum(n_cells)) 
                    # drop_na()
                # coord flip if samples have long names
                flip_axes <- max(nchar(as.character(plot_dat[[gv]] ))) > 10
                
                p1 <- ggplot(plot_dat,  aes_string(x=gv, y="n_cells", fill=cat_v)) + 
                    geom_col() + 
                    # ggtitle(plt_title) + 
                    scale_y_origin() + 
                    panpipes_theme(font_size = 10) 

                p2 <- ggplot(plot_dat,  aes_string(x=gv, y="proportion", fill=cat_v)) + 
                    geom_col() + 
                    # ggtitle(plt_title) + 
                    scale_y_origin() + 
                    panpipes_theme(font_size = 10)
                
                if (flip_axes){
                    p1 <- p1 + coord_flip()
                    p2 <- p2 + coord_flip()
                }
                pg <- cowplot::plot_grid(p1, p2)
                save_plot(filename=file.path(mod, gv, paste0("stackedbar_", cat_v, ".png")),
                                    plot=pg, ncol=2, 
                                    nrow=1, base_height = 5, 
                                    base_asp = ifelse(flip_axes, 1.5, 1.1)) 
            }
        }
        
        
    }
}

# Violin by grouping var ----------------------------------------------------------------
if( PARAMS$do_plots$continuous_violin){


    for( mod in unique(cont_vars$mod)){
        print(mod)
        uniq_cont_vars = cont_vars[cont_vars$mod == mod, 'variable']
        
        for( gv in grp_vars$variable){
            print(gv)
            if (!dir.exists(file.path(mod, gv))) dir.create(file.path(mod, gv), recursive=T)

            stacked_barplots <- list()
            
            for (cont_v in uniq_cont_vars){
                
                if( gv == cont_v) next
                
                # plot_dat = cmtd %>% group_by( !!sym(gv), !!sym(cont_v)) %>% 
                #     summarise(n_cells=n()) %>%
                #     mutate(proportion= n_cells / sum(n_cells)) %>% 
                #     drop_na()
                # coord flip if samples have long names
                plot_dat = cmtd %>% filter(!is.na(!!sym(cont_v)) & !is.na(!!sym(gv)))
                flip_axes <- max(nchar(as.character(plot_dat[[gv]] ))) > 10
                
                p1 <- ggplot(plot_dat,  aes_string(x=gv, y=cont_v, fill=gv)) + 
                    geom_violin() + 
                    # ggtitle(plt_title) + 
                    scale_y_origin() + 
                    panpipes_theme(font_size = 10) 
                
                
                if (flip_axes){
                    p1 <- p1 + coord_flip()
                }

                save_plot(filename=file.path(mod, gv, paste0("violin_", cont_v, ".png")),
                                    plot=p1, ncol=1, 
                                    nrow=1, base_height = 5, 
                                    base_asp = ifelse(flip_axes, 1.5, 1.1)) 
            }
        }
        
        
    }

}
