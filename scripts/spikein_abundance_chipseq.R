library(tidyverse)
library(forcats)
library(gridExtra)
library(ggthemes)

main = function(in_path, sample_list, controls, conditions, plot_out, stats_out){
    df = read_tsv(in_path) %>%
        filter(sample %in% sample_list) %>%
        mutate(abundance = (experimental_counts_IP / spikein_counts_IP ) *
                                    (spikein_counts_input / experimental_counts_input)) %>%
        mutate_at(vars(sample, group), funs(fct_inorder(., ordered=TRUE))) %>%
        group_by(group) %>%
        mutate(outlier= ifelse(abundance >
                                   2.5*quantile(abundance, .75) -
                                   1.5*quantile(abundance, .25) |
                                   abundance <
                                   -2.5*quantile(abundance, .25) -
                                   1.5*quantile(abundance, .75),
                               TRUE, FALSE))

    n_samples = nrow(df)
    n_groups = df %>% pull(group) %>% n_distinct()

    barplot = ggplot(data=df, aes(x=sample, fill=group, y=abundance)) +
        geom_col() +
        geom_text(aes(label=round(abundance, 2)), size=12/75*25.4,
                  position=position_stack(vjust=0.9)) +
        scale_fill_ptol(guide=FALSE) +
        ylab("spike-in normalized\nabundance vs. input") +
        theme_light() +
        theme(axis.text = element_text(size=10, color="black"),
              axis.text.x = element_text(angle=30, hjust=0.9),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size=10, color="black",
                                          angle=0, vjust=0.5, hjust=1))

    boxplot = ggplot(data = df, aes(x=group, y=abundance, fill=group)) +
        geom_boxplot(outlier.shape=16, outlier.size=1.5, outlier.color="red", outlier.stroke=0) +
        geom_point(shape=16, size=1, stroke=0) +
        scale_fill_ptol(guide=FALSE) +
        scale_y_continuous(name = "spike-in normalized\nabundance vs. input",
                           limits = c(0, NA)) +
        theme_light() +
        theme(axis.text = element_text(size=10, color="black"),
              axis.text.x = element_text(angle=30, hjust=0.9),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size=10, color="black"))

    stats_table = df %>%
        add_count(group) %>%
        group_by(group) %>%
        mutate(median = median(abundance)) %>%
        ungroup() %>%
        filter(!outlier) %>%
        add_count(group) %>%
        group_by(group) %>%
        summarise(n = first(n),
                  median = first(median),
                  n_no_outlier = first(nn),
                  mean_no_outlier = mean(abundance),
                  sd_no_outlier = sd(abundance)) %>%
        write_tsv(path = stats_out, col_names=TRUE)

    #set width
    wl = 1+1.6*n_samples
    wr = 1+1.8*n_groups
    th = 0
    if (!(is.null(conditions) || is.null(controls))){
        levels_df = tibble(condition=conditions, control=controls) %>%
            left_join(stats_table %>% select(group, mean_no_outlier),
                      by=c("condition"="group")) %>%
            rename(condition_abundance=mean_no_outlier) %>%
            left_join(stats_table %>% select(group, mean_no_outlier),
                      by=c("control"="group")) %>%
            rename(control_abundance=mean_no_outlier) %>%
            mutate(levels = condition_abundance/control_abundance)

        levels_table = levels_df %>%
            select(condition, control, levels) %>%
            mutate_at("levels", funs(round(., digits=3)))
        levels_draw = tableGrob(levels_table,
                                rows=NULL,
                                cols=c("condition","control","relative levels"),
                                ttheme_minimal(base_size=10))

        th = 1+length(conditions)/2
        page = arrangeGrob(barplot, boxplot, levels_draw,
                           layout_matrix=rbind(c(1,2),c(3,3)),
                           widths=unit(c(wl, wr), "cm"),
                           heights=unit(c(9,th),"cm"))
    } else {
        page = arrangeGrob(barplot, boxplot,
                           widths=unit(c(wl, wr), "cm"),
                           heights=unit(c(9,th),"cm"))

    }
    ggsave(plot_out, page, width = wl+wr, height=9+th+.5, units = "cm")
}

main(in_path = snakemake@input[[1]],
     sample_list = snakemake@params[["samplelist"]],
     controls = snakemake@params[["controls"]],
     conditions = snakemake@params[["conditions"]],
     plot_out = snakemake@output[["plot"]],
     stats_out = snakemake@output[["stats"]])

