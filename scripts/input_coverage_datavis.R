library(tidyverse)
library(magrittr)
library(ggthemes)

main = function(
    input_path = "Spn1_chipseq_union-bedgraph-libsizenorm-midpoint-window-200-allsamples.tsv.gz",
    input_samples = c("wt-30C-input-1",
                       "wt-30C-input-2",
                       "spt6-YW-30C-input-1",
                       "spt6-YW-30C-input-2",
                       "pob3-E154K-30C-input-1",
                       "pob3-E154K-30C-input-2",
                       "spt5-QS-30C-input-1",
                       "spt5-QS-30C-input-2",
                       "spt6-YW-pob3-E154K-30C-input-1",
                       "spt6-YW-pob3-E154K-30C-input-2",
                       "spt6-YW-spt5-QS-30C-input-1",
                       "spt6-YW-spt5-QS-30C-input-2"),
    input_groups = c("wt-30C",
                     "wt-30C",
                     "spt6-YW-30C",
                     "spt6-YW-30C",
                     "pob3-E154K-30C",
                     "pob3-E154K-30C",
                     "spt5-QS-30C",
                     "spt5-QS-30C",
                     "spt6-YW-pob3-E154K-30C",
                     "spt6-YW-pob3-E154K-30C",
                     "spt6-YW-spt5-QS-30C",
                     "spt6-YW-spt5-QS-30C"),
    haploid_quantile = 0.1,
    coverage_dotplot_out = "coverage_dotplot.pdf",
    violin_facet_sample_out = "violin_facet_sample.pdf",
    violin_facet_chrom_out = "violin_facet_chrom.pdf"){

    group_mapping = tibble(sample = input_samples,
                           group = input_groups) %>%
        group_by(group) %>%
        mutate(replicate = row_number())
    max_replicates = max(group_mapping[["replicate"]])
    n_groups = n_distinct(group_mapping[["group"]])

    df = read_tsv(input_path) %>%
        pivot_longer(-1,
                     names_to="sample",
                     values_to="value") %>%
        filter(sample %in% input_samples) %>%
        left_join(group_mapping,
                  by="sample") %>%
        separate(name,
                 into=c("chrom", "start", "end"),
                 sep="-",
                 convert=TRUE) %>%
        mutate(group=fct_inorder(group, ordered=TRUE),
               chrom=fct_inorder(chrom, ordered=TRUE),
               value=value / (end-start))

    haploid_value = df %>%
        group_by(group, replicate, chrom) %>%
        summarize(median_value = median(value)) %>%
        ungroup %>%
        summarize(haploid_value = quantile(median_value, haploid_quantile)) %>%
        pull(haploid_value)

    chrom_sizes = df %>%
        group_by(chrom) %>%
        summarize(size = max(end)) %>%
        ungroup() %>%
        mutate(cum_start = cumsum(size) - size)
    n_chroms = nrow(chrom_sizes)

    df %<>%
        mutate(value=scales::rescale(value, from=c(0, haploid_value))) %>%
        left_join(chrom_sizes,
                  by="chrom") %>%
        mutate(position=(start+end)/2 + cum_start)

    coverage_dotplot = ggplot(data=df,
                              aes(x=position,
                                  y=value,
                                  color=chrom)) +
        geom_hline(yintercept=1,
                   color="grey70",
                   size=0.3,
                   alpha=0.5) +
        geom_point(size=0.1,
                   shape=1,
                   alpha=0.2) +
        # geom_rug(sides="b",
                 # size=0.01,
                 # alpha=0.05) +
        facet_grid(group~replicate,
                   switch="y") +
        scale_x_continuous(breaks=chrom_sizes[["cum_start"]] + (chrom_sizes[["size"]]/ 2),
                           labels=chrom_sizes[["chrom"]],
                           expand=c(0.01,0),
                           sec.axis=dup_axis()) +
        scale_y_log10(limits=c(0.25,4),
                           oob=scales::squish,
                      breaks=c(0.25, 0.5,1,2,4),
                      labels=c("n/4", "n/2", "n", "2n", "4n"),
                      name=NULL,
                      sec.axis=dup_axis()) +
        scale_color_manual(values=rep(viridis::viridis(2, end=0.3),
                                      ceiling(n_chroms/2)),
                           guide=FALSE) +
        theme_light() +
        theme(axis.title.x=element_blank(),
              axis.text.x.bottom=element_text(angle=30, hjust=1, vjust=1),
              axis.text.x.top=element_text(angle=30, hjust=0, vjust=0),
              panel.grid=element_blank(),
              # panel.grid.major.x=element_blank(),
              # panel.grid.minor.x=element_blank(),
              # panel.grid.major.y=element_line(),
              # panel.grid.minor.y=element_blank(),
              axis.text=element_text(color="black"),
              strip.background=element_blank(),
              strip.text=element_text(color="black"),
              strip.text.y=element_text(angle=-180, vjust=0.5, hjust=1),
              strip.placement ="outside")

    ggsave(coverage_dotplot_out,
           plot=coverage_dotplot,
           width= 16 * max_replicates,
           height= 9 * n_groups / 2,
           units="cm",
           limitsize=FALSE)

    violin_facet_sample = ggplot(data=df,
           aes(x=chrom,
               y=value,
               fill=chrom)) +
        geom_hline(yintercept=1,
                   size=0.5,
                   color="grey70") +
        geom_violin(draw_quantiles=c(0.5),
                    scale="width") +
        scale_y_log10(limits=c(0.2, 5),
                      breaks=c(0.25, 0.5, 1, 2, 4),
                      labels=c("n/4", "n/2", "n", "2n", "4n"),
                      sec.axis=dup_axis()) +
        facet_grid(group~replicate,
                   switch="y") +
        scale_fill_viridis_d(begin=0.3, guide=FALSE) +
        theme_light() +
        theme(legend.position="none",
              strip.placement="outside",
              strip.background=element_blank(),
              strip.text=element_text(color="black"),
              strip.text.y=element_text(angle=-180, vjust=0.5, hjust=1),
              axis.text=element_text(color="black"),
              axis.text.x=element_text(angle=30, hjust=1, vjust=1),
              axis.text.x.top=element_text(color="black"),
              axis.title=element_blank(),
              panel.grid.minor=element_blank(),
              panel.grid.major.x=element_blank(),
              panel.grid.major.y=element_line(size=0.4, linetype="dotted"))

    ggsave(violin_facet_sample_out,
           plot=violin_facet_sample,
           width=16 * max_replicates,
           height=9 * n_groups / 3,
           units="cm",
           limitsize=FALSE)

    violin_facet_chrom = ggplot(data=df,
           aes(x=interaction(replicate, group),
               y=value,
               fill=group)) +
        geom_hline(yintercept=1,
                   size=0.5,
                   color="grey70") +
        geom_violin(size=0.2,
                    draw_quantiles=c(0.5),
                    position=position_dodge(width=1)) +
        scale_y_log10(limits=c(0.2, 5),
                      breaks=c(0.25, 0.5, 1, 2, 4),
                      labels=c("n/4", "n/2", "n", "2n", "4n"),
                      sec.axis=dup_axis()) +
        facet_wrap(~chrom) +
        scale_fill_manual(values=rep(ptol_pal()(min(n_groups, 12)), ceiling(n_groups/12))) +
        theme_light() +
        theme_light() +
        theme(legend.position="bottom",
              strip.background=element_blank(),
              strip.text=element_text(color="black",
                                      hjust=0),
              axis.ticks.x=element_blank(),
              axis.text=element_text(color="black"),
              axis.text.x=element_blank(),
              axis.title=element_blank(),
              panel.grid.minor=element_blank(),
              panel.grid.major.x=element_blank(),
              panel.grid.major.y=element_line(size=0.4, linetype="dotted"))

    ggsave(violin_facet_chrom_out,
           plot=violin_facet_chrom,
           width=16*1.5,
           height=9*1.5,
           units="cm",
           limitsize=FALSE)
}

main(input_path = snakemake@input[["tsv"]],
     input_samples = snakemake@params[["input_samples"]],
     input_groups = snakemake@params[["input_groups"]],
     haploid_quantile = snakemake@params[["haploid_quantile"]],
     coverage_dotplot_out = snakemake@output[["coverage_dotplot"]],
     violin_facet_sample_out = snakemake@output[["violin_facet_sample"]],
     violin_facet_chrom_out = snakemake@output[["violin_facet_chrom"]])
