library(tidyverse)
library(magrittr)
library(DESeq2)
# library(ashr)
library(gridExtra)

get_countdata = function(path="Spn1_chipseq_non-depleted-counts-midpoints-window-10.tsv.gz",
                         samples){
    df = read_tsv(path) %>%
        select(samples) %>%
        rowid_to_column(var="index") %>%
        column_to_rownames(var="index") %>%
        as.data.frame()
    df = df[rowSums(df)>1,]
    return(df)
}

initialize_dds = function(data_path,
                          samples,
                          rna_sources){
    dds = DESeqDataSetFromMatrix(countData = get_countdata(path=data_path,
                                                           samples=samples),
                                 colData = data.frame(rna_source = factor(rna_sources,
                                                                          levels = c("input",
                                                                                     "ChIP")),
                                                      row.names = samples),
                                 design = ~ rna_source)
    return(dds)
}

extract_normalized_counts = function(dds){
    dds %>%
        counts(normalized=TRUE) %>%
        as.data.frame() %>%
        rownames_to_column(var="index") %>%
        as_tibble() %>%
        return()
}

extract_rlog_counts = function(dds){
    dds %>%
        rlog(blind=FALSE) %>%
        assay() %>%
        as.data.frame() %>%
        rownames_to_column(var="index") %>%
        as_tibble() %>%
        return()
}

build_mean_sd_df_pre = function(dds){
     dds %>%
        normTransform() %>%
        assay() %>%
        as_tibble() %>%
        rowid_to_column(var="index") %>%
        gather(sample, signal, -index) %>%
        group_by(index) %>%
        summarise(mean = mean(signal),
                  sd = sd(signal)) %>%
        mutate(rank = min_rank(dplyr::desc(mean))) %>%
        return()
}

build_mean_sd_df_post = function(counts){
    counts %>%
        gather(sample, signal, -index) %>%
        group_by(index) %>%
        summarise(mean = mean(signal),
                  sd = sd(signal)) %>%
        mutate(rank = min_rank(dplyr::desc(mean))) %>%
        return()
}

reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    scales::trans_new(paste0("reverselog-", format(base)), trans, inv,
              scales::log_breaks(base = base),
              domain = c(1e-100, Inf))
}

mean_sd_plot = function(df, ymax, title){
    ggplot(data = df,
           aes(x=rank,
               y=sd)) +
        geom_hex(aes(fill=..count..,
                     color=..count..),
                 bins=100,
                 size=0) +
        geom_smooth(color="#4292c6") +
        scale_fill_viridis_c(option="inferno",
                             name=expression(log[10](count)),
                             guide=FALSE) +
        scale_color_viridis_c(option="inferno",
                              guide=FALSE) +
        scale_x_continuous(trans = reverselog_trans(10),
                           name="rank(mean enrichment)",
                           expand = c(0,0)) +
        scale_y_continuous(limits = c(NA, ymax),
                           name = "SD") +
        theme_light() +
        ggtitle(title) +
        theme(text = element_text(size=8))
}

extract_deseq_results = function(dds,
                                 annotations){
    lfc_shrunk = lfcShrink(dds,
                           coef="rna_source_ChIP_vs_input",
                           type="ashr") %>%
        as_tibble(rownames="index")

    results(dds,
            tidy=TRUE,
            cooksCutoff=FALSE,
            independentFiltering=FALSE) %>%
        as_tibble() %>%
        left_join(annotations, .,
                  by=c("index"="row")) %>%
        left_join(lfc_shrunk,
                  by=c("index", "baseMean"),
                  suffix=c("_original", "_shrunken")) %>%
        return()
}

write_counts_table = function(annotations,
                              counts_df,
                              output_path){
    annotations %>%
        left_join(counts_df,
                  by="index") %>%
        select(-index) %>%
        write_tsv(output_path) %>%
        return()
}

main = function(exp_table="Spn1_chipseq_non-depleted-counts-midpoints-window-20.tsv.gz",
                spike_table="Spn1_chipseq_non-depleted-sicounts-midpoints-window-20.tsv.gz",
                samples=read_tsv(exp_table) %>% select(-1) %>% names(),
                rna_sources=c(rep("input",2), rep("ChIP", 2)),
                norm="spikenorm",
                counts_norm_out="counts_norm.tsv.gz",
                counts_rlog_out="counts_rlog.tsv.gz",
                results_all_out="results_all.tsv.gz",
                bedgraph_out = "test.bedgraph",
                qc_plots_out="qcplots.png"){

    annotations = read_tsv(exp_table) %>%
        select(name) %>%
        rownames_to_column(var="index") %>%
        separate(name,
                 into=c("chrom", "start", "end"),
                 sep="-",
                 convert=TRUE)

    dds = initialize_dds(data_path=exp_table,
                         samples,
                         rna_sources)

    if (norm=="spikenorm"){
        dds_spike = initialize_dds(data_path=spike_table,
                                   samples,
                                   rna_sources) %>%
            estimateSizeFactors()
        sizeFactors(dds) = sizeFactors(dds_spike)
    } else {
        dds %<>% estimateSizeFactors()
    }
    dds %<>% estimateDispersions() %>% nbinomWaldTest()

    #extract normalized counts and write to file
    counts_norm = extract_normalized_counts(dds = dds)
    counts_rlog = extract_rlog_counts(dds = dds)

    mean_sd_df_pre = build_mean_sd_df_pre(dds)
    mean_sd_df_post = build_mean_sd_df_post(counts_rlog)

    sd_max = max(c(mean_sd_df_pre[["sd"]],
                   mean_sd_df_post[["sd"]]),
                 na.rm=TRUE)*1.01

    mean_sd_plot_pre = mean_sd_plot(df = mean_sd_df_pre,
                                    ymax = sd_max,
                                    title = expression(log[2] ~ "counts," ~ "pre-shrinkage"))
    mean_sd_plot_post = mean_sd_plot(df = mean_sd_df_post,
                                     ymax = sd_max,
                                     title = expression(regularized ~ log[2] ~ "counts"))

    results_df = extract_deseq_results(dds = dds,
                                       annotations = annotations)

    write_counts_table(annotations = annotations,
                       counts_df = counts_norm,
                       output_path = counts_norm_out)
    write_counts_table(annotations = annotations,
                       counts_df = counts_rlog,
                       output_path = counts_rlog_out)

    results_df %>%
        select(-index) %>%
        write_tsv(results_all_out) %>%
        # select(chrom, start, end, log2FoldChange_shrunken) %>%
        # replace_na(list("log2FoldChange_shrunken"=0)) %>%
        select(chrom, start, end, log2FoldChange_original) %>%
        replace_na(list("log2FoldChange_original"=0)) %>%
        arrange(chrom, start) %>%
        write_tsv(bedgraph_out,
                  col_names=FALSE)

    # shrinkage_plot = ggplot(data=results_df,
    #        aes(x=log2FoldChange_original,
    #            y=log2FoldChange_shrunken,
    #            color=log10(baseMean))) +
    #     geom_hline(yintercept = 0,
    #                size=0.2,
    #                color="gray70") +
    #     geom_vline(xintercept = 0,
    #                size=0.2,
    #                color="gray70") +
    #     geom_abline(slope=1,
    #                 intercept=0,
    #                 size=0.2,
    #                 color="gray70") +
    #     geom_point(alpha=0.5,
    #                shape=16,
    #                size=0.5) +
    #     scale_color_viridis_c(name=expression("log"[10]("mean counts"))) +
    #     scale_x_continuous(name=bquote("log"[2] ~
    #                                        textstyle(frac(.(nfactor), .(dfactor))) ~ ", original"),
    #                        breaks=scales::pretty_breaks(5)) +
    #     scale_y_continuous(name=bquote("log"[2] ~
    #                                        textstyle(frac(.(nfactor), .(dfactor))) ~ ", shrunken"),
    #                        breaks=scales::pretty_breaks(5)) +
    #     theme_light() +
    #     theme(panel.grid=element_blank(),
    #           axis.text=element_text(color="black"),
    #           legend.position=c(0.01,0.99),
    #           legend.background = element_blank(),
    #           legend.justification=c(0,1),
    #           axis.title.y=element_text(angle=0,
    #                                     vjust=0.5,
    #                                     hjust=1))

    qc_plots = arrangeGrob(mean_sd_plot_pre,
                           mean_sd_plot_post,
                           # shrinkage_plot,
                           grid::nullGrob(),
                           layout_matrix=rbind(c(1,3),
                                               c(2,3)),
                           widths=c(0.5,1))

    ggsave(qc_plots_out,
           plot = qc_plots,
           width = 16*1.5,
           height = 9*1.5,
           units="cm")
}

main(exp_table = snakemake@input[["exp_table"]],
     spike_table = snakemake@input[["spike_table"]],
     samples = snakemake@params[["samples"]],
     rna_sources = snakemake@params[["rna_sources"]],
     norm = snakemake@wildcards[["norm"]],
     counts_norm_out = snakemake@output[["counts_norm"]],
     counts_rlog_out = snakemake@output[["counts_rlog"]],
     results_all_out = snakemake@output[["tsv"]],
     bedgraph_out = snakemake@output[["bedgraph"]],
     qc_plots_out = snakemake@output[["qc_plots"]])
