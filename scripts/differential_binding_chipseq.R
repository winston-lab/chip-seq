library(tidyverse)
library(magrittr)
library(DESeq2)
library(viridis)
library(scales)
library(gridExtra)

get_countdata = function(path, samples){
    df = read_tsv(path) %>%
        select(samples) %>%
        rowid_to_column(var="index") %>%
        column_to_rownames(var="index") %>%
        as.data.frame()
    df = df[rowSums(df)>1,]
    return(df)
}

initialize_dds = function(comparison = "enrichment_vs_enrichment",
                          data_path,
                          samples,
                          conditions,
                          rna_sources,
                          condition_id,
                          control_id){
    if (comparison=="enrichment_vs_enrichment"){
        dds = DESeqDataSetFromMatrix(countData = get_countdata(data_path, samples),
                                     colData = data.frame(condition = factor(conditions,
                                                                             levels = c(control_id,
                                                                                        condition_id)),
                                                          rna_source = factor(rna_sources,
                                                                              levels = c("input",
                                                                                         "ChIP")),
                                                          row.names = samples),
                                     design = ~ rna_source + condition + rna_source:condition)
    } else {
        indices = (conditions == ifelse(comparison=="condition", condition_id, control_id))
        dds = DESeqDataSetFromMatrix(countData = get_countdata(data_path, samples[indices]),
                                     colData = data.frame(rna_source = factor(rna_sources[indices],
                                                                              levels = c("input",
                                                                                         "ChIP")),
                                                          row.names = samples[indices]),
                                     design = ~ rna_source)
    }
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
    trans_new(paste0("reverselog-", format(base)), trans, inv,
              log_breaks(base = base),
              domain = c(1e-100, Inf))
}

mean_sd_plot = function(df, ymax, title){
    ggplot(data = df, aes(x=rank, y=sd)) +
        geom_hex(aes(fill=log10(..count..), color=log10(..count..)), bins=100, size=0) +
        geom_smooth(color="#4292c6") +
        scale_fill_viridis(option="inferno", name=expression(log[10](count)), guide=FALSE) +
        scale_color_viridis(option="inferno", guide=FALSE) +
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
                                 annotations,
                                 polyenrichment,
                                 alpha,
                                 lfc){
    results(dds,
            alpha=alpha,
            lfcThreshold=lfc,
            altHypothesis="greaterAbs",
            tidy=TRUE) %>%
        as_tibble() %>%
        left_join(annotations, ., by=c("index"="row")) %>%
        left_join(polyenrichment, ., by=c("row"="index")) %>%
        arrange(padj) %>%
        mutate(name = if_else(name==".",
                              paste0("peak_", row_number()),
                              name),
               score = as.integer(pmin(-125*log2(padj), 1000))) %>%
        mutate_at(vars(pvalue, padj), ~(-log10(.))) %>%
        mutate_if(is.double, round, 3) %>%
        select(index=row, chrom, start, end, name, score, strand,
               log2FC_enrichment=log2FoldChange, lfc_SE=lfcSE,
               stat, log10_pval=pvalue, log10_padj=padj, mean_expr=baseMean,
               condition_enrichment, condition_enrichment_SE,
               control_enrichment, control_enrichment_SE) %>%
        return()
}

write_counts_table = function(results_df,
                              annotations,
                              counts_df,
                              output_path){
    results_df %>%
        select(1:7) %>%
        right_join(annotations %>% select(-c(name, score)),
                   by = c("index", "chrom", "start", "end", "strand")) %>%
        left_join(counts_df, by="index") %>%
        select(-index) %>%
        write_tsv(output_path) %>%
        return()
}

plot_ma = function(df_sig = results_df_filtered_significant,
                   df_nonsig = results_df_filtered_nonsignificant,
                   xvar = mean_expr,
                   yvar = log2_enrichment,
                   lfc,
                   condition,
                   control){
    xvar = enquo(xvar)
    yvar = enquo(yvar)
    ggplot() +
        geom_hline(yintercept = 0, color="black", linetype="dashed") +
        geom_hline(yintercept = c(-lfc, lfc), color="grey70", linetype="dashed") +
        stat_bin_hex(data = df_nonsig,
                     geom="point",
                     aes(x=!!xvar, y=!!yvar, alpha=..count..),
                     binwidth = c(.01, 0.01),
                     color="black", stroke=0, size=0.7) +
        stat_bin_hex(data = df_sig,
                     geom="point",
                     aes(x=!!xvar, y=!!yvar, alpha=..count..),
                     binwidth = c(.01, 0.01),
                     color="red", stroke=0, size=0.7) +
        scale_x_log10(name="mean of normalized counts") +
        scale_alpha_continuous(range = c(0.5, 1)) +
        ylab(bquote(log[2]~frac("enrichment in " * .(condition),
                                "enrichment in " * .(control)))) +
        theme_light() +
        theme(text = element_text(size=8, color="black"),
              axis.text = element_text(color = "black"),
              axis.title.y = element_text(angle=0, hjust=1, vjust=0.5),
              legend.position = "none")
}

plot_volcano = function(df = results_df_filtered,
                        xvar = log2_enrichment,
                        yvar = log10_padj,
                        lfc,
                        alpha,
                        condition,
                        control){
    xvar = enquo(xvar)
    yvar = enquo(yvar)
    ggplot() +
        geom_vline(xintercept = 0, color="black", linetype="dashed") +
        geom_vline(xintercept = c(-lfc, lfc), color="grey70", linetype="dashed") +
        stat_bin_hex(data = df,
                     geom = "point",
                     aes(x = !!xvar, y = !!yvar, color=log10(..count..)),
                     binwidth = c(0.01, 0.1),
                     alpha=0.8, stroke=0, size=0.7) +
        geom_hline(yintercept = -log10(alpha), color="red", linetype="dashed") +
        xlab(bquote(log[2] ~ frac("enrichment in " * .(condition),
                                  "enrichment in " * .(control)))) +
        ylab(expression(-log[10] ~ FDR)) +
        scale_color_viridis(option="inferno") +
        theme_light() +
        theme(text = element_text(size=8),
              axis.title.y = element_text(angle=0, hjust=1, vjust=0.5),
              legend.position = "none")
}

main = function(exp_table="high-affinity-ZF-v-reporter-only_allsamples-experimental-ZF-chipseq-counts-peaks.tsv.gz",
                spike_table=NULL,
                samples=read_tsv(exp_table) %>% select(-c(1:6)) %>% names(),
                conditions=rep(c(rep("reporter-only",2), rep("high-affinity-ZF",2)), 2),
                rna_sources=c(rep("input",4), rep("ChIP", 4)),
                norm="libsizenorm",
                condition="high-affinity-ZF",
                control="reporter-only",
                alpha=0.1,
                lfc=0,
                counts_norm_out="counts_norm.tsv",
                counts_rlog_out="counts_rlog.tsv",
                results_all_out="results_all.tsv",
                results_up_out="results_up.tsv",
                results_down_out="results_down.tsv",
                results_unchanged_out="results_unch.tsv",
                # bed_all_out="all.bed",
                # bed_up_out="up.bed",
                # bed_down_out="down.bed",
                # bed_unchanged_out="nonsignificant.bed",
                qc_plots_out="qcplots.png"){

    annotations = read_tsv(exp_table) %>%
        select(1:6) %>%
        rownames_to_column(var="index") %>%
        mutate(chrom = str_replace(chrom, "-minus$|-plus$", ""))

    # initialize 3 DESeq datasets,
    # for control polysome vs. total,
    # condition polysome vs. total,
    # and condition polysome enrichment vs control polysome enrichment
    #
    # the two polysome vs. total comparisons are only used to
    # get regularized enrichment estimates for the output tables
    dds = initialize_dds(comparison = "enrichment_vs_enrichment",
                         data_path = exp_table,
                         samples = samples,
                         conditions = conditions,
                         rna_sources = rna_sources,
                         condition_id = condition,
                         control_id = control)
    dds_condition = initialize_dds(comparison = "condition",
                                   data_path = exp_table,
                                   samples = samples,
                                   conditions = conditions,
                                   rna_sources = rna_sources,
                                   condition_id = condition,
                                   control_id = control)
    dds_control = initialize_dds(comparison = "control",
                                 data_path = exp_table,
                                 samples = samples,
                                 conditions = conditions,
                                 rna_sources = rna_sources,
                                 condition_id = condition,
                                 control_id = control)

    if (norm=="spikenorm"){
        dds_spike = initialize_dds(comparison = "enrichment_vs_enrichment",
                                   data_path = spike_table,
                                   samples = samples,
                                   conditions = conditions,
                                   rna_sources = rna_sources,
                                   condition_id = condition,
                                   control_id = control) %>%
            estimateSizeFactors()
        dds_spike_condition = initialize_dds(comparison = "condition",
                                             data_path = spike_table,
                                             samples = samples,
                                             conditions = conditions,
                                             rna_sources = rna_sources,
                                             condition_id = condition,
                                             control_id = control) %>%
            estimateSizeFactors()
        dds_spike_control = initialize_dds(comparison = "control",
                                           data_path = spike_table,
                                           samples = samples,
                                           conditions = conditions,
                                           rna_sources = rna_sources,
                                           condition_id = condition,
                                           control_id = control) %>%
            estimateSizeFactors()

        sizeFactors(dds) = sizeFactors(dds_spike)
        sizeFactors(dds_condition) = sizeFactors(dds_spike_condition)
        sizeFactors(dds_control) = sizeFactors(dds_spike_control)
    } else {
        dds %<>% estimateSizeFactors()
        dds_condition %<>% estimateSizeFactors()
        dds_control %<>% estimateSizeFactors()
    }
    dds %<>% estimateDispersions() %>% nbinomWaldTest()
    dds_condition %<>% estimateDispersions() %>% nbinomWaldTest()
    dds_control %<>% estimateDispersions() %>% nbinomWaldTest()

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

    polyenrichment = results(dds_condition, tidy=TRUE) %>%
        as_tibble() %>%
        select(row,
               condition_enrichment = log2FoldChange,
               condition_enrichment_SE = lfcSE) %>%
        full_join(results(dds_control, tidy=TRUE) %>%
                      as_tibble() %>%
                      select(row,
                             control_enrichment = log2FoldChange,
                             control_enrichment_SE = lfcSE),
                  by = "row")

    results_df = extract_deseq_results(dds = dds,
                                       annotations = annotations,
                                       polyenrichment = polyenrichment,
                                       alpha = alpha,
                                       lfc = lfc) %>%
        mutate(chrom = str_replace(chrom, "-minus$|-plus$", ""))

    write_counts_table(results_df = results_df,
                       annotations = annotations,
                       counts_df = counts_norm,
                       output_path = counts_norm_out)
    write_counts_table(results_df = results_df,
                       annotations = annotations,
                       counts_df = counts_rlog,
                       output_path = counts_rlog_out)

    results_df %<>%
        select(-index) %>%
        write_tsv(results_all_out)
    # results_df %>%
    #     select(1:6) %>%
    #     write_tsv(bed_all_out, col_names=FALSE)

    results_df_significant = results_df %>%
        filter(log10_padj > -log10(alpha))
    results_df_nonsignificant = results_df %>%
        filter(log10_padj <= -log10(alpha)) %>%
        write_tsv(results_unchanged_out)
    # results_df_nonsignificant %>%
    #     select(1:6) %>%
    #     write_tsv(bed_unchanged_out, col_names=FALSE)

    results_df_significant %>%
        filter(log2FC_enrichment >= 0) %>%
        write_tsv(results_up_out)
        # write_tsv(results_up_out) %>%
        # select(1:6) %>%
        # write_tsv(bed_up_out, col_names=FALSE)

    results_df_significant %>%
        filter(log2FC_enrichment < 0) %>%
        write_tsv(results_down_out)
        # write_tsv(results_down_out) %>%
        # select(1:6) %>%
        # write_tsv(bed_down_out, col_names=FALSE)

    maplot = plot_ma(df_sig = results_df_significant,
                     df_nonsig = results_df_nonsignificant,
                     xvar = mean_expr,
                     yvar = log2FC_enrichment,
                     lfc = lfc,
                     condition = condition,
                     control = control)

    volcano = plot_volcano(df = results_df,
                           xvar = log2FC_enrichment,
                           yvar = log10_padj,
                           lfc = lfc,
                           alpha = alpha,
                           condition = condition,
                           control = control)

    qc_plots = arrangeGrob(mean_sd_plot_pre,
                           mean_sd_plot_post,
                           maplot,
                           volcano,
                           ncol=2)

    ggsave(qc_plots_out,
           plot = qc_plots,
           width = 16*1.5,
           height = 9*1.5,
           units="cm")
}

main(exp_table = snakemake@input[["exp_counts"]],
     spike_table = snakemake@input[["spike_counts"]],
     samples = snakemake@params[["samples"]],
     conditions = snakemake@params[["conditions"]],
     rna_sources = snakemake@params[["sampletypes"]],
     norm = snakemake@wildcards[["norm"]],
     condition = snakemake@wildcards[["condition"]],
     control = snakemake@wildcards[["control"]],
     alpha = snakemake@params[["alpha"]],
     lfc = snakemake@params[["lfc"]],
     counts_norm_out = snakemake@output[["counts_norm"]],
     counts_rlog_out = snakemake@output[["counts_rlog"]],
     results_all_out = snakemake@output[["results_all"]],
     results_up_out = snakemake@output[["results_up"]],
     results_down_out = snakemake@output[["results_down"]],
     results_unchanged_out = snakemake@output[["results_nonsig"]],
     # bed_all_out = snakemake@output[["bed_all"]],
     # bed_up_out = snakemake@output[["bed_up"]],
     # bed_down_out = snakemake@output[["bed_down"]],
     # bed_unchanged_out = snakemake@output[["bed_nonsig"]],
     qc_plots_out = snakemake@output[["qc_plots"]])

