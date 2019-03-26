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

initialize_dds = function(data_path,
                          samples,
                          groups,
                          condition_id,
                          control_id){
    DESeqDataSetFromMatrix(countData = get_countdata(data_path, samples),
                           colData = data.frame(condition = factor(groups,
                                                                   levels = c(control_id,
                                                                              condition_id)),
                                                row.names = samples),
                           design = ~ condition) %>%
        return()
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
        mutate(rank = min_rank(desc(mean))) %>%
        return()
}

build_mean_sd_df_post = function(counts){
    counts %>%
        gather(sample, signal, -index) %>%
        group_by(index) %>%
        summarise(mean = mean(signal),
                  sd = sd(signal)) %>%
        mutate(rank = min_rank(desc(mean))) %>%
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
                           name="rank(mean expression)",
                           expand = c(0,0)) +
        scale_y_continuous(limits = c(NA, ymax),
                           name = "SD") +
        theme_light() +
        ggtitle(title) +
        theme(text = element_text(size=8))
}

get_mean_counts = function(counts_table,
                           samples,
                           groups,
                           condition_id){
    counts_table %>%
        # select(-c(1:6)) %>%
        # rownames_to_column(var = "index") %>%
        gather(sample, value, -index) %>%
        mutate(group = if_else(sample %in% samples[groups==condition_id],
                               "condition_occupancy",
                               "control_occupancy")) %>%
        group_by(index, group) %>%
        summarise(mean = mean(value)) %>%
        spread(group, mean) %>%
        ungroup() %>%
        return()
}

extract_deseq_results = function(dds,
                                 annotations,
                                 mean_counts_table,
                                 alpha,
                                 lfc){
    results(dds,
            alpha=alpha,
            lfcThreshold=lfc,
            altHypothesis="greaterAbs") %>%
        as.data.frame() %>%
        rownames_to_column(var = "index") %>%
        as_tibble() %>%
        left_join(annotations, ., by="index") %>%
        left_join(mean_counts_table, ., by="index") %>%
        arrange(padj) %>%
        mutate(name = if_else(name==".",
                              paste0("peak_", row_number()),
                              name),
               score = as.integer(pmin(-125*log2(padj), 1000))) %>%
        mutate_at(vars(pvalue, padj), ~(-log10(.))) %>%
        mutate_if(is.double, round, 3) %>%
        select(index, chrom, start, end, name, score, strand,
               log2_foldchange=log2FoldChange, lfc_SE=lfcSE,
               stat, log10_pval=pvalue, log10_padj=padj, mean_occupancy=baseMean,
               condition_occupancy, control_occupancy) %>%
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
                   xvar = mean_occupancy,
                   yvar = log2_foldchange,
                   exp_type = "IP",
                   lfc,
                   condition,
                   control){
    xvar = enquo(xvar)
    yvar = enquo(yvar)
    ggplot() +
        geom_hline(yintercept = 0, color="black", linetype="dashed") +
        geom_hline(yintercept = c(-lfc, lfc), color="grey70", linetype="dashed") +
        geom_point(data = df_nonsig,
                   aes(x=!!xvar, y=!!yvar),
                   color="black", alpha=0.3, stroke=0, size=0.7) +
        geom_point(data = df_sig,
                   aes(x=!!xvar, y=!!yvar),
                   color="red", alpha=0.3, stroke=0, size=0.7) +
        scale_x_log10(name="mean of normalized counts") +
        ylab(bquote(.(exp_type) ~ log[2]~frac(.(condition),.(control)))) +
        theme_light() +
        theme(text = element_text(size=8, color="black"),
              axis.text = element_text(color = "black"),
              axis.title.y = element_text(angle=0, hjust=1, vjust=0.5))
}

plot_volcano = function(df = results_df_filtered,
                        xvar = log2_foldchange,
                        yvar = log10_padj,
                        exp_type = "IP",
                        lfc,
                        alpha,
                        condition,
                        control){
    xvar = enquo(xvar)
    yvar = enquo(yvar)
    ggplot() +
        geom_vline(xintercept = 0, color="black", linetype="dashed") +
        geom_vline(xintercept = c(-lfc, lfc), color="grey70", linetype="dashed") +
        geom_point(data = df,
                   aes(x = !!xvar, y = !!yvar),
                   alpha=0.3, stroke=0, size=0.7) +
        geom_hline(yintercept = -log10(alpha), color="red", linetype="dashed") +
        xlab(bquote(.(exp_type) ~ log[2] ~ frac(.(condition), .(control)))) +
        ylab(expression(-log[10] ~ FDR)) +
        theme_light() +
        theme(text = element_text(size=8),
              axis.title.y = element_text(angle=0, hjust=1, vjust=0.5))
}

main = function(exp_table,
                spike_table,
                chip_samples,
                input_samples,
                chip_groups,
                input_groups,
                norm,
                condition,
                control,
                alpha,
                lfc,
                counts_norm_chip_out,
                counts_norm_input_out,
                counts_rlog_chip_out,
                counts_rlog_input_out,
                results_all_out,
                results_up_out,
                results_down_out,
                results_unchanged_out,
                bed_all_out,
                bed_up_out,
                bed_down_out,
                bed_unchanged_out,
                qc_plots_out){

    annotations = read_tsv(exp_table) %>%
        select(1:6) %>%
        rownames_to_column(var="index")

    # initialize DESeq datasets
    dds_chip = initialize_dds(data_path = exp_table,
                              samples = chip_samples,
                              groups = chip_groups,
                              condition_id = condition,
                              control_id = control)

    dds_input = initialize_dds(data_path = exp_table,
                               samples = input_samples,
                               groups = input_groups,
                               condition_id = condition,
                               control_id = control)

    if (norm=="spikenorm"){
        dds_chip_spike = initialize_dds(data_path = spike_table,
                                        samples = chip_samples,
                                        groups = chip_groups,
                                        condition_id = condition,
                                        control_id = control) %>%
            estimateSizeFactors()
        dds_input_spike = initialize_dds(data_path = spike_table,
                                         samples = input_samples,
                                         groups = input_groups,
                                         condition_id = condition,
                                         control_id = control) %>%
            estimateSizeFactors()

        sizeFactors(dds_chip) = sizeFactors(dds_chip_spike)
        sizeFactors(dds_input) = sizeFactors(dds_input_spike)
    } else {
        dds_chip %<>% estimateSizeFactors()
        dds_input %<>% estimateSizeFactors()
    }
    dds_chip %<>% estimateDispersions() %>% nbinomWaldTest()
    dds_input %<>% estimateDispersions() %>% nbinomWaldTest()

    #extract normalized counts and write to file
    counts_norm_chip = extract_normalized_counts(dds = dds_chip)
                                                 # annotations = annotations,
                                                 # output_path = counts_norm_chip_out)
    counts_norm_input = extract_normalized_counts(dds = dds_input)
                                                  # annotations = annotations,
                                                  # output_path = counts_norm_input_out)
    counts_rlog_chip = extract_rlog_counts(dds = dds_chip)
                                           # annotations = annotations,
                                           # output_path = counts_rlog_chip_out)
    counts_rlog_input = extract_rlog_counts(dds = dds_input)
                                           # annotations = annotations,
                                           # output_path = counts_rlog_input_out)

    mean_sd_df_pre_chip = build_mean_sd_df_pre(dds_chip)
    mean_sd_df_pre_input = build_mean_sd_df_pre(dds_input)
    mean_sd_df_post_chip = build_mean_sd_df_post(counts_rlog_chip)
    mean_sd_df_post_input = build_mean_sd_df_post(counts_rlog_input)

    sd_max = max(c(mean_sd_df_pre_chip[["sd"]],
                   mean_sd_df_pre_input[["sd"]],
                   mean_sd_df_post_chip[["sd"]],
                   mean_sd_df_post_input[["sd"]]),
                 na.rm=TRUE)*1.01

    mean_sd_plot_pre_chip = mean_sd_plot(df = mean_sd_df_pre_chip,
                                         ymax = sd_max,
                                         title = expression(IP ~ log[2] ~ "counts," ~ "pre-shrinkage"))
    mean_sd_plot_pre_input = mean_sd_plot(df = mean_sd_df_pre_input,
                                         ymax = sd_max,
                                         title = expression(input ~ log[2] ~ "counts," ~ "pre-shrinkage"))
    mean_sd_plot_post_chip = mean_sd_plot(df = mean_sd_df_post_chip,
                                         ymax = sd_max,
                                         title = expression(regularized ~ IP ~ log[2] ~ "counts"))
    mean_sd_plot_post_input = mean_sd_plot(df = mean_sd_df_post_input,
                                         ymax = sd_max,
                                         title = expression(regularized ~ input ~ log[2] ~ "counts"))

    mean_counts_norm_chip = get_mean_counts(counts_table = counts_norm_chip,
                                            samples = chip_samples,
                                            groups = chip_groups,
                                            condition_id = condition)
    mean_counts_norm_input = get_mean_counts(counts_table = counts_norm_input,
                                             samples = input_samples,
                                             groups = input_groups,
                                             condition_id = condition)

    results_df_chip = extract_deseq_results(dds = dds_chip,
                                            annotations = annotations,
                                            mean_counts_table = mean_counts_norm_chip,
                                            alpha = alpha,
                                            lfc = lfc)
    # set annotation names in input to names assigned from the ChIP analysis
    results_df_input = extract_deseq_results(dds = dds_input,
                                             annotations = annotations,
                                             mean_counts_table = mean_counts_norm_input,
                                             alpha = alpha,
                                             lfc = lfc) %>%
        select(-name) %>%
        left_join(results_df_chip %>%
                      select(index, name),
                  by = "index") %>%
        select(index, chrom, start, end, name, score, strand, everything())

    write_counts_table(results_df = results_df_chip,
                       annotations = annotations,
                       counts_df = counts_norm_chip,
                       output_path = counts_norm_chip_out)
    write_counts_table(results_df = results_df_chip,
                       annotations = annotations,
                       counts_df = counts_rlog_chip,
                       output_path = counts_rlog_chip_out)
    write_counts_table(results_df = results_df_input,
                       annotations = annotations,
                       counts_df = counts_norm_input,
                       output_path = counts_norm_input_out)
    write_counts_table(results_df = results_df_input,
                       annotations = annotations,
                       counts_df = counts_rlog_input,
                       output_path = counts_rlog_input_out)

    #blacklist differentially bound regions that are also differentially bound in the input
    results_df_filtered = results_df_chip %>%
        left_join(results_df_input,
                  by = c("index", "chrom", "start", "end", "name", "strand"),
                  suffix = c("", "_input")) %>%
        select(-index) %>%
        mutate_at(vars(c(5, 10, 11)),
                  ~(if_else((log10_padj > -log10(alpha) &
                                    log10_padj_input > -log10(alpha)) &
                                   (sign(log2_foldchange)==sign(log2_foldchange_input)),
                                     -1, as.numeric(.), missing=as.numeric(.)))) %>%
        mutate(score = as.integer(score)) %>%
        arrange(desc(log10_padj)) %>%
        write_tsv(results_all_out)

    results_df_filtered_significant = results_df_filtered %>%
        filter(log10_padj > -log10(alpha))
    results_df_filtered_nonsignificant = results_df_filtered %>%
        filter(log10_padj <= -log10(alpha)) %>%
        write_tsv(results_unchanged_out)

    results_df_filtered_significant %>%
        filter(log2_foldchange >= 0) %>%
        write_tsv(results_up_out) %>%
        select(1:6) %>%
        write_tsv(bed_up_out, col_names=FALSE)

    results_df_filtered_significant %>%
        filter(log2_foldchange < 0) %>%
        write_tsv(results_down_out) %>%
        select(1:6) %>%
        write_tsv(bed_down_out, col_names=FALSE)

    results_df_filtered %>%
        select(1:6) %>%
        write_tsv(bed_all_out, col_names=FALSE)
    results_df_filtered_nonsignificant %>%
        select(1:6) %>%
        write_tsv(bed_unchanged_out, col_names=FALSE)


    maplot_chip = plot_ma(df_sig = results_df_filtered_significant,
                          df_nonsig = results_df_filtered_nonsignificant,
                          xvar = mean_occupancy,
                          yvar = log2_foldchange,
                          exp_type = "IP",
                          lfc = lfc,
                          condition = condition,
                          control = control)
    maplot_input = plot_ma(df_sig = results_df_filtered_significant,
                           df_nonsig = results_df_filtered_nonsignificant,
                           xvar = mean_occupancy_input,
                           yvar = log2_foldchange_input,
                           exp_type = "input",
                           lfc = lfc,
                           condition = condition,
                           control = control)

    volcano_chip = plot_volcano(df = results_df_filtered,
                                xvar = log2_foldchange,
                                yvar = log10_padj,
                                exp_type = "IP",
                                lfc = lfc,
                                alpha = alpha,
                                condition = condition,
                                control = control)
    volcano_input = plot_volcano(df = results_df_filtered,
                                 xvar = log2_foldchange_input,
                                 yvar = log10_padj_input,
                                 exp_type = "input",
                                 lfc = lfc,
                                 alpha = alpha,
                                 condition = condition,
                                 control = control)

    qc_plots = arrangeGrob(mean_sd_plot_pre_chip,
                           mean_sd_plot_post_chip,
                           maplot_chip,
                           volcano_chip,
                           mean_sd_plot_pre_input,
                           mean_sd_plot_post_input,
                           maplot_input,
                           volcano_input,
                           ncol=4,
                           widths = c(0.7, 0.7, 1, 1))
    ggsave(qc_plots_out,
           plot = qc_plots,
           width = 16*2.5,
           height = 9*2.5,
           units="cm")
}

main(exp_table = snakemake@input[["exp_counts"]],
     spike_table = snakemake@input[["spike_counts"]],
     chip_samples = snakemake@params[["chip_samples"]],
     input_samples = snakemake@params[["input_samples"]],
     chip_groups = snakemake@params[["chip_groups"]],
     input_groups = snakemake@params[["input_groups"]],
     norm = snakemake@wildcards[["norm"]],
     condition = snakemake@wildcards[["condition"]],
     control = snakemake@wildcards[["control"]],
     alpha = snakemake@params[["alpha"]],
     lfc = snakemake@params[["lfc"]],
     counts_norm_chip_out = snakemake@output[["counts_norm_chip"]],
     counts_norm_input_out = snakemake@output[["counts_norm_input"]],
     counts_rlog_chip_out = snakemake@output[["counts_rlog_chip"]],
     counts_rlog_input_out = snakemake@output[["counts_rlog_input"]],
     results_all_out = snakemake@output[["results_all"]],
     results_up_out = snakemake@output[["results_up"]],
     results_down_out = snakemake@output[["results_down"]],
     results_unchanged_out = snakemake@output[["results_unchanged"]],
     bed_all_out = snakemake@output[["bed_all"]],
     bed_up_out = snakemake@output[["bed_up"]],
     bed_down_out = snakemake@output[["bed_down"]],
     bed_unchanged_out = snakemake@output[["bed_unchanged"]],
     qc_plots_out = snakemake@output[["qc_plots"]])

