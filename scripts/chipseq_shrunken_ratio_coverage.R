library(tidyverse)
library(magrittr)
library(DESeq2)
library(ashr)
library(gridExtra)

get_countdata = function(path_denominator="H3_chipseq_allsamples-counts-midpoints-window-200.tsv.gz",
                         path_numerator="H3K36me3_chipseq_allsamples-counts-midpoints-window-200.tsv.gz",
                         samples_denominator,
                         samples_numerator){
    df = read_tsv(path_denominator) %>%
        select(1, samples_denominator) %>%
        magrittr::set_colnames(c("name",
                                 paste0(samples_denominator, "_D"))) %>%
        full_join(read_tsv(path_numerator) %>%
                      select(1, samples_numerator) %>%
                      magrittr::set_colnames(c("name",
                                               paste0(samples_numerator, "_N"))),
                  by=c("name")) %>%
        select(-c(1)) %>%
        mutate_all(~replace_na(., 0)) %>%
        rowid_to_column(var="index") %>%
        column_to_rownames(var="index") %>%
        as.data.frame()
    # df = df[rowSums(df)>1,]
    return(df)
}

initialize_dds = function(data_path_denominator,
                          data_path_numerator,
                          samples_denominator,
                          samples_numerator,
                          rna_sources_denominator,
                          rna_sources_numerator){
    dds = DESeqDataSetFromMatrix(countData = get_countdata(path_denominator=data_path_denominator,
                                                           path_numerator=data_path_numerator,
                                                           samples_denominator=samples_denominator,
                                                           samples_numerator=samples_numerator),
                                 colData = data.frame(rna_source = factor(c(rna_sources_denominator,
                                                                            rna_sources_numerator),
                                                                          levels = c("input",
                                                                                     "ChIP")),
                                                      chip_factor = factor(c(rep("denominator",
                                                                                 length(samples_denominator)),
                                                                             rep("numerator",
                                                                                 length(samples_numerator))),
                                                                           levels=c("denominator",
                                                                                    "numerator")),
                                                      row.names = c(paste0(samples_denominator, "_D"),
                                                                    paste0(samples_numerator, "_N"))),
                                 design = ~ rna_source * chip_factor)
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
                           coef="rna_sourceChIP.chip_factornumerator",
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

main = function(exp_table_denominator="H3_chipseq_allsamples-counts-midpoints-window-200.tsv.gz",
                exp_table_numerator="H3K36me3_chipseq_allsamples-counts-midpoints-window-200.tsv.gz",
                spike_table_denominator="depleted-v-non-depleted_allsamples-spikein-H3-chipseq-counts-verified-coding-genes.tsv.gz",
                spike_table_numerator="depleted-v-non-depleted_allsamples-spikein-H3K36me3-chipseq-counts-verified-coding-genes.tsv.gz",
                samples_denominator=read_tsv(exp_table_denominator) %>% names() %>% extract(c(2:5, 10:13)),
                samples_numerator=read_tsv(exp_table_numerator) %>% names() %>% extract(c(2,3,6,7)),
                rna_sources_denominator=c(rep("input",4), rep("ChIP", 4)),
                rna_sources_numerator=c(rep("input",2), rep("ChIP", 2)),
                nfactor="H3K36me3",
                dfactor="H3",
                norm="libsizenorm",
                counts_norm_out="counts_norm.tsv.gz",
                counts_rlog_out="counts_rlog.tsv.gz",
                results_all_out="results_all.tsv.gz",
                bedgraph_out = "test.bedgraph",
                qc_plots_out="qcplots.png"){

    annotations = read_tsv(exp_table_denominator) %>%
        select(name) %>%
        full_join(read_tsv(exp_table_numerator) %>%
                      select(name),
                  by="name") %>%
        rownames_to_column(var="index") %>%
        separate(name,
                 into=c("chrom", "start", "end"),
                 sep="-",
                 convert=TRUE)

    dds = initialize_dds(data_path_denominator=exp_table_denominator,
                         data_path_numerator=exp_table_numerator,
                         samples_denominator,
                         samples_numerator,
                         rna_sources_denominator,
                         rna_sources_numerator)

    if (norm=="spikenorm"){
        dds_spike = initialize_dds(data_path_denominator=spike_table_denominator,
                                   data_path_numerator=spike_table_numerator,
                                   samples_denominator,
                                   samples_numerator,
                                   rna_sources_denominator,
                                   rna_sources_numerator) %>%
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
        select(chrom, start, end, log2FoldChange_shrunken) %>%
        replace_na(list("log2FoldChange_shrunken"=0)) %>%
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

main(exp_table_denominator = snakemake@input[["exp_table_denominator"]],
     exp_table_numerator = snakemake@input[["exp_table_numerator"]],
     spike_table_denominator = snakemake@input[["spike_table_denominator"]],
     spike_table_numerator = snakemake@input[["spike_table_numerator"]],
     samples_denominator = snakemake@params[["samples_denominator"]],
     samples_numerator = snakemake@params[["samples_numerator"]],
     rna_sources_denominator = snakemake@params[["rna_sources_denominator"]],
     rna_sources_numerator = snakemake@params[["rna_sources_numerator"]],
     nfactor=snakemake@wildcards[["nfactor"]],
     dfactor=snakemake@wildcards[["dfactor"]],
     norm = snakemake@wildcards[["norm"]],
     counts_norm_out = snakemake@output[["counts_norm"]],
     counts_rlog_out = snakemake@output[["counts_rlog"]],
     results_all_out = snakemake@output[["tsv"]],
     bedgraph_out = snakemake@output[["bedgraph"]],
     qc_plots_out = snakemake@output[["qc_plots"]])

