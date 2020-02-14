library(tidyverse)
library(magrittr)
library(ggridges)
library(ggthemes)

main = function(annotation_paths = c("Scer_transcripts_w_verifiedORFs-nonoverlapping_slopR300.bed"),
                denominator_path = "verified-transcripts-nonoverlapping-slopR300-TSS-allsamples-allannotations-chipseq-spikenorm-midpoints-input-subtracted.tsv.gz",
                ratio_path = "verified-transcripts-nonoverlapping-slopR300-TSS-allsamples-allannotations-chipseq-spikenorm-ratio.tsv.gz",
                condition = "depleted",
                control = "non-depleted",
                heatmap_colormap = "inferno",
                ridgeline_colormap = "viridis",
                colorbar_quantile_low = 0.02,
                colorbar_quantile_high = 0.99999,
                length_sort = TRUE,
                annotations_out = c("anno.bed"),
                nfactor = "H3K36me3",
                dfactor = "H3",
                distance_type = "absolute",
                scaled_length = 0,
                refpt_label = "TSS",
                endpt_label = "HAIL SATAN",
                ratio_heatmap_out = "ratio_heatmap.pdf",
                ratio_metagene_out = "ratio_metagene.pdf",
                ridgelines_out = "ridgelines.pdf"){

    df_ratio = read_tsv(ratio_path,
                        col_names=c("group", "group_duplicate", "annotation", "index", "position", "ratio")) %>%
        select(-group_duplicate) %>%
        filter(group %in% c(control, condition)) %>%
        mutate(annotation = fct_inorder(annotation, ordered=TRUE),
               group = ordered(group,
                               levels=c(control,
                                        condition))) %>%
        group_by(annotation) %>%
        mutate(annotation_withcount = fct_inorder(paste(n_distinct(index),
                                                        annotation),
                                                  ordered=TRUE)) %>%
        ungroup()

    annotation_labels = unique(df_ratio[["annotation"]])
    n_annotations = length(annotation_labels)
    annotations = tibble()

    for (i in 1:n_annotations){
        annotations = read_tsv(annotation_paths[i],
                               col_names=c("chrom",
                                           "start",
                                           "end",
                                           "name",
                                           "score",
                                           "strand")) %>%
            mutate(annotation = annotation_labels[i]) %>%
            rowid_to_column(var = "index") %>%
            bind_rows(annotations, .)
    }

    if (length_sort){
        annotations %<>%
            group_by(annotation) %>%
            arrange(end - start,
                    .by_group=TRUE) %>%
            rowid_to_column(var = "sorted_index") %>%
            mutate(sorted_index = as.integer(sorted_index + 1 - min(sorted_index)))

        df_ratio %<>%
            left_join(select(annotations,
                             annotation,
                             index,
                             sorted_index),
                      by=c("annotation",
                           "index")) %>%
            select(-index) %>%
            rename(index = sorted_index)
    }

    for (i in 1:n_annotations){
        annotations %>%
            filter(annotation == annotation_labels[i]) %>%
            ungroup() %>%
            select(chrom, start, end, name, score, strand) %>%
            write_tsv(annotations_out[i], col_names=FALSE)
        }


    heatmap_colorbar_limits = quantile(df_ratio[["ratio"]],
                                       probs=c(colorbar_quantile_low,
                                               colorbar_quantile_high))


    ratio_heatmap = ggplot(data=df_ratio,
           aes(x=position,
               y=index,
               fill=ratio)) +
        geom_raster(interpolate=TRUE) +
        facet_grid(annotation_withcount ~ group,
                   switch="y") +
        scale_y_reverse(expand=c(0, 0),
                        name=NULL) +
        scale_fill_viridis_c(option=heatmap_colormap,
                             limits=heatmap_colorbar_limits,
                             breaks=scales::pretty_breaks(3),
                             oob=scales::squish,
                             guide=guide_colorbar(title.position="top",
                                                  title.hjust=0.5,
                                                  barwidth=12,
                                                  barheight=0.5),
                             name=bquote("relative enrichment:" ~
                                             textstyle(frac(.(nfactor),
                                                            .(dfactor))))) +
        theme_light() +
        theme(text=element_text(size=12),
              strip.background=element_blank(),
              strip.text=element_text(color="black"),
              strip.placement="outside",
              axis.text=element_text(color="black"),
              axis.text.y=element_blank(),
              panel.border=element_blank(),
              panel.spacing.x=unit(12, "pt"),
              legend.position="top",
              legend.margin=margin(0, 0, -10, 0, "pt"),
              legend.box.margin=margin(0, 0, 0, 0),
              legend.title=element_text(size=12,
                                        margin=margin(b=-4, unit="pt")),
              legend.text=element_text(margin=margin(t=-2, unit="pt")),
              plot.margin=margin(0, 12, 0, 0, "pt"))
    if (distance_type == "absolute"){
        ratio_heatmap = ratio_heatmap +
            scale_x_continuous(expand=c(0, 0),
                               breaks=scales::pretty_breaks(3),
                               labels=function(x) ifelse(x==0,
                                                         refpt_label,
                                                         paste(x, "kb")),
                               name=NULL)
    } else {
        ratio_heatmap = ratio_heatmap +
            scale_x_continuous(expand=c(0, 0),
                               breaks=c(0, scaled_length/2, scaled_length),
                               labels=c(refpt_label, "", endpt_label),
                               name=NULL)
    }

    ggsave(ratio_heatmap_out,
           plot=ratio_heatmap,
           width=16,
           height=9,
           units="cm")

    df_ratio_metagene = df_ratio %>%
        group_by(group,
                 annotation,
                 annotation_withcount,
                 position) %>%
        summarize(low = quantile(ratio, 0.25),
                  mid = median(ratio),
                  high = quantile(ratio, 0.75))

    ratio_metagene = ggplot(data=df_ratio_metagene,
           aes(x=position,
               y=mid,
               ymin=low,
               ymax=high,
               color=group,
               fill=group))
    if (distance_type == "absolute"){
        ratio_metagene = ratio_metagene +
            geom_vline(xintercept=0,
                       color="grey70",
                       size=0.2) +
            scale_x_continuous(expand=c(0, 0),
                               breaks=scales::pretty_breaks(3),
                               labels=function(x) ifelse(x==0,
                                                         refpt_label,
                                                         paste(x, "kb")),
                               name=NULL)
    } else {
        ratio_metagene = ratio_metagene +
            geom_vline(xintercept=c(0, scaled_length),
                       color="grey70",
                       size=0.2) +
            scale_x_continuous(expand=c(0, 0),
                               breaks=c(0, scaled_length/2, scaled_length),
                               labels=c(refpt_label, "", endpt_label),
                               name=NULL)
    }
    ratio_metagene = ratio_metagene +
        geom_ribbon(alpha=0.1,
                    linetype="blank") +
        geom_line(alpha=0.9,
                  size=1) +
        facet_wrap(~annotation_withcount) +
        scale_y_continuous(breaks=scales::pretty_breaks(2),
                           name=bquote(atop("relative enrichment:",
                                             textstyle(frac(.(nfactor),
                                                            .(dfactor)))))) +
        scale_fill_ptol() +
        scale_color_ptol() +
        theme_light() +
        theme(strip.background = element_blank(),
              strip.text = element_text(color="black",
                                        hjust=0),
              panel.grid=element_blank(),
              legend.title=element_blank(),
              legend.spacing.x=unit(2, "pt"),
              axis.title.y=element_text(angle=0,
                                        vjust=0.66,
                                        hjust=1),
              axis.text=element_text(color="black"),
              legend.position=c(-0.35,0.33),
              legend.justification=c(0,0.5),
              plot.margin=margin(0, 12, 0, 0, "pt"))

    ggsave(ratio_metagene_out,
           plot=ratio_metagene,
           width=16,
           height=9,
           units="cm")

    df_denominator = read_tsv(denominator_path,
                              col_names=c("group",
                                          "sample",
                                          "annotation",
                                          "index",
                                          "position",
                                          "denominator")) %>%
        filter(group %in% c(control, condition)) %>%
        group_by(group, annotation, index, position) %>%
        summarize(denominator = mean(denominator, na.rm=TRUE)) %>%
        group_by(group, annotation, position) %>%
        summarize(low=quantile(denominator, 0.25),
                  mid=median(denominator),
                  high=quantile(denominator, 0.75)) %>%
        ungroup() %>%
        mutate(group = ordered(group,
                               levels=c(control,
                                        condition)))

    df_ridges = df_ratio_metagene %>%
        left_join(df_denominator,
                  by=c("group",
                       "annotation",
                       "position"),
                  suffix=c("_ratio", "_denominator")) %>%
        ungroup() %>%
        mutate(mid = scales::rescale(mid_denominator))

    ridgelines = ggplot(data=df_ridges,
           aes(x=position,
               y=fct_rev(group),
               height=mid,
               fill=mid_ratio))
    if (distance_type == "absolute"){
        ridgelines = ridgelines +
            geom_vline(xintercept=0,
                       color="grey70",
                       size=0.2) +
            scale_x_continuous(expand=c(0, 0),
                               breaks=scales::pretty_breaks(3),
                               labels=function(x) ifelse(x==0,
                                                         refpt_label,
                                                         paste(x, "kb")),
                               name=NULL)
    } else {
        ridgelines = ridgelines +
            geom_vline(xintercept=c(0, scaled_length),
                       color="grey70",
                       size=0.2) +
            scale_x_continuous(expand=c(0, 0),
                               breaks=c(0, scaled_length/2, scaled_length),
                               labels=c(refpt_label, "", endpt_label),
                               name=NULL)
    }
    ridgelines = ridgelines +
        geom_ridgeline_gradient(scale=2,
                                # color=viridis::viridis(3, option=ridgeline_colormap)[2],
                                size=0.8) +
        facet_wrap(~annotation_withcount) +
        scale_fill_viridis_c(option=ridgeline_colormap,
                             breaks=scales::pretty_breaks(3),
                             oob=scales::squish,
                             guide=guide_colorbar(title.position="top",
                                                  title.hjust=0.5,
                                                  barwidth=12,
                                                  barheight=0.5),
                             name=bquote("relative enrichment:" ~
                                             textstyle(frac(.(nfactor),
                                                            .(dfactor))))) +
        scale_y_discrete(expand=c(0,0),
                         name=bquote("relative" ~ .(dfactor) ~ "enrichment")) +
        theme_light() +
        theme(axis.text=element_text(color="black"),
              axis.text.y=element_text(vjust=0),
              strip.background = element_blank(),
              strip.text=element_text(color="black",
                                      hjust=0),
              panel.border=element_blank(),
              legend.position="top",
              legend.margin=margin(0, 0, -10, 0, "pt"),
              legend.box.margin=margin(0, 0, 0, 0),
              legend.title=element_text(size=12,
                                        margin=margin(b=-4, unit="pt")),
              legend.text=element_text(margin=margin(t=-2, unit="pt")),
              plot.margin=margin(0, 12, 0, 0, "pt"))

    ggsave(ridgelines_out,
           plot=ridgelines,
           width=16,
           height=9,
           units="cm")
}

main(annotation_paths = snakemake@input[["annotations"]],
     denominator_path = snakemake@input[["denominator"]],
     ratio_path = snakemake@input[["ratio"]],
     condition = snakemake@wildcards[["condition"]],
     control = snakemake@wildcards[["control"]],
     heatmap_colormap = snakemake@params[["heatmap_colormap"]],
     ridgeline_colormap = snakemake@params[["ridgeline_colormap"]],
     colorbar_quantile_low = snakemake@params[["colorbar_quantile_low"]],
     colorbar_quantile_high = snakemake@params[["colorbar_quantile_high"]],
     length_sort = snakemake@params[["length_sort"]],
     nfactor = snakemake@wildcards[["nfactor"]],
     dfactor = snakemake@wildcards[["dfactor"]],
     distance_type = snakemake@params[["distance_type"]],
     scaled_length = snakemake@params[["scaled_length"]] / 1000,
     refpt_label = snakemake@params[["refpt_label"]],
     endpt_label = snakemake@params[["endpt_label"]],
     annotations_out = snakemake@params[["annotations_out"]],
     ratio_heatmap_out = snakemake@output[["ratio_heatmap"]],
     ratio_metagene_out = snakemake@output[["ratio_metagene"]],
     ridgelines_out = snakemake@output[["ridgelines"]])

