.gg.plot <- function(...) {
    ggplot(...) +
        theme_classic() +
        theme(axis.title = element_text(size=8)) +
        theme(axis.text = element_text(size=6)) +
        theme(legend.spacing = unit(.1, "lines"),
              legend.key.size = unit(.5, "lines"),
              legend.text = element_text(size=5),
              legend.title = element_text(size=6),
              panel.background = element_rect(fill='transparent'),
              plot.background = element_rect(fill='transparent', color=NA),
              legend.background = element_rect(fill='transparent', size=0.05),
              legend.box.background = element_rect(fill='transparent'))
}

row.order <- function(mat) {
    require(cba)
    require(proxy)

    if(nrow(mat) < 3) {
        return(1:nrow(mat))
    }

    D = proxy::dist(mat, method <- function(a,b) 1 - cor(a,b, method = 'spearman'))
    D[!is.finite(D)] = 0
    h.out = hclust(D)
    o.out = cba::order.optimal(D, h.out$merge)
    return(o.out$order)
}

col.order <- function(pair.tab, .ro, ret.tab = FALSE) {
    require(tidyr)
    require(dplyr)

    M = pair.tab %>%
        dplyr::select(row, col, weight) %>%
        dplyr::mutate(row = factor(row, .ro)) %>%
        tidyr::spread(key = col, value = weight, fill = 0)

    co = order(apply(M[, -1], 2, which.max), decreasing = TRUE)
    .co = colnames(M)[-1][co]
    if(ret.tab) {
        ret = pair.tab %>%
            dplyr::mutate(row = factor(row, .ro)) %>%
            dplyr::mutate(col = factor(col, .co))
    } else {
        ret = .co
    }
    return(ret)
}

order.pair <- function(pair.tab, ret.tab=FALSE) {

    require(tidyr)
    require(dplyr)

    .tab = pair.tab %>% dplyr::select(row, col, weight)

    M = .tab %>% tidyr::spread(key = col, value = weight, fill = 0)
    rr = M[, 1] %>% unlist(use.names = FALSE)
    cc = colnames(M)[-1] %>% unlist(use.names = FALSE)

    ## log.msg('Built the Mat: %d x %d', nrow(M), ncol(M))
    ro = row.order(M %>% dplyr::select(-row) %>% as.matrix())

    ## log.msg('Sort the rows: %d', length(ro))
    co = order(apply(M[ro, -1], 2, which.max), decreasing = TRUE)

    ## co = row.order(t(M %>% dplyr::select(-row) %>% as.matrix()))
    ## log.msg('Sort the columns: %d', length(co))

    if(ret.tab){
        ret = pair.tab %>%
            dplyr::mutate(row = factor(row, rr[ro])) %>%
            dplyr::mutate(col = factor(col, cc[co]))
    } else {
        ret = list(rows = rr[ro], cols = cc[co], M = M)
    }

    return(ret)
}