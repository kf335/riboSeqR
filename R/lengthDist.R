lengthDist <- function(riboDat, add = FALSE, legend = NULL, ...)
  {
    lens <- lapply(riboDat@riboGR, function(gr) table(width(gr)) / length(gr))
    if(!add) {
#      xlim <- range(as.numeric(unlist(lapply(lens, names))))
#      if("ylim" %in% names(list(...))) ylim <- list(...)$ylim else ylim <- c(0, 1)
      do.call("plot", modifyList(list(x = NA, y = NA, xlim = range(as.numeric(unlist(lapply(lens, names)))), ylim = c(0, 1), ylab = "Proportion of aligned reads", xlab = "Length of read"), list(...)))
    }
    args <- modifyList(list(col = rainbow(length(riboDat@riboGR)), lwd = 1, lty = 1), list(...))
    args$col <- rep(args$col, length(riboDat@riboGR))
    args$lty <- rep(args$lty, length(riboDat@riboGR))

    for(ii in 1:length(lens)) {
      cargs <- args
      cargs$col <- cargs$col[ii]
      cargs$lty <- cargs$lty[ii]
      cargs[!names(cargs) %in% c("x", "y")]
      do.call(lines, c(list(y = as.numeric(lens[[ii]]), x = as.numeric(names(lens[[ii]]))), cargs))
    }
  
    if(!is.null(legend)) do.call(graphics::legend, modifyList(args[!names(args) %in% c("ylim", "cex.lab", "cex.axis")], list(x = "topright", legend = legend)))
  }
