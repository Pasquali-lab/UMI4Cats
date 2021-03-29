#' Make windows merging restriction fragments
#' 
#' Use a set of continuous restriction fragments to generate windows containing
#' a fixed number of fragments (n_frags).
#' @param input Input object containing the restriction fragments. Should be 
#' class UMI4C (rowRanges will be extracted) or class GRanges.
#' @param n_frags Number of fragments to use for generating the windows. This 
#' should include restriction fragmetns with 0 counts.
#' @param sliding Numeric indicating the factor for generating sliding windows.
#' If set to 1 (default) will use fixed windows. If set to > 0 and < 1 will use
#' n_frags * sliding fragments to generate sliding windows.
#' @return A GRanges object containing the windows of merged restriction 
#' fragments.
#' @export
#' @examples
#' data("ex_ciita_umi4c")
#' 
#' # Without sliding windows
#' win_frags <- makeWindowFragments(ex_ciita_umi4c, n_frags=30, sliding=1)
#' win_frags
#' 
#' # With sliding windows (n_frags*sliding)
#' win_frags <- makeWindowFragments(ex_ciita_umi4c, n_frags=30, sliding=0.5)
#' win_frags
makeWindowFragments <- function(input, 
                                n_frags=30,
                                sliding=1) {
  if (isClass(input, "UMI4C")) frags <- rowRanges(input)
  else if (isClass(input, "GRanges")) frags <- input
  else stop("Input object should be of class 'UMI4C' or 'GRanges'")

  frags <- frags[order(frags)]
  start_upstream <- frags$id_contact[frags$position=="upstream" & start(frags) == max(start(frags[frags$position=="upstream"]))]
  start_downstream <- frags$id_contact[frags$position=="downstream" & start(frags) == min(start(frags[frags$position=="downstream"]))]
  
  if (sliding==1) window_frags <- .makeFixedWindow(frags, n_frags, start_upstream, start_downstream)
  else if (sliding < 1 & sliding > 0) window_frags <- .makeSlidingWindow(frags, n_frags, start_upstream, start_downstream, sliding)
  else stop("Sliding value should be 1 (no sliding) or > 0 & < 1 for sliding.")
  
  return(window_frags)
}

.makeFixedWindow <- function(frags,
                             n_frags,
                             start_upstream,
                             start_downstream) {
  ## Upstream
  up1 <- seq(grep(start_upstream, frags$id_contact), 1, by=-n_frags)
  up2 <- c(up1[-1]+1, 1)
  frags_upstream <- GRanges(seqnames = unique(seqnames(frags)),
                            ranges = IRanges(start=start(frags[up2]), end=end(frags[up1])),
                            mcols= data.frame("id"=paste0("window_UP_", seq_len(length(up1))),
                                              "position"="upstream"))
  
  ## Downstream
  dwn1 <- seq(grep(start_downstream, frags$id_contact), length(frags), by=n_frags)
  dwn2 <- c(dwn1[-1]-1, length(frags))
  frags_downstream <- GRanges(seqnames = unique(seqnames(frags)),
                              ranges = IRanges(start=start(frags[dwn1]), end=end(frags[dwn2])),
                              mcols= data.frame("id"=paste0("window_DOWN_", seq_len(length(dwn1))),
                                                "position"="downstream"))
  
  ## Merge
  window_frags <- c(frags_upstream, frags_downstream)
  window_frags <- window_frags[order(window_frags)]
  
  return(window_frags)
}

.makeSlidingWindow <- function(frags,
                               n_frags,
                               start_upstream,
                               start_downstream,
                               sliding) {
  
  ## Upstream
  up1 <- seq(grep(start_upstream, frags$id_contact), 1, by=-(n_frags*sliding))
  up2 <- up1 - n_frags + 1
  up2[up2<=0] <- 1
  
  frags_upstream <- GRanges(seqnames = unique(seqnames(frags)),
                            ranges = IRanges(start=start(frags[up2]), end=end(frags[up1])),
                            mcols= data.frame("id"=paste0("window_UP_", seq_len(length(up1))),
                                              "position"="upstream"))
  
  ## Downstream
  dwn1 <- seq(grep(start_downstream, frags$id_contact), length(frags), by=n_frags*sliding)
  dwn2 <-dwn1 + n_frags - 1
  dwn2[dwn2 > length(frags)] <- length(frags)
  
  frags_downstream <- GRanges(seqnames = unique(seqnames(frags)),
                              ranges = IRanges(start=start(frags[dwn1]), end=end(frags[dwn2])),
                              mcols= data.frame("id"=paste0("window_DOWN_", seq_len(length(dwn1))),
                                                "position"="downstream"))
  
  ## Merge
  window_frags <- c(frags_upstream, frags_downstream)
  window_frags <- window_frags[order(window_frags)]
  
  return(window_frags)
}