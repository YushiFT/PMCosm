#' Traditional categories of microbial taxa based on the relative abundance
#'
#' This function classify microbial taxa into six categories -
#' 1) RT: rare taxa;
#' 2) IT: intermediate taxa;
#' 3) AT: abundant taxa;
#' 4) CRT: conditionally rare taxa;
#' 5) CAT: conditionally abundant taxa;
#' 6) CRAT: conditionally rare or abundant taxa -
#' by following four popular methods summerized in Tang et al. (2022) including -
#' 1) Method 0 (Dai et al. 2016);
#' 2) Method a (Campbell et al. 2011);
#' 3) Method b (Logares et al. 2014);
#' 4) Method c (Dai et al. 2016)
#'
#'
#' @param x  The n by m matrix of microbial abundance.
#' @param t_lower The lower threshold to separate rare and intermediate taxa.
#' @param t_upper The upper threshold to separate intermediate and abundant taxa.
#' @param method The method index which can be 'a', 'b', 'c', or '0'.
#'
#' @return A list of microbial taxa categories for a given method.
#'
#' @examples
#' # import data
#' data(hzmicrobe)
#' # try method 0 where the lower and upper thresholds can be adjusted
#' category_method_0 <- classify_taxa(mic_bay, t_l=0.0001, t_u=0.01, method='0')
#' # print abundant taxa id
#' print(category_method_0$AT)
#' # similar for printing out other category groups
#' # try method a where the lower and upper thresholds are fixed
#' category_method_a <- classify_taxa(mic_bay, method='a')
#' print(category_method_a$AT)
#' # try method b
#' category_method_b <- classify_taxa(mic_bay, method='b')
#' print(category_method_b$AT)
#' # try method c
#' category_method_c <- classify_taxa(mic_bay, method='c')
#' print(category_method_c$AT)
#'
#'
#' @export
#'
classify_taxa <- function(x, t_lower=0.0001, t_upper=0.01, method='0'){
  if (!(method %in% c('0','a','b','c'))){
    stop("method unavailable, should be either one of 0, a, b, or c.")
  }

  if(method=='0'){
    # extract rare taxa
    # rare taxa has relative abundance no greater than t_upper in all samples
    x_prop <- x/colSums(x) # transform to relative abundance
    f_rt <- function(v){
      record_rt <- sum(v<=t_lower)
      return(record_rt==length(v))
    }
    is_rt_vec <- apply(x_prop, 1, f_rt)
    id_rt <- rownames(x)[is_rt_vec==TRUE]

    # extract intermediate taxa
    # intermediate taxa has relative abundance between t_lower and t_upper in all samples
    f_it <- function(v){
      record_it <- sum((v>t_lower)&(v<=t_upper))
      return(record_it==length(v))
    }
    is_it_vec <- apply(x_prop, 1, f_it)
    id_it <- rownames(x)[is_it_vec==TRUE]

    # extract abundant taxa
    # abundant taxa has relative abundance greater than t_upper in all samples
    f_at <- function(v){
      record_at <- sum(v>t_upper)
      return(record_at==length(v))
    }
    is_at_vec <- apply(x_prop, 1, f_at)
    id_at <- rownames(x)[is_at_vec==TRUE]

    # extract conditionally rare taxa
    # conditionally rare taxa has relative abundance no greater than t_upper in all samples
    # and no greater than 0.01% in some samples
    f_crt <- function(v){
      record_1 <- sum(v<=t_upper)
      record_2 <- sum((v>t_lower)&(v<=t_upper))
      return( (record_1==length(v))&(record_2>0)&(record_2<length(v)) )
    }
    is_crt_vec <- apply(x_prop, 1, f_crt)
    id_crt <- rownames(x)[is_crt_vec==TRUE]

    # extract conditionally abundant taxa
    # conditionally abundant taxa has relative abundance greater than t_lower in all samples
    # and greater than t_upper in some samples
    f_cat <- function(v){
      record_1 <- sum(v>t_lower)
      record_2 <- sum(v>t_upper)
      return( (record_1==length(v))&(record_2>0)&(record_2<length(v)) )
    }
    is_cat_vec <- apply(x_prop, 1, f_cat)
    id_cat <- rownames(x)[is_cat_vec==TRUE]

    # extract conditionally rare or abundant taxa
    # conditionally rare or abundant taxa has relative abundance varying from no greater than t_lower
    # to greater than t_upper
    f_crat <- function(v){
      record_1 <- sum(v<=t_lower)
      record_2 <- sum(v>t_upper)
      return( (record_1>0)&(record_2>0) )
    }
    is_crat_vec <- apply(x_prop, 1, f_crat)
    id_crat <- rownames(x)[is_crat_vec==TRUE]

    taxa_list <- list(id_rt, id_it, id_at, id_crt, id_cat, id_crat)
    names(taxa_list) <- c('RT','IT','AT','CRT','CAT','CRAT')
    return(taxa_list)

  } else if (method=='a'){
    x_prop <- x/colSums(x) # transform to relative abundance

    # extract rare taxa according to method a
    # rare taxa has relative abundance less than 1% in all samples
    f_rt <- function(v){
      record_rt <- sum(v<0.01)
      return(record_rt==length(v))
    }
    is_rt_vec <- apply(x_prop, 1, f_rt)
    id_rt <- rownames(x)[is_rt_vec==TRUE]

    # extract rare taxa according to method a
    # rare taxa has relative abundance no less than 1% in any sample
    f_at <- function(v){
      record_at <- sum(v>=0.01)
      return(record_at>0)
    }
    is_at_vec <- apply(x_prop, 1, f_at)
    id_at <- rownames(x)[is_at_vec==TRUE]

    taxa_list <- list(id_rt, id_at)
    names(taxa_list) <- c('RT','AT')
    return(taxa_list)

  } else if (method=='b'){

    x_prop <- x/colSums(x) # transform to relative abundance
    v_prop <- rowMeans(x_prop)

    # extract regional rare taxa according to method b
    # regional rare taxa has mean relative abundance less than 0.001%
    id_rt <- rownames(x)[v_prop<0.00001]

    # extract regional abundant taxa according to method b
    # regional abundant taxa has mean relative abundance greater than 0.1%
    id_at <- rownames(x)[v_prop>0.001]

    # extract regional intermediate taxa according to method b
    # regional intermediate taxa has mean relative abundance between 0.001% and 0.1%
    id_it <- rownames(x)[(v_prop>=0.00001)&(v_prop<=0.001)]

    taxa_list <- list(id_rt, id_it, id_at)
    names(taxa_list) <- c('RT','IT','AT')
    return(taxa_list)

  } else if (method=='c') {
    x_prop <- x/colSums(x) # transform to relative abundance

    # extract rare taxa according to method c
    # rare taxa has relative abundance less than 0.1% in all samples
    f_rt <- function(v){
      record_rt <- sum(v<0.001)
      return(record_rt==length(v))
    }
    is_rt_vec <- apply(x_prop, 1, f_rt)
    id_rt <- rownames(x)[is_rt_vec==TRUE]

    # extract abundant taxa according to method c
    # abundant taxa has relative abundance no less than 0.1% in half or more samples
    # and occurred in more than 80% samples
    f_at <- function(v){
      record_at_1 <- sum(v>=0.001)
      record_at_2 <- sum(v>0)
      return((record_at_1>0.5*length(v)) & (record_at_2>0.8*length(v)))
    }
    is_at_vec <- apply(x_prop, 1, f_at)
    id_at <- rownames(x)[is_at_vec==TRUE]

    # extract intermediate taxa according to method c
    # intermediate taxa is neither rt or at
    id_it <- rownames(x)[(is_rt_vec==FALSE)&(is_at_vec==FALSE)]

    taxa_list <- list(id_rt, id_it, id_at)
    names(taxa_list) <- c('RT','IT','AT')
    return(taxa_list)
  }
}
