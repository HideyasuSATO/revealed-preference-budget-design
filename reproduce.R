#!/usr/bin/env Rscript
# Run: Rscript --vanilla reproduce.R [--representative|--all]

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) == 0L) "--all" else args[[1L]]
if (!mode %in% c("--all", "--representative")) {
  stop("Usage: Rscript --vanilla reproduce.R [--representative|--all]")
}
cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
if (length(file_arg) != 1L) stop("Run this file with Rscript.")
root <- dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE))
for (pkg in c("gmp", "jsonlite")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing package ", pkg, ". See README.txt; no package is installed automatically.")
  }
}
suppressPackageStartupMessages(library(gmp))
Q <- function(x, y = NULL) {
  if (is.null(y)) as.bigq(x) else as.bigq(as.bigz(x), as.bigz(y))
}
Z <- as.bigz
q_floor <- function(x) numerator(x) %/% denominator(x)
q_ceil <- function(x) -((-numerator(x)) %/% denominator(x))
q_text <- function(x) as.character(x)
z_text <- function(x) as.character(x)

isqrt_z <- function(n) {
  if (n < 0) stop("isqrt requires a nonnegative integer")
  if (n < 2) return(n)
  x <- Z(10)^ceiling(nchar(z_text(n)) / 2)
  repeat {
    y <- (x + n %/% x) %/% 2L
    if (y >= x) {
      stopifnot(x*x <= n, (x+1L)*(x+1L) > n)
      return(x)
    }
    x <- y
  }
}

# Round to nearest, ties to even.
decimal_text <- function(x, digits = 80L) {
  if (x == 0) return("0")
  if (x < 0) return(paste0("-", decimal_text(-x, digits)))
  p <- numerator(x); d <- denominator(x)
  exponent <- nchar(z_text(p)) - nchar(z_text(d))
  ten_power <- function(k) if (k >= 0L) Q(Z(10)^k) else Q(1L, Z(10)^(-k))
  if (x < ten_power(exponent)) exponent <- exponent - 1L
  while (x >= ten_power(exponent + 1L)) exponent <- exponent + 1L
  scaled <- x * ten_power(digits - 1L - exponent)
  rounded <- q_floor(scaled)
  rem <- numerator(scaled) - rounded * denominator(scaled)
  twice <- 2L * rem
  if (twice > denominator(scaled) ||
      (twice == denominator(scaled) && rounded %% 2L == 1L)) rounded <- rounded + 1L
  txt <- z_text(rounded)
  if (nchar(txt) > digits) {
    rounded <- rounded %/% 10L
    exponent <- exponent + 1L
    txt <- z_text(rounded)
  }
  txt <- paste0(strrep("0", digits - nchar(txt)), txt)
  paste0(substr(txt, 1L, 1L), ".", substr(txt, 2L, digits), "e", exponent)
}

q_positive <- function(s, n, R, eps) {
  if (!(s >= 2L && n >= 1L && n <= s %/% 2L && R > 1 && eps >= 0)) {
    stop("Invalid model parameters")
  }
  a <- R - 1L; b <- 1L - 1L / R
  if (eps >= b) return(Q(0L))
  A <- a + eps; B <- b - eps; m <- s - 2L*n
  x <- A / (A+B); y <- eps / b
  l_term <- Q(1L)
  prefix <- list(Q(1L))
  if (n > 1L) for (ell in seq_len(n - 1L)) {
    l_term <- if (m > 0L) l_term * Q(m + ell - 1L, ell) * y else Q(0L)
    prefix[[ell + 1L]] <- prefix[[ell]] + l_term
  }
  total <- Q(0L); j_term <- Q(1L)
  for (j in 0:(n - 1L)) {
    if (j > 0L) j_term <- j_term * Q(n + j - 1L, j) * x
    total <- total + j_term * prefix[[n - j]]
  }
  (B / (A+B))^n * (B / b)^m * total
}

log_interval <- function(x, digits) {
  if (x < 1L) stop("Positive log series requires x >= 1")
  if (x == 1L) return(list(lo=Q(0L), hi=Q(0L), terms=0L))
  z <- (x-1L)/(x+1L)
  stopifnot(z > 0L, z < 1L)
  z2 <- z*z; term <- z; total <- Q(0L); k <- 0L
  tolerance <- Q(1L, Z(10)^(digits + 8L))
  repeat {
    total <- total + 2L*term / (2L*k + 1L)
    k <- k + 1L
    term <- term*z2
    remainder <- 2L*term / ((2L*k+1L)*(1L-z2))
    if (remainder <= tolerance) {
      return(list(lo=total, hi=total+remainder, terms=k))
    }
    if (k > 10000L) stop("Series did not converge within the safety limit")
  }
}

outward <- function(lo, hi, digits=80L) {
  scale <- Z(10)^digits
  list(lo=Q(q_floor(lo*scale), scale), hi=Q(q_ceil(hi*scale), scale))
}

target_interval <- function(alpha, digits=80L) {
  if (alpha != Q(4L,5L)) stop("Only the frozen alpha = 4/5 is implemented")
  scale <- Z(10)^(digits+12L)
  root_floor <- isqrt_z(5L*scale*scale)
  sqrt_lo <- Q(root_floor,scale); sqrt_hi <- Q(root_floor+1L,scale)
  stopifnot(sqrt_lo*sqrt_lo <= 5L, sqrt_hi*sqrt_hi > 5L)
  log2 <- log_interval(Q(2L), digits+4L)
  reduced_lo <- log_interval((5L+2L*sqrt_lo)/8L, digits+4L)
  reduced_hi <- log_interval((5L+2L*sqrt_hi)/8L, digits+4L)
  bounds <- outward(3L*log2$lo+reduced_lo$lo,
                    3L*log2$hi+reduced_hi$hi, digits)
  c(bounds, list(terms=list(log2=log2$terms, reduced_lower=reduced_lo$terms,
                           reduced_upper=reduced_hi$terms)))
}

minimum_M <- function(q, target, digits=80L) {
  stopifnot(q > 0L, q < 1L)
  logs <- log_interval(1L/(1L-q), digits+4L)
  D <- outward(logs$lo, logs$hi, digits)
  if (D$lo <= 0L) stop("Fixed precision cannot certify this case")
  candidate <- q_ceil(target$hi/D$lo)
  achieved_margin <- candidate*D$lo-target$hi
  previous_margin <- target$lo-(candidate-1L)*D$hi
  if (!(achieved_margin > 0L && previous_margin > 0L)) {
    stop("80-place bounds do not certify the minimum. Do not report success.")
  }
  list(M=candidate, D_lower=q_text(D$lo), D_upper=q_text(D$hi),
       D_series_terms=logs$terms, M_margin_lower=q_text(achieved_margin),
       previous_margin_lower=q_text(previous_margin), certified=TRUE)
}

beta_integer <- function(n, x) {
  result <- Q(0L)
  for (j in n:(2L*n-1L)) {
    result <- result + Q(chooseZ(2L*n-1L,j))*x^j*(1L-x)^(2L*n-1L-j)
  }
  result
}

spec <- jsonlite::read_json(file.path(root,"inputs","fixed_spec.json"), simplifyVector=FALSE)
expected <- jsonlite::read_json(file.path(root,"expected","exact_reference.json"), simplifyVector=FALSE)
alpha <- Q(spec$alpha); fixed_M <- as.integer(spec$fixed_M)
stopifnot(length(spec$cases)==12L, alpha==Q(4L,5L), fixed_M==10L)
expected_index <- new.env(parent=emptyenv())
for (row in expected$evaluations) {
  assign(paste(row$case_id,row$n,sep=":"),row,envir=expected_index)
}
target <- target_interval(alpha)
out_name <- if (mode=="--all") "all" else "representative"
out_dir <- file.path(root,"output",out_name)
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

rows <- list(); exact <- list(); enumerations <- list(); checked <- list()
endpoint_checks <- list(); row_id <- 0L
cases_to_run <- if (mode=="--all") spec$cases else
  Filter(function(x) x$id %in% c("C02","C08","C11"),spec$cases)

for (case in cases_to_run) {
  cid <- case$id; s <- as.integer(case$s); R <- Q(case$R); eps <- Q(case$epsilon)
  a <- R-1L; b <- 1L-1L/R; d <- (R-1L)^2L/R
  kstar <- eps/(d*(1L-eps)); if (kstar>Q(1L,2L)) kstar <- Q(1L,2L)
  n_rule <- min(s%/%2L,max(1L,as.integer(q_floor(kstar*s))))
  rule_names <- c("two_changed_goods","rounded_asymptotic_rule","all_goods_changed")
  rule_n <- c(1L,n_rule,s%/%2L)
  do_enum <- mode=="--all" && isTRUE(case$enumerate_all_n)
  ns <- if (do_enum) seq_len(s%/%2L) else sort(unique(rule_n))
  cache <- list()
  for (n in ns) {
    q <- q_positive(s,n,R,eps)
    p <- (1L-(1L-q)^fixed_M)^2L
    cert <- minimum_M(q,target)
    key <- paste(cid,n,sep=":")
    if (!exists(key,envir=expected_index,inherits=FALSE)) stop("Unfrozen evaluation: ",key)
    ref <- get(key,envir=expected_index,inherits=FALSE)
    stopifnot(q==Q(ref$q), z_text(cert$M)==ref$minimum_M)
    checked[[length(checked)+1L]] <- list(case_id=cid,n=n,q_exact_match=TRUE,
                                        M_exact_match=TRUE,minimum_certified=TRUE)
    ev <- list(case_id=cid,s=s,R=case$R,epsilon=case$epsilon,n=n,
               q=q_text(q),q_decimal_80=decimal_text(q),
               P_fixed_M=q_text(p),P_fixed_M_decimal_80=decimal_text(p),
               minimum_M=z_text(cert$M),minimum_T=z_text(2L*cert$M),
               certificate=cert[setdiff(names(cert),"M")])
    exact[[length(exact)+1L]] <- ev
    cache[[as.character(n)]] <- list(q=q,M=cert$M,ev=ev)
  }
  for (r in seq_along(rule_names)) {
    ev <- cache[[as.character(rule_n[[r]])]]$ev
    row_id <- row_id+1L
    rows[[row_id]] <- data.frame(case_id=cid,role=case$role,R=case$R,
      epsilon=case$epsilon,s=s,rule=rule_names[[r]],n=rule_n[[r]],
      changed_goods=2L*rule_n[[r]],q=ev$q_decimal_80,
      P_M10=ev$P_fixed_M_decimal_80,minimum_M=ev$minimum_M,
      minimum_T=ev$minimum_T,stringsAsFactors=FALSE)
  }
  q1 <- q_positive(s,1L,R,eps)
  n1_match <- q1 == (1L-eps/b)^(s-1L)/(R+1L)
  stopifnot(n1_match,s%%2L==0L)
  full_match <- q_positive(s,s%/%2L,R,eps)==beta_integer(s%/%2L,(b-eps)/(a+b))
  stopifnot(full_match)
  endpoint_checks[[length(endpoint_checks)+1L]] <- list(case_id=cid,
    n1_closed_form_exact=n1_match,no_neutral_goods_beta_cdf_exact=full_match)
  if (do_enum) {
    bestq <- cache[[1L]]$q; bestM <- cache[[1L]]$M
    for (v in cache) {
      if (v$q>bestq) bestq <- v$q
      if (v$M<bestM) bestM <- v$M
    }
    q_argmax <- ns[vapply(cache,function(v) isTRUE(v$q==bestq),logical(1))]
    M_argmin <- ns[vapply(cache,function(v) isTRUE(v$M==bestM),logical(1))]
    er <- Filter(function(x) x$case_id==cid,expected$enumerations)[[1L]]
    stopifnot(identical(as.integer(q_argmax),as.integer(unlist(er$q_argmax_n))),
              identical(as.integer(M_argmin),as.integer(unlist(er$M_argmin_n))))
    enumerations[[length(enumerations)+1L]] <- list(case_id=cid,
      q_argmax_n=as.list(q_argmax),M_argmin_n=as.list(M_argmin),
      minimum_M=z_text(bestM),rounded_rule_n=n_rule)
  }
  message(cid,": exact q, exact M, certificates and endpoints passed.")
}

q_double_sum <- function(s,n,R,eps) {
  a<-R-1L; b<-1L-1L/R
  if (eps>=b) return(Q(0L))
  A<-a+eps; B<-b-eps; m<-s-2L*n; total<-Q(0L)
  for (j in 0:(n-1L)) {
    for (ell in 0:(n-1L-j)) {
      if (m==0L && ell>0L) next
      coeff_l <- if (m==0L) Q(1L) else Q(chooseZ(m+ell-1L,ell))
      total <- total + Q(chooseZ(n+j-1L,j))*(A/(A+B))^j*
        coeff_l*(eps/b)^ell
    }
  }
  (B/(A+B))^n*(B/b)^m*total
}
small <- list()
for (pair in list(c("6/5","1/200"),c("3/2","1/20"))) {
  R<-Q(pair[[1L]]); eps<-Q(pair[[2L]])
  for (s in c(2L,4L,5L,6L)) for (n in seq_len(s%/%2L)) {
    stopifnot(q_positive(s,n,R,eps)==q_double_sum(s,n,R,eps))
    small[[length(small)+1L]] <- list(s=s,n=n,R=pair[[1L]],epsilon=pair[[2L]],
                                      direct_sum_exact_match=TRUE)
  }
}
stopifnot(q_positive(2L,1L,Q("6/5"),Q("1/6"))==0L)

table <- do.call(rbind,rows)
write.csv(table,file.path(out_dir,"comparison.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")
jsonlite::write_json(list(definitions=list(q="strict one-direction edge probability",
  probability="[1-(1-q)^M]^2",M="equal repeats at each budget",T="2M",
  fixed_M=fixed_M,alpha=spec$alpha),evaluations=exact,enumerations=enumerations),
  file.path(out_dir,"exact_results.json"),pretty=TRUE,auto_unbox=TRUE)
validation <- list(mode=mode,R=R.version.string,
  packages=list(gmp=as.character(packageVersion("gmp")),jsonlite=as.character(packageVersion("jsonlite"))),
  cases=length(cases_to_run),rule_rows=nrow(table),unique_evaluations=length(exact),
  enumerated_cases=length(enumerations),reference_matches=checked,endpoints=endpoint_checks,
  small_direct_sum_checks=small,strict_threshold_endpoint_zero=TRUE,
  rational_log_certificate_decimal_places=80L,
  candidate_method="ceil(C_upper/D_lower), using exact rationals; not floating logarithms",
  target_C_lower=q_text(target$lo),target_C_upper=q_text(target$hi),
  target_C_terms=target$terms,all_checks_passed=TRUE,random_draws=0L)
jsonlite::write_json(validation,file.path(out_dir,"validation.json"),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out_dir,"session_info.txt"))
message("Completed ",length(exact)," exact evaluations and ",nrow(table)," rule rows; all checks passed.")
