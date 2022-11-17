setwd(base_dir);
if(.Platform$OS.type == "windows") {
  # source("C:/Users/nikhi/Dropbox/research/earlier_code_files/FF_repro/FF_repro_functions.R");
} else {
  # Only source it once!
  if(!exists(as.character(substitute(R_PROFILE_SOURCED)))) {
    source("~/.Rprofile");
    setwd(base_dir);
  }
  # source("FF_repro_functions.R");
}
setwd(base_dir);
source(paste0(base_dir, "sorts.R"));









DOWNLOAD_OVERRIDE = F; # This will force data download irrespective of other factors
DOWNLOAD_NEW_DATA = F;

data_begin = "1962-01-01";
data_end = "2020-12-31";


# this option is for generating visually different lot colors in a greyscale image. Keep it FALSE unless you know what you are doing. Also the below is only tested/used for a few 
GREY_SCALE_PLOTS = F;



find_52_week_vars = T; # leave for some time


INCLUDE_DAVIS_BOOK_EQUITY_DATA = F;

EXCLUDE_NASDAQ = F;

EVALUATE__ROLLING_REGRESSIONS = F;

# Evaluation indicators (TRUE/FALSE)
EVALUATE__data.comp = T;
EVALUATE__data.crsp = T;
EVALUATE__data.crsp.cln = T;
EVALUATE__data.both = T;
EVALUATE__data.both.a = T;


EVALUATE__funda.anomaly = T; # fundamental anomalies using data.both.a. Anomalies are computed each month.


EVALUATE__tech.daily = F; # technical indicators from crps.dsf
EVALUATE__tech.monthly.anomaly = F; # technical anomalies using tech.daily. Daily signals are aggregated into monthly anomalies.


EVALUATE__tech2.daily = T; # technical indicators (alternative definition) using crsp.dsf
EVALUATE__tech2.monthly.anomaly = F; # technical anomalies using tech2.daily. Daily signals are aggregated into monthly anomalies.


EVALUATE__tech3.monthly = T; # evaluates technical indicators at month-end. Uses tech2.daily
EVALUATE__tech3.monthly.anomaly = F; # technical anomalies using tech3.monthly


EVALUATE__tech.patterns.monthly = F; # technical patterns (will take days)
EVALUATE__tech2.monthly.anomaly = F; # technical anomalies based on tech.patterns.monthly






# default decile and tercile names
DECILE_NAMES = c(paste0("DEC_", 1:9), "DEC_Ten");
TERCILE_NAMES = paste0("TERC_", 1:3);
QUINTILE_NAMES = paste0("QUNT_", 1:5);







coalesce = data.table::fcoalesce;

# this function should be called for each cusip separately
na_locf_until = function(x, n) {
  # in time series data, fill in na's untill indicated n
  l <- cumsum(!is.na(x)); # Basically tells which are NA and which aren't
  c(NA, x[!is.na(x)])[replace(l, ave(l, l, FUN=seq_along) > (n+1), 0) + 1]
}








# read function with windows and linux directories
ind_readRDS = function(win_file = NULL, linux_file = NULL, refhook = NULL) {
  if(is.null(win_file) & is.null(linux_file)) {
    stop("Provide atleast one valid location in win_file or linux_file");
  }
  
  if(.Platform$OS.type == "windows" & !is.null(win_file)) {
    return(readRDS(win_file, refhook = refhook));
  } else if(.Platform$OS.type == "unix" & !is.null(linux_file)) {
    return(readRDS(linux_file, refhook = refhook));
  } else {
    stop("Please provide a consistent file (w.r.t. the underlying OS)");
  }
  
}




# text evaluation
EVL = function(text, print = F) {
  if(print) {
    cat("Executing --> ", text, "\n");
  }
  eval(parse(text = text), envir=.GlobalEnv);
}


# to make query building cleaner
filter_query_by_date = function(query, date_var) {
  
  # if there is a 'where' already in the query then no need to use the first and
  if(EXISTS(data_begin)) {
    query = paste0(query,
                   ifelse(stringr::str_detect(query, " where "), " and ", " where "),
                   date_var, " >= '", data_begin, "'");
  }
  
  # if there is a 'where' already in the query then no need to use the first and
  if(EXISTS(data_end)) {
    query = paste0(query,
                   ifelse(stringr::str_detect(query, " where "), " and ", " where "),
                   date_var, " <= '", data_end, "'");
  }
  
  return(query);
}



# list all the functions in a R source file
is_function = function (expr) {

  is_assign = is.call(expr) && as.character(expr[[1]]) %in% c('=', '<-', 'assign');
  
  if (!is_assign) {
    return(FALSE);
  }
  value = expr[[3]];
  is.call(value) && as.character(value[[1]]) == 'function';
}

list_functions <- function(filename = paste0(base_dir, "signal_dispersion_ESSENTIALS.R"),
                           sort = F) {
  
  function_name = function(expr) as.character(expr[[2]]);
  
  file_parsed = parse(filename);
  functions = Filter(is_function, file_parsed);
  function_names = unlist(Map(function_name, functions));
  if(sort == T) {
    return(sort(function_names));
  } else {
    return(function_names);
  }
}





get_chunks = function(n, c) {
  1:n %>% split(., ceiling(c * seq_along(.) / (length(.)))) %>% return;
}


# Run until pass (in try catch)

TRY_CATCH = gtools::defmacro(expression, print.attempts = T, max.attempts = 10, ret_val = NULL, sleep = 0, expr = {
  passed__ = F;
  ret__ = ret_val;
  attempt__ = 0;
  while(passed__ == F & attempt__ < max.attempts) {
    Sys.sleep(sleep*1e-3); # in milli seconds
    attempt__ = attempt__ + 1;
    if(print.attempts == T) {
      cat("Attempt: ", attempt__, "\n");
    }
    tryCatch({
      ret__ = eval(expression);
      passed__ = T;
    }, error = function(e){
    }, finally = {})
  }
  remove(passed__, attempt__);
  return(ret__);
})


# proper way to call exists
EXISTS = function(x, ...) {
  exists(as.character(substitute(x)), ...);
}


# The below is a macro (short-hand) which will read the file from disk in case it is not already loaded
GET = gtools::defmacro(dt, name = NULL, expr = {
  dt_name__ = as.character(substitute(dt));
  if(!exists(dt_name__)) {
    file_name__ = ifelse(is.null(name), paste0(dt_name__, ".RData"), name);
    dt = readRDS(file_name__);
  }
})



# The below is a macro (short-hand) for the most used case of merge
MRG = gtools::defmacro(dt_1, dt_2, by = c("cusip", "Date"), expr = {
  cols_1 = names(dt_1);
  cols_2 = names(dt_2);
  by_cols = by;
  mrg_cols = cols_2[which(!(cols_2 %in% cols_1))];
  if(len(mrg_cols) == 0) {
    # no need to merge
    dt_1 = dt_1;
  } else {
    mrg_cols = c(by_cols, mrg_cols);
    dt_1 = merge(dt_1, dt_2[, mrg_cols, with = F], by = by, all.x = T);
  }
})



# This function reads a variable from file: f and return it. If the file doesn't exists then it returns 0
read_or_0 = function(f) {
  if(file.exists(f)) {
    return( readRDS(f) );
  } else {
    return(0);
  }
}

# the function lag is present both in dplyr and data.table. We need to use lag from dplyr
LAG = function(x, n = 1L, default = NA, order_by = NULL, ...) {
  dplyr::lag(x, n, default, order_by, ...);
}

# the function lag is present both in dplyr and data.table. We need to use lag from dplyr
LEAD = function(x, n = 1L, default = NA, order_by = NULL, ...) {
  dplyr::lead(x, n, default, order_by, ...);
}

DELTA = function(x) {
  x - dplyr::lag(x);
}

# AVG(x) is same as (x_t + x_{t-1}) / 2
AVG = function(x, times = 2, lag = 0) {
  shift(x, seq(lag, lag + times - 1), type = "lag") %>% Reduce(`+`, .) / times;
}

# AVG1(x) is same as (x_{t-1} + x_{t-2}) / 2
AVG1 = function(x, times = 2) {
  AVG(x, times = times, lag = 1);
}


# Cumulative mean/sum etc with NAs
cumsum = function(x, na.rm = F) {
  if(na.rm == F) {
    base::cumsum(x);
  } else {
    "[<-"(x, !is.na(x), base::cumsum(na.omit(x)));
  }
}

cumprod = function(x, na.rm = F) {
  if(na.rm == F) {
    base::cumprod(x);
  } else {
    "[<-"(x, !is.na(x), base::cumprod(na.omit(x)));
  }
}

cummin = function(x, na.rm = F) {
  if(na.rm == F) {
    base::cummin(x);
  } else {
    "[<-"(x, !is.na(x), base::cummin(na.omit(x)));
  }
}

cummax = function(x, na.rm = F) {
  if(na.rm == F) {
    base::cummax(x);
  } else {
    "[<-"(x, !is.na(x), base::cummax(na.omit(x)));
  }
}

cummean = function(x, na.rm = F) {
  if(na.rm == F) {
    dplyr::cummean(x);
  } else {
    "[<-"(x, !is.na(x), dplyr::cummean(na.omit(x)));
  }
}






# CL0(x,y) is a short-hand for coalesce(x,0) + coalesce(y,0) with the additional requirement that both x and y are not NA.
CL0 = function(...) {
  x = 0;
  for(i in list(...)) {
    x = x + coalesce(i, as(0, typeof(i)));
  }
  
  # make sure all are not NAs
  y = 0;
  for(i in list(...)) {
    y = y | !is.na(i);
  }
  y[y == 0] = NA;
  
  x = x*y;
  
  return(x);
}


# CL0_prod(x,y) is a short-hand for coalesce(x,1) * coalesce(y,1) with the additional requirement that both x and y are not NA.
CL0_prod = function(...) {
  x = 1;
  for(i in list(...)) {
    x = x * coalesce(i, as(1, typeof(i)));
  }
  
  # make sure all are not NAs
  y = 0;
  for(i in list(...)) {
    y = y | !is.na(i);
  }
  y[y == 0] = NA;
  
  x = x*y;
  
  return(x);
}


COR = function(x, y = NULL, use = "pairwise.complete.obs", method = c("pearson", "kendall", "spearman"), round_digits = 2, ts_avg = NULL, ts_exclude_NA = T) {
  
  if(is.null(ts_avg)) {
    C = cor(x, y, use, method);
  } else {
    if(is.na(match(ts_avg, names(x)))) {
      stop("Please make sure ", ts_avg, " is part of your data.frame!")
    }
    if(!is.null(y)) {
      stop("Please provide entire data.frame as single argument!");
    }
    # find C for each ts_avg then average it
    cols = setdiff(names(x), ts_avg);
    unq_ts = unique(x[, ts_avg, with = F]) %>% unlist %>% sort %>% na.omit;
    C = matrix(0, nrow = len(cols), ncol = len(cols));
    C_zeros = C + diag(len(cols));
    rownames(C) = colnames(C) = cols;
    for(u in unq_ts) {
      C__ = cor(x[get(ts_avg) == u, cols, with = F], y, use, method);
      if(ts_exclude_NA == T) {
        C__ = coalesce(C__, C_zeros);
      }
      C = C + C__;
    }
    C = C / len(unq_ts);
  }
  
  if(is.null(round_digits) || round_digits == 0) {
    # leave C as it is
  } else {
    C = round(C, round_digits);
  }
  
  return(C);
}





COR_2 = function(x, y = NULL, method = c("pearson", "spearman"), round_digits = 2) {
  
  if(is.null(y)) {
    C = Hmisc::rcorr(x %>% as.matrix, type = method);
  } else {
    C = Hmisc::rcorr(x %>% as.matrix, y %>% as.matrix, type = method);
  }
  
  if(is.null(round_digits) || round_digits == 0) {
    # leave C as it is
  } else {
    C$r = round(C$r, round_digits);
  }
  
  return(C);
}


# Cross-sectional correlations. Takes time when both is set to T.
cross_COR = function(dt, cols = NULL, both = F, stars = NULL, round_digits = 3, use_rank_vars = F) {
  
  if(is.null(cols)) {
    cols = names(dt);
  }
  
  N = length(cols);
  
  unq_Dates = dt[, unique(Date)] %>% sort;
  
  C = array(NA, dim = c(length(unq_Dates), N, N));
  
  for(i in 1:length(unq_Dates)) {
    C[i,,] = dt[Date == unq_Dates[i], cols, with = F] %>% COR(method = "pearson", round_digits = NULL);
  }
  
  # basically need to find both pearson and spearman correlation. The former shall go in lower triangle.
  if(both == T) {
    
    C2 = array(NA, dim = c(length(unq_Dates), N, N));
    
    # takes time
    for(i in 1:length(unq_Dates)) {
      if(use_rank_vars == T) {
        C2[i,,] = dt[Date == unq_Dates[i], paste0(cols, "_Date_RANK"), with = F] %>% COR(method = "pearson", round_digits = NULL);
      } else {
        C2[i,,] = dt[Date == unq_Dates[i], cols, with = F] %>% COR(method = "spearman", round_digits = NULL);
      }
    }
    
    # Need to take lower.tri of C and upper.tri of C2 and add them!
    lw = lower.tri(matrix(N*N, N, N), diag = F);
    up = upper.tri(matrix(N*N, N, N), diag = F);
    
    for(i in 1:length(unq_Dates)) {
      
      C[i,,] = lw*C[i,,] + up*C2[i,,];
      
      diag(C[i,,]) = NA;
      
    }
    
  }
  
  
  corr = apply(C, c(2,3), mean, na.rm = T);
  corr[is.nan(corr)] = NA;
  corr = round(corr, round_digits);
  
  sd = apply(C, c(2,3), sd, na.rm = T);
  n = apply(!is.na(C), c(2,3), sum);
  se = sd / sqrt(n);
  abs.t.stat = abs(corr) / se;
  p.val = 2*pt(abs.t.stat, Inf, lower.tail = F);
  
  if(is.null(stars)) {
    # stars = stars(p.val); # see the function stars()
    stars = ifelse(p.val > .05, "#", " ");
  }
  
  final = matrix(paste0(corr, stars), N, N);
  diag(final) = NA;
  
  return(final);
  
}



# both correlations combined
COR_both = function(dt, cols = NULL, round_digits = 2) {
  
  if(is.null(cols)) {
    cols = names(dt);
  }
  
  N = length(cols);
  
  C1 = COR(dt, method = "pearson", round_digits = round_digits);
  C2 = COR(dt, method = "spearman", round_digits = round_digits);
  
  # Need to take lower.tri of C and upper.tri of C2 and add them!
  lw = lower.tri(matrix(N*N, N, N), diag = F);
  up = upper.tri(matrix(N*N, N, N), diag = F);
  
  C = lw*C1 + up*C2;
  diag(C) = NA;
  
  return(C);

}




COV = function(x, y = NULL, use = "pairwise.complete.obs", method = c("pearson", "kendall", "spearman")) {
  cov(x, y, use, method)
}

VAR = function(x, y = NULL, use = "pairwise.complete.obs", ...) {
  var(x, y, use, method, ...)
}


is.POS = function(x) {
  return(x > 0)
}

is.NEG = function(x) {
  return(x < 0)
}

is.ZERO = function(x) {
  return(x == 0)
}


ann_ret = function(x, FUN = mean, rnd_digits = 2, shift = 0, scale = 1200, ...) {
  round(scale*FUN(x, ...) + shift, rnd_digits)
}



# The below function will assign lags of specified vars into data (in place. No need to copy!)
assign_lags = function(data, vars, lags = 1, group_by = "cusip", time_var = "Date") {
  setorderv(data, c(group_by, time_var));
  for(v in vars) {
    for(l in lags) {
      var_name = paste0("LAG_", l, "_", v);
      if(is.element(var_name, colnames(data))) {
        data[, eval(var_name) := NULL];
      }
      data[, eval(var_name) := LAG(get(v), l), by = get(group_by)];
    }
  }
}


# The below function will assign dispersion measures of all the splits in the specified data.table in place. No need to copy!
assign_dispersion_measures = function(data, splits = c("50_50", "80_20", "70_30"),
                                      disp_measures = c("std_dev", "abs_dev", "entropy", "hhi"),
                                      keywords = NULL, nomenclature = "", min_data_rows = 10) {
  
  if(is.null(keywords)) {
    keywords = "signal";
  }
  
  if(nomenclature != "") {
    nomenclature = paste0("_", nomenclature);
  }
  
  which_cols = colnames(data); # start with all columns and then successively apply filter rules
  for(k in keywords) {
    which_cols = grep(k, which_cols, value = T);
  }
  
  new_cols = NULL;

  
  for(s in splits) {
    
    suffix = paste0(nomenclature, if(s != "") "_", s);

    diff_sig = paste0("num_diff", suffix); # difference between number of buys and sells
    abs_diff_sig = paste0("abs_num_diff", suffix); # difference between number of buys and sells
    mean_sig = paste0("mean_sig", suffix); # `mean_signal` for NO split; `mean_signal_70_30` for 70_30 split AND `mean_signal_funda_70_30` if nomenclature is `funda`
    
    data[, eval(diff_sig) := NA_real_];
    data[, eval(abs_diff_sig) := NA_real_];
    data[, eval(mean_sig) := NA_real_];

    new_cols = c(new_cols, diff_sig, abs_diff_sig, mean_sig);
    
    data[ get(paste0("num_tot", suffix)) >= min_data_rows, eval(diff_sig) := get(paste0("num_buy", suffix)) - get(paste0("num_sell", suffix))];
    data[, eval(abs_diff_sig) := abs(get(diff_sig))];
    data[, eval(mean_sig) := get(diff_sig) / get(paste0("num_tot", suffix))];

        
    cols = grep(s, which_cols, value = T);

    for(d in disp_measures) {
      
      dev_sig = paste0(d, suffix);
      data[, eval(dev_sig) := NA_real_];
      new_cols = c(new_cols, dev_sig);

      data[, n_tot := get(paste0("num_tot", suffix))];
      
      if(d == "std_dev") {
        data[, sq_sum := rowSums(.SD * .SD, na.rm = T), .SDcols = cols];
        
        data[ n_tot >= min_data_rows, eval(dev_sig) := sqrt( (sq_sum - n_tot*get(mean_sig)^2) / (n_tot - 1) )];
        
        data[, sq_sum := NULL];
      }
      
      if(d == "abs_dev") {
        data[n_tot >= min_data_rows, eval(dev_sig) := rowMeans(abs(get(mean_sig) - .SD), na.rm = T), .SDcols = cols];
      }
      
      if(d == "entropy") {
        
        data[, pr_min_1 := get(paste0("num_sell", suffix)) / n_tot];
        data[, pr_0     := get(paste0("num_hold", suffix)) / n_tot];
        data[, pr_1     := get(paste0("num_buy", suffix))  / n_tot];
        
        data[n_tot >= min_data_rows,
             eval(dev_sig) := ifelse(pr_min_1 > 0, -pr_min_1*log2(pr_min_1), 0) +
               ifelse(pr_0 > 0, -pr_0*log2(pr_0), 0) +
               ifelse(pr_1 > 0, -pr_1*log2(pr_1), 0) ];
        
        data[, `:=`(pr_min_1 = NULL, pr_0 = NULL, pr_1 = NULL)];
      }
      
      if(d == "hhi") {
        
        data[, pr_min_1 := get(paste0("num_sell", suffix)) / n_tot];
        data[, pr_0     := get(paste0("num_hold", suffix)) / n_tot];
        data[, pr_1     := get(paste0("num_buy", suffix))  / n_tot];
        
        data[n_tot >= min_data_rows,
             eval(dev_sig) := pr_min_1^2 + pr_0^2 + pr_1^2];
        
        data[, `:=`(pr_min_1 = NULL, pr_0 = NULL, pr_1 = NULL)];
      }
      
      data[, n_tot := NULL];
      
    }
  }
  
  # Just making sure that the new columns appear in the end  
  setcolorder(data, c(setdiff(colnames(data), new_cols), new_cols));
  
}



# The below function assigns buy/sell/hold/NA counts for each split. No need to copy; in place assignment.
assign_counts = function(data, splits = c("50_50", "80_20", "70_30"), keywords = NULL, nomenclature = "") {
  
  if(is.null(keywords)) {
    keywords = "signal";
  }
  
  if(nomenclature != "") {
    nomenclature = paste0("_", nomenclature);
  }
  
  which_cols = colnames(data); # start with all columns and then successively apply filter rules
  for(k in keywords) {
    which_cols = grep(k, which_cols, value = T);
  }
  
  new_cols = NULL;
  
  for(s in splits) {
    
    cols = grep(s, which_cols, value = T);
    
    vars = c("count_NA", "num_buy", "num_sell", "num_hold") %>% paste0(., nomenclature);
    funs = c(is.na, is.POS, is.NEG, is.ZERO);
    
    for(v in 1:length(vars)) {
      
      var_name = paste0(vars[v], if(s != "") "_", s);
      
      data[, eval(var_name) := rowSums( funs[[v]](.SD), na.rm = T ), .SDcols = cols];

      new_cols = c(new_cols, var_name);
      
    }
    
    suffix = paste0(nomenclature, if(s != "") "_", s);
    
    var_name = paste0("num_tot", suffix);
    
    data[, eval(var_name) := get(paste0("num_buy", suffix)) + 
                             get(paste0("num_sell", suffix)) +
                             get(paste0("num_hold", suffix))];
    
    new_cols = c(new_cols, var_name);
    
  }

  # Just making sure that the new columns appear in the end  
  setcolorder(data, c(setdiff(colnames(data), new_cols), new_cols));
  
}


# Remove trailing split (6 chars) from a list of character names
remove_trailing_splits = function(cols, splits = "") {
  if(splits == "") {
    splits = "50_50$|70_30$|80_20$";
  } else {
    splits = paste0(splits, "$");
  }
  
  idx = grep(splits, cols);
  if(length(idx) > 0) {
    cols[idx] = cols[idx] %>% substr(., 1, nchar(.)-6);
  }
  
  return(cols);
}




# calculate dependent variable(s)
calc_dep_vars = function(dt, vars = c("turn"), output = c("delta", "detrended"), log = T) {
  
  # Example
  # setorder(dt, cusip, Date);
  # dt[, log_turn := log(turn)];
  # dt[is.na(log_turn) | is.nan(log_turn) | is.infinite(log_turn), log_turn := NA_real_];
  # dt[, delta_log_turn := DELTA(log_turn), by = cusip];
  # dt[, log_turn_detrended := NA_real_];
  # valid_cusips = dt[, sum(!is.na(log_turn)), by = cusip][V1 > 0]$cusip;
  # dt[cusip %in% valid_cusips, log_turn_detrended := log_turn - lm(log_turn ~ Date)$coef %*% rbind(rep(1, .N), Date), by = cusip];
  
  if( match(c("cusip", "Date"), colnames(dt)) %>% is.na %>% sum > 0 ) {
    stop("The data.table doesn't have Date and cusip vars!")
  }
  
  setorder(dt, cusip, Date);
  
  for(v in vars) {
    
    var_name = v;
    if(log == T) {
      var_name = paste0("log_", v);
      dt[get(v) > 0, eval(var_name) := log( get(v) )]; # log is only defined for +ve valued vars
    } else
    
    # sanitize
    dt[is.na(get(var_name)) | is.nan(get(var_name)) | is.infinite(get(var_name)), eval(var_name) := NA_real_];
    
    for(ot in output) {
      
      if(ot == "delta") {
        dt[, eval(paste0(ot, "_", var_name)) := DELTA(get(var_name)), by = cusip];
      }
      
      if(ot == "detrended") {
        dt[, eval(paste0(var_name, "_", ot)) := NA_real_];
        valid_cusips = dt[, sum(!is.na(get(var_name))), by = cusip][V1 > 0]$cusip;
        # below takes 30 sec. Reason to not take residuals directly is that if NA is present, length of residuals vector will be less.
        dt[cusip %in% valid_cusips, eval(paste0(var_name, "_", ot)) := get(var_name) - lm(get(var_name) ~ Date)$coef %*% rbind(rep(1, .N), Date), by = cusip];
      }
      
    }
    
  }
  
}






optimal_layout = function(w, h, N, alpha = 0.5) {
  
  # forces alhpah to be between 0 and 1
  alpha = max(0, min(alpha, 1));
  
  best_obj = Inf;
  best_n1 = NA_integer_;
  best_n2 = NA_integer_;
  
  for(n1 in 1:N) {
    for(n2 in 1:N) {
      
      if(n1 * n2 < N) {
        next;
      }
      
      obj = alpha*((n1/n2)/(h/w) - 1)^2 + (1-alpha)*(n1*n2/N - 1)^2;
      
      if(obj < best_obj) {
        best_obj = obj;
        best_n1 = n1;
        best_n2 = n2;
      }
      
    }
  }
  
  return(c(best_n1, best_n2));
  
}





# useful for plotting huge data columns
PLT = function(x, y = NULL, layout = NULL, layout_tuning = 0.5, gui_plot = F, width = 7.17, height = 10, pdf_file = "del.pdf",
               grid = F, pts = NULL, figure_margin = NULL, title = NULL, title_size = 2,
               h_line = NA, v_line = NA, h_line_col = "grey25", v_line_col = "grey75",
               type = "l", lwd = 2, col = "blue",  xlim = NULL, ylim = NULL,
               pre_expr = NA, post_expr = NA, pre_loop_expr = NULL, post_loop_expr = NULL,
               log = "", main = NULL, sub = NULL, xlab = NULL, ylab = NULL, lty = NULL,
               ann = par("ann"), axes = TRUE, frame.plot = axes,
               panel.first = NULL, panel.last = NULL, asp = NA, ...) {
  
  if(is.data.frame(x) == F) {
    if(is.null(y)) {
      # only x provided. Assume the other as counting from 1 to nrow(x)
      x = data.frame(x = 1:length(x), y = x);
    } else {
      x = data.frame(x = x, y = y);
    }
    print(x);
  }

  if(is.null(pts)) {
    pts = ifelse(is.data.frame(x), nrow(x), length(x));
  } else if(pts == "all") {
    pts = min(ifelse(is.data.frame(x), nrow(x), length(x)), 1000);
  } else if(is.numeric(pts) | is.integer(pts)) {
    pts = pts;
  } else {
    stop("pts must either be: 1) number/integer, 2) all OR 3) NULL")
  }
  
  pts = min(ifelse(is.data.frame(x), nrow(x), length(x)), ifelse(is.null(pts), 1000, pts));
  
  probabilities = seq(from = 0, to = 1, length.out = pts);
  
  is_DT = is.data.table(x);
  

  if(is.data.frame(x)) {
    
    if(ncol(x) < 2) {
      stop("data frame must have at-least 2 columns");
    }
    
    setDF(x);

    xx = quantile(x[,1], probs = probabilities, na.rm = T, names = F, type = 1);
    
    yy = x[match(xx, x[,1] %>% unlist), 2:ncol(x)];
    
    if(is.data.frame(yy) == F) {
      yy = data.frame(yy);
    }
    
    if(is_DT) {
      setDT(x);
    }
    
    if(is.null(xlab)) {
      xlab =  colnames(x[0,1]);
    } else if(length(xlab) == 1) {
      xlab = rep(xlab, ncol(x)-1);
    }
    
    if(is.null(ylab)) {
      ylab =  colnames(x[0, 2:ncol(x)]);
    } else if(length(ylab) == 1) {
      ylab = rep(ylab, ncol(x)-1);
    }
    
    # vectorize the below args
    if(!is.null(main) & length(main) == 1) { main = rep(main, ncol(x) - 1) }
    if(!is.null(sub) & length(sub) == 1) { sub = rep(sub, ncol(x) - 1) }
    if(!is.null(log) & length(log) == 1) { log = rep(log, ncol(x) - 1) }
    if(!is.null(type) & length(type) == 1) { type = rep(type, ncol(x) - 1) }
    if(!is.null(lty) & length(lty) == 1) { lty = rep(lty, ncol(x) - 1) }
    if(!is.null(col) & length(col) == 1) { col = rep(col, ncol(x) - 1) }
    if(!is.null(lwd) & length(lwd) == 1) { lwd = rep(lwd, ncol(x) - 1) }
    
    if(length(h_line) == 1) {
      h_line = rep(h_line, ncol(x) - 1);
    }
    
    if(length(v_line) == 1) {
      v_line = rep(v_line, ncol(x) - 1);
    }
    
    if(length(pre_expr) == 1) {
      pre_expr = rep(pre_expr, ncol(x) - 1);
    }
    
    if(length(post_expr) == 1) {
      post_expr = rep(post_expr, ncol(x) - 1);
    }
    
  }
  
  if(is.null(pdf_file)) {
    pdf_file = paste0(Sys.time() %>% gsub(" |:", "-", .), ".pdf");
  }
  
  
  if(gui_plot == F) {
    pdf(pdf_file, width, height);
  }
  
  # multi_plot(xx, yy, layout, type, xlim, ylim, log, main, sub, xlab, ylab, ann, axes, frame.plot, panel.first, panel.last, asp, ...);
  
  org_mfrow = par("mfrow");
  org_mar = par("mar");
  
  if(!is.null(figure_margin)) {
    par(mar = figure_margin);
  }
  
  if(is.null(layout)) {
    par(mfrow = optimal_layout(width, height, length(yy), layout_tuning)); # set the layout to be optimal. It mayn't be perfect
  } else {
    if(layout[1]*layout[2] < length(yy)) {
      stop("Wrong layout provided!");
    }
    par(mfrow = layout);
  }
  
  if(length(grid) == 1) {
    grid = rep(grid, length(yy));
  }

  if(!is.null(pre_loop_expr)) {
    for(p in 1:length(pre_loop_expr)) {
      eval(parse(text = pre_loop_expr[p]))
    }
  }
  
  for(i in 1:length(yy)) {

    if(!is.na(pre_expr[[i]][1])) {
      for(p in 1:length(pre_expr[[i]])) {
        eval(parse(text = pre_expr[[i]][p]))
      }
    }
    
    plot(xx, yy[, i] %>% unlist,
         type = type[i], xlim = xlim, ylim = ylim, log = log[i], main = main[i], sub = sub[i],
         xlab = xlab[i], ylab = ylab[i], lwd = lwd[i], col = col[i], lty = lty[i],
         ann = ann, axes = axes, frame.plot = frame.plot, panel.first = panel.first, panel.last = panel.last, asp = asp, ...);
    
    if(!is.na(post_expr[[i]][1])) {
      for(p in 1:length(post_expr[[i]])) {
        eval(parse(text = post_expr[[i]][p]))
      }
    }
    
    # h_line and v_line are list of vectors. One list element for each line respectively.
    if(!is.na(h_line[[i]][1])) {
      abline(h = h_line[[i]], lwd = 2, col = h_line_col, lty = 4);
    }
    
    if(!is.na(v_line[[i]][1])) {
      abline(v = v_line[[i]], lwd = 2, col = v_line_col, lty = 4);
    }
    
    if(grid[i] == T) {
      grid(NULL, NULL, lty = 1, col = "cornsilk2");
    }
  }
  
  if(!is.null(post_loop_expr)) {
    for(p in 1:length(post_loop_expr)) {
      eval(parse(text = post_loop_expr[p]))
    }
  }
  
  if(!is.null(title)) {
    title(title, line = -2, outer = T, cex.main = title_size);
  }
  
  if(gui_plot == F) {
    dev.off();
  }
  
  par(mar = org_mar);
  par(mfrow = org_mfrow);
  
}





# epsilon range
EPS = function(x = NULL, n = 1) {
  if(is.null(x)) {
    return( .Machine$double.eps * n );
  } else {
    return(c(x - .Machine$double.eps*n, x + .Machine$double.eps*n));
  }
}





# Winsorize "cols" within data.table dt at the level of "level". Do it separately for each every_col variable. If exclude == T, then extreme data points are removed from the table rather than winsoring it
# NOTE: since dt is being modified here, call by reference will not work! So the best way to use winsorization would be to apply it on a small set of variables before a regression.
# In case the winsorization needs to be done on entire panel use every_col = "ones" where ones is a column of all 1's
winsorize = function(dt, cols, level = 0.01, every_col = "Date", exclude = F, quantile_type = 7) {
  
  if(is.na(match(every_col, colnames(dt)))) {
    stop(paste0("The data.table doesn't have ", every_col, " and cusip vars!"));
  }
  
  pct_lo = dt[, lapply(.SD, quantile, p = c(level), na.rm = T, type = quantile_type), .SDcols = cols, by = every_col];
  pct_hi = dt[, lapply(.SD, quantile, p = c(1-level), na.rm = T, type = quantile_type), .SDcols = cols, by = every_col];
  
  # Rename column names for differentiation
  colnames(pct_lo) = paste0("pct_lo.", colnames(pct_lo));
  colnames(pct_hi) = paste0("pct_hi.", colnames(pct_hi));
  # Don't change name of every_col variable
  setnames(pct_lo, paste0("pct_lo.", every_col), every_col);
  setnames(pct_hi, paste0("pct_hi.", every_col), every_col);
  
  dt = merge(dt, pct_lo, by = every_col, all.x = T, allow.cartesian = T);
  dt = merge(dt, pct_hi, by = every_col, all.x = T, allow.cartesian = T);
  
  for(c in cols) {
    col_class = dt[, class(get(c))];
    c_lo = paste0("pct_lo.", c);
    c_hi = paste0("pct_hi.", c);
    if(exclude == T) {
      # exclude at "level"
      # setting NA is equivalent to dropping it (they won't be used in rgeressions)
      dt = dt[get(c) < get(c_lo), eval(c) := as(NA, col_class)];
      dt = dt[get(c) > get(c_hi), eval(c) := as(NA, col_class)];
    } else {
      # winsorize at "level"
      dt[get(c) < get(c_lo), eval(c) := as(get(c_lo), col_class)];
      dt[get(c) > get(c_hi), eval(c) := as(get(c_hi), col_class)];
    }
  }
  
  # remove pct_lo and pct_hi vars
  del_cols = c(colnames(pct_lo), colnames(pct_hi)) %>% setdiff(every_col);
  dt[, (del_cols) := NULL];

  remove(pct_lo, pct_hi, c_lo, c_hi);
  
  return(dt);
}






# Finding ranks and then normalizing it to fall between 0 and 1
find_ranks = function(dt, vars, by_col = NULL) {
  
  options(warn = 2); # Need errors so that TRY_CATCH can return NA
  
  for(v in vars) {
    
    new_var_name = paste0(v, "_", paste0(by_col, collapse = "_"), if(!is.null(by_col)) "_", "RANK");
    # Do nothing if the variable already exists
    if(is.element(names(dt), new_var_name) %>% sum > 0) {
      cat("variable: ", new_var_name, " already exists!", "\n", sep = "");
      next;
    }
    
    print(v);
    
    # ranks start from 1
    
    dt[, new_var__ := frank(get(v), na.last = "keep", ties.method = "dense"), by = by_col];
    
    dt[, temp_min := TRY_CATCH(min(new_var__, na.rm = T), F, 1, NA_integer_), by = by_col];
    dt[, temp_max := TRY_CATCH(max(new_var__, na.rm = T), F, 1, NA_integer_), by = by_col];
    
    # normalize to lie betweeon 0 and 1    
    dt[, new_var__ := (new_var__ - temp_min) / (temp_max - temp_min)];
    
    dt[, temp_min := NULL];
    dt[, temp_max := NULL];
    
    setnames(dt, "new_var__", paste0(v, "_", paste0(by_col, collapse = "_"), if(!is.null(by_col)) "_", "RANK"));
    
  }
  
  options(warn = 0);
  
}





# Once ranks are known finding quintiles is extremely easy and fast
find_quintiles_from_ranks = function(dt, vars, boundaries = 10*(1:9), var_suffix = "_decile",
                                     quintile_names = NULL, by_col = NULL, print = T, custom_rank_vars = NULL) {
  
  pct = c(0, boundaries, 100) / 100;
  
  n = length(pct) - 1; # total number of quintiles
  
  if(is.null(quintile_names)) {
    quintile_names = paste0("QUINTILE_", 1:n);
  }
  
  cnt = 0;
  for(v in vars) {
    cnt = cnt + 1;
    
    if(is.null(custom_rank_vars)) {
      rank_var = paste0(v, "_", paste0(by_col, collapse = "_"), if(!is.null(by_col)) "_", "RANK");
    } else {
      rank_var = custom_rank_vars[cnt];
    }
    
    
    if(!is.element(rank_var, colnames(dt))) {
      stop(paste0("Var: ", rank_var, " not present in data.table! ",
                  "Run find_ranks(dt, \"", v, "\", ", ifelse(is.null(by_col), "NULL", paste0("\"", by_col, "\"")), ") first!"));
    }
    
    if(print) {
      print(v);
    }
    
    final_var = paste0(v, var_suffix);
    for(i in 1:n) {
      dt[pct[i] <= get(rank_var) & get(rank_var) < pct[i+1], eval(final_var) := quintile_names[i]];
    }
    
  }
  
}







# PLEASE consider using find_ranks() followed by find_quintiles_from_ranks() for speed!
# Find quintiles at boundaries for variables: vars. There will be length(boundaries) + 1 quintiles. Default is equally placed deciles
find_quintiles = function(dt, vars, boundaries = 10*(1:9), var_suffix = "_decile", quintile_names = NULL, by_col = NULL, print = T) {
  
  pct = c(0, boundaries, 100);
  
  n = length(pct) - 1; # total number of quintiles
  
  if(is.null(quintile_names)) {
    quintile_names = paste0("QUINTILE_", 1:n);
  }
  
  for(d in vars) {
    
    if(print) {
      print(d);
    }
    
    dt[, eval(paste0(d, ".pct.", 0)) := -Inf];
    for(b in boundaries) {
      var_name = paste0(d, ".pct.", b);
      if(is.null(by_col)) {
        dt[, eval(var_name) := quantile(get(d), p = b/100, na.rm = T)];
      } else {
        dt[, eval(var_name) := quantile(get(d), p = b/100, na.rm = T), by = by_col];
      }
    }
    dt[, eval(paste0(d, ".pct.", 100)) := Inf];
    
    final_var = paste0(d, var_suffix);
    for(i in 1:n) {
      lt_pct = paste0(d, ".pct.", pct[i]);
      rt_pct = paste0(d, ".pct.", pct[i+1]);
      dt[get(d) >= get(lt_pct) & get(d) < get(rt_pct), eval(final_var) := quintile_names[i]];
    }
    
    for(i in 1:length(pct)) {
      dt[, eval(paste0(d, ".pct.", pct[i])) := NULL];
    }
  }
  
}








# Rolling Fama-French type regressions. It is assummed that all return data is available
# The below finds alpha and various betas. I do not store residuals. Although it can be found by retadj_rf - alpha - beta*mkt_rf
rolling_regression = function(dt, y_col, x_cols, print = F, min_obs = 24, window_len = 60, col_names = NULL) {
  
  k = length(x_cols) + 1; # num of beta estimates (i.e. regression parameters)
  
  if(is.null(col_names)) {
    col_names = c(paste("b", 0:(k-1)), "ivol", "R2");
  } else {
    if(length(col_names) != (k+2)) {
      stop("You must provide the right number of columns. It is num of parameters + 2.");
    }
  }
  
  formula = paste0(y_col, " ~ ", paste0(x_cols, collapse = " + "));
  
  is_present = is.element(col_names, colnames(dt));
  # only need to be added if all columns are not already present
  if(sum(!is_present) != 0) {
    dt[, (col_names[!is_present]) := NA_real_]; # initialize to NA
  }
  
  setorder(dt, cusip, Date);
  unq_cusip = dt[, unique(cusip)];
  data_lens = dt[, .N, by = cusip][, N];
  
  st_time = Sys.time();
  
  for(c in seq_along(unq_cusip)) {
    
    dt_c = copy( dt[cusip == unq_cusip[c]] );
    
    N = data_lens[c];
    
    # reg_result = matrix(NA_real_, nrow = N, ncol = k + 2) %>% as.data.frame;
    reg_result = dt_c[, col_names, with = F];
    setDF(reg_result);
    colnames(reg_result) = col_names;
    
    if(N >= window_len) {
      
      for(i in window_len:N) {
        # only if no new column was added, can we hope for an already existent data on alpha on beta
        if(sum(!is_present) == 0) {
          # check whether alpha and beta are already calculated and present. This can be achieved by checking whether these values are non-NA before each regression. This will reduce the incremental burden. Below already_present == 0 represents that all values are non-NA.
          already_present = dt_c[i, col_names, with = F] %>% unlist %>% is.na %>% sum;
          if(already_present == 0) {
            next;
          }
        }
        
        idx = (i-(window_len-1)):i;
        
        dt_c_non_NA = na.omit(dt_c[idx, c(x_cols, y_col), with = F]);
        if(nrow(dt_c_non_NA) < min_obs) {
          # Not enough data
          if(FALSE & print) {
            cat(c, " out of total ", length(unq_cusip), "  --  ", formula,
                " -- Not enough observations to estimate. Need ", min_obs,
                " but present only ", nrow(dt_c_non_NA), "\n");
          }
          next;
        }
        
        f = lm(formula, data = dt_c[idx]) %>% summary;
        
        reg_result[i, 1:k]  = f$coefficients[, "Estimate"];
        reg_result[i, k+1]  = f$residuals %>% sd;
        reg_result[i, k+2]  = f$r.squared;
        
      }
      
    }
    
    dt[cusip == unq_cusip[c], (col_names) := reg_result];
    
    curr_time = Sys.time();
    time_used = as.numeric(curr_time) - as.numeric(st_time);
    frac_done = (sum(data_lens[1:c]) - (window_len-1)*c) / (nrow(dt) - (window_len-1)*length(unq_cusip));
    est_completion_time = time_used / frac_done;
    time_left = est_completion_time - time_used;
    time_left = chron::times(time_left/(24*60*60)) %>% as.character;
    
    if(print) {
      cat(c, " out of total ", length(unq_cusip), "  --  ", formula, " -- Time left: ", time_left, "\n");
    }

  }
  
}






# wrapper function for stargazer
REG = function(..., type = "text", out = NULL, omit_ser = T,
               title = "", style = "default", summary = NULL, out.header = FALSE, column.labels = NULL, column.separate = NULL,
               covariate.labels = NULL, dep.var.caption = NULL, dep.var.labels = NULL, dep.var.labels.include = TRUE, align = FALSE,
               coef = NULL, se = NULL, t = NULL, p = NULL, t.auto = TRUE, p.auto = TRUE, ci = FALSE, ci.custom = NULL, ci.level = 0.95,
               ci.separator = NULL, add.lines = NULL, apply.coef = NULL, apply.se = NULL, apply.t = NULL, apply.p = NULL, apply.ci = NULL,
               colnames = NULL, column.sep.width = "5pt", decimal.mark = NULL, df = TRUE, digit.separate = NULL, digit.separator = NULL, 
               digits = NULL, digits.extra = NULL, flip = FALSE, float = TRUE, float.env = "table", font.size = NULL, header = TRUE, 
               initial.zero = NULL, intercept.bottom = TRUE, intercept.top = FALSE, keep = NULL, keep.stat = NULL, label = "", model.names = NULL, 
               model.numbers = NULL, multicolumn = TRUE, no.space = NULL, notes = NULL, notes.align = NULL, notes.append = TRUE, notes.label = NULL, 
               object.names = FALSE, omit = NULL, omit.labels = NULL, omit.stat = NULL, omit.summary.stat = NULL, omit.table.layout = NULL, omit.yes.no = c("Yes", "No"),
               order = NULL, ord.intercepts = FALSE, perl = FALSE, report = NULL, rownames = NULL, rq.se = "nid", selection.equation = FALSE, 
               single.row = FALSE, star.char = NULL, star.cutoffs = NULL, suppress.errors = FALSE, table.layout = NULL, table.placement = "!htbp", 
               zero.component = FALSE, summary.logical = TRUE, summary.stat = NULL, nobs = TRUE, mean.sd = TRUE, min.max = TRUE, median = FALSE, iqr = FALSE) {
  
  
  
  if(is.null(out)) {
    if(type == "text") {
      # ok! file name not mandatory for text printing
      if(is.null(omit.stat) & omit_ser) {
        omit.stat = "ser";
      }
    } else {
      stop("Please provide output file name!")
    }
  } else {
    # Do both html and latex
    out_html  = paste0(out, ".html");
    out_latex = paste0(out, ".tex");
  }
  
  
  final = list();
  for(f in list(...)) {
    if(class(f) != "list") {
      f = list(f);
    }
    f = f[which(!sapply(f, is.null))];
    final = append(final, f);
  }
  
  
  if(type == "both") {
    
    stargazer::stargazer(final, type = "html", out = out_html,
                         title = title, style = style, summary = summary, out.header = out.header, column.labels = column.labels, column.separate = column.separate,
                         covariate.labels = covariate.labels, dep.var.caption = dep.var.caption, dep.var.labels = dep.var.labels, dep.var.labels.include = dep.var.labels.include, align = align,
                         coef = coef, se = se, t = t, p = p, t.auto = t.auto, p.auto = p.auto, ci = ci, ci.custom = ci.custom, ci.level = ci.level,
                         ci.separator = ci.separator, add.lines = add.lines, apply.coef = apply.coef, apply.se = apply.se, apply.t = apply.t, apply.p = apply.p, apply.ci = apply.ci,
                         colnames = colnames, column.sep.width = column.sep.width, decimal.mark = decimal.mark, df = df, digit.separate = digit.separate, digit.separator = digit.separator,
                         digits = digits, digits.extra = digits.extra, flip = flip, float = float, float.env = float.env, font.size = font.size, header = header,
                         initial.zero = initial.zero, intercept.bottom = intercept.bottom, intercept.top = intercept.top, keep = keep, keep.stat = keep.stat, label = label, model.names = model.names,
                         model.numbers = model.numbers, multicolumn = multicolumn, no.space = no.space, notes = notes, notes.align = notes.align, notes.append = notes.append, notes.label = notes.label,
                         object.names = object.names, omit = omit, omit.labels = omit.labels, omit.stat = omit.stat, omit.summary.stat = omit.summary.stat, omit.table.layout = omit.table.layout, omit.yes.no = omit.yes.no,
                         order = order, ord.intercepts = ord.intercepts, perl = perl, report = report, rownames = rownames, rq.se = rq.se, selection.equation = selection.equation,
                         single.row = single.row, star.char = star.char, star.cutoffs = star.cutoffs, suppress.errors = suppress.errors, table.layout = table.layout, table.placement = table.placement,
                         zero.component = zero.component, summary.logical = summary.logical, summary.stat = summary.stat, nobs = nobs, mean.sd = mean.sd, min.max = min.max, median = median, iqr = iqr);

    
    stargazer::stargazer(final, type = "latex", out = out_latex,
                         title = title, style = style, summary = summary, out.header = out.header, column.labels = column.labels, column.separate = column.separate,
                         covariate.labels = covariate.labels, dep.var.caption = dep.var.caption, dep.var.labels = dep.var.labels, dep.var.labels.include = dep.var.labels.include, align = align,
                         coef = coef, se = se, t = t, p = p, t.auto = t.auto, p.auto = p.auto, ci = ci, ci.custom = ci.custom, ci.level = ci.level,
                         ci.separator = ci.separator, add.lines = add.lines, apply.coef = apply.coef, apply.se = apply.se, apply.t = apply.t, apply.p = apply.p, apply.ci = apply.ci,
                         colnames = colnames, column.sep.width = column.sep.width, decimal.mark = decimal.mark, df = df, digit.separate = digit.separate, digit.separator = digit.separator,
                         digits = digits, digits.extra = digits.extra, flip = flip, float = float, float.env = float.env, font.size = font.size, header = header,
                         initial.zero = initial.zero, intercept.bottom = intercept.bottom, intercept.top = intercept.top, keep = keep, keep.stat = keep.stat, label = label, model.names = model.names,
                         model.numbers = model.numbers, multicolumn = multicolumn, no.space = no.space, notes = notes, notes.align = notes.align, notes.append = notes.append, notes.label = notes.label,
                         object.names = object.names, omit = omit, omit.labels = omit.labels, omit.stat = omit.stat, omit.summary.stat = omit.summary.stat, omit.table.layout = omit.table.layout, omit.yes.no = omit.yes.no,
                         order = order, ord.intercepts = ord.intercepts, perl = perl, report = report, rownames = rownames, rq.se = rq.se, selection.equation = selection.equation,
                         single.row = single.row, star.char = star.char, star.cutoffs = star.cutoffs, suppress.errors = suppress.errors, table.layout = table.layout, table.placement = table.placement,
                         zero.component = zero.component, summary.logical = summary.logical, summary.stat = summary.stat, nobs = nobs, mean.sd = mean.sd, min.max = min.max, median = median, iqr = iqr);
    
  } else {
    
    stargazer::stargazer(final, type = type, out = out,
                         title = title, style = style, summary = summary, out.header = out.header, column.labels = column.labels, column.separate = column.separate,
                         covariate.labels = covariate.labels, dep.var.caption = dep.var.caption, dep.var.labels = dep.var.labels, dep.var.labels.include = dep.var.labels.include, align = align,
                         coef = coef, se = se, t = t, p = p, t.auto = t.auto, p.auto = p.auto, ci = ci, ci.custom = ci.custom, ci.level = ci.level,
                         ci.separator = ci.separator, add.lines = add.lines, apply.coef = apply.coef, apply.se = apply.se, apply.t = apply.t, apply.p = apply.p, apply.ci = apply.ci,
                         colnames = colnames, column.sep.width = column.sep.width, decimal.mark = decimal.mark, df = df, digit.separate = digit.separate, digit.separator = digit.separator,
                         digits = digits, digits.extra = digits.extra, flip = flip, float = float, float.env = float.env, font.size = font.size, header = header,
                         initial.zero = initial.zero, intercept.bottom = intercept.bottom, intercept.top = intercept.top, keep = keep, keep.stat = keep.stat, label = label, model.names = model.names,
                         model.numbers = model.numbers, multicolumn = multicolumn, no.space = no.space, notes = notes, notes.align = notes.align, notes.append = notes.append, notes.label = notes.label,
                         object.names = object.names, omit = omit, omit.labels = omit.labels, omit.stat = omit.stat, omit.summary.stat = omit.summary.stat, omit.table.layout = omit.table.layout, omit.yes.no = omit.yes.no,
                         order = order, ord.intercepts = ord.intercepts, perl = perl, report = report, rownames = rownames, rq.se = rq.se, selection.equation = selection.equation,
                         single.row = single.row, star.char = star.char, star.cutoffs = star.cutoffs, suppress.errors = suppress.errors, table.layout = table.layout, table.placement = table.placement,
                         zero.component = zero.component, summary.logical = summary.logical, summary.stat = summary.stat, nobs = nobs, mean.sd = mean.sd, min.max = min.max, median = median, iqr = iqr);
    
  }
  
  
  
}




# wrapper function for grep
GRP = function(x, pattern, ignore.case = F, perl = F, value = T, fixed = F, useBytes = F, invert = F, non_char = F) {
  if(!non_char & !("character" %in% class(x))) {
    stop("You're trying to grep a non-character object. To proceed, please set non_char to TRUE.");
  }
  
  return( grep(pattern = pattern, x = x, ignore.case = ignore.case, perl = perl, value = value, fixed = fixed, useBytes = useBytes, invert = invert) );
}





# takes around 2 mins to run for each variable
GRT_adjustment = function(dt, vars, min_obs = 20, print = T) {
  
  setorder(dt, cusip, Date);
  dt[, counts__ := 1:.N];
  
  for(v in vars) {
    
    st_time = Sys.time();
    cat("GRT adjustment started for var: ", v, "\n");
    
    new_var = paste0(v, "_GRT");
    
    dt[, eval(new_var) := NA_real_];
    
    if(is.element(new_var, colnames(dt)) == FALSE) {
      stop("New column not added in data.table");
    }
    
    unq_cusip = dt[, unique(cusip)] %>% sort;
    n = length(unq_cusip);
    
    # this loop takes some time (20+ mins for each v)
    cnt = 0;
    for(c in unq_cusip) {
      
      cnt = cnt + 1;
      
      dt_c = copy( dt[cusip == c, c("Date", v, "counts__"), with = F] );
      dt_c = dt_c %>% na.omit;
      
      names(dt_c)[2] = "grt_var";
      
      pre_N = nrow(dt_c);
      # residuals will be all 0 if there are very few data points in the time series (of each firm)
      if(pre_N < min_obs) {
        next;
      }
      
      pre_mu = dt_c[, mean(grt_var, na.rm = T)];
      pre_sigma = dt_c[, sd(grt_var, na.rm = T)];
      
      dt_c[, Date2 := Date^2];
      dt_c[, month := month(Date)];
      
      # mean equation
      fit = TRY_CATCH( lm(grt_var ~ Date + Date2 + factor(month), data = dt_c),
                       print.attempts = F, max.attempts = 1, ret_val = NULL);
      
      if(is.null(fit)) {
        next;
      }
      
      dt_c[, first_stage := log(fit$res^2)];
      
      # variance equation
      fit_2 = TRY_CATCH( lm(first_stage ~ Date + Date2 + factor(month), data = dt_c),
                         print.attempts = F, max.attempts = 1, ret_val = NULL);
      
      if(is.null(fit_2)) {
        next;
      }
      
      # because I am regressing log(u^2) in fit_2. fit_2$res are errors of dep. var. log(fit$res^2). Therefore the errors of fit$res are exp(fit_2$res/2)
      second_stage = exp(fit_2$res/2);
      second_stage = second_stage - mean(second_stage, na.rm = T);
      second_stage = second_stage / sd(second_stage, na.rm = T);

      # dt[cusip == c & !is.na(get(v)), eval(new_var) := second_stage*pre_sigma + pre_mu];
      dt[dt_c$counts__, eval(new_var) := second_stage*pre_sigma + pre_mu];
      
      if(print & (cnt %% round(n/100) == 0)) {
        
        curr_time = Sys.time();
        time_used = as.numeric(curr_time) - as.numeric(st_time);
        pcnt_done = cnt / round(n/100);
        est_completion_time = time_used / (pcnt_done / 100);
        time_left = est_completion_time - time_used;
        
        cat(pcnt_done, "% of the process is complete. Est. Time left: ", chron::times(time_left/(24*60*60)) %>% as.character, "\n");
      }
      
    }
    
  }
  
  dt[, counts__ := NULL];
  
  return(dt);
  
}





# 2-way dependent sort
two_way_dependent_sort = function(data, num_ports = c(2,2), col_1, col_2, new_two_way_var, col_1_id = "LEFT", col_2_id = "RIGHT") {
  
  data[, (new_two_way_var) := NA_character_];
  
  temp = form_breakpoints(data, char_list = c(col_1, col_2),
                          brk_pts = list(seq(0, 100, length.out = num_ports[1] + 1)[2:num_ports[1]],
                                         seq(0, 100, length.out = num_ports[2] + 1)[2:num_ports[2]]),
                          sort_order = c(1,2), is_perc = T, infinite_extremes = T, return_indices = T, summary_vector = "retadj");
  
  for(i in 1:num_ports[1]) {
    for(j in 1:num_ports[2]) {
      
      data[ temp$indices[[i]][[j]], eval(new_two_way_var) := paste(col_1_id, i, col_2_id, j, sep = "_")];
      
    }
  }
  
}





# Make a decile view with col_1 going top to bottom and col_2 going left to right. fun is applied on the var
decile_view = function(dt, col_1, col_2, grid = c(10, 10), small_calc = T, var = "log_turn", fun = mean, print = F, dim_names = T, switch_dimensions = F, ...) {
  
  if(switch_dimensions) {
    temp__ = col_1;
    col_1 = col_2;
    col_2 = temp__;
    remove(temp__);
  }
  
  var_names =  var %>%
    strsplit("+", fixed = T) %>% unlist %>% trimws %>%
    strsplit("-", fixed = T) %>% unlist %>% trimws %>%
    strsplit("*", fixed = T) %>% unlist %>% trimws %>%
    strsplit("/", fixed = T) %>% unlist %>% trimws %>%
    strsplit("(", fixed = T) %>% unlist %>% trimws %>%
    strsplit(")", fixed = T) %>% unlist %>% trimws %>%
    strsplit("\\^\\d+", fixed = F) %>% unlist %>% trimws;
  var_names = stringr::str_replace_all(var_names, "^\\d+$", "");
  var_names = var_names[nchar(var_names) > 0];
  
  temp = dt[, lapply(.SD, fun, ...), .SDcols = var_names, by = c(col_1, col_2)] %>% copy;
  temp = temp[!is.na(get(col_1)) & !is.na(get(col_2)) & get(col_1) != "" & get(col_2) != ""];

  eval(parse(text = paste0("temp[, final := ", var, "]")));

  # temp = temp[, c(col_1, col_2, "final"), with = F];
  temp[, (var_names) := NULL];
  
  temp = temp %>% dcast(get(col_1) ~ get(col_2), value.var = "final");
  row_names = temp[,1] %>% unlist;
  temp = temp[, -1]; # remove the names column
  
  if(is.null(grid)) {
    grid = dim(temp);
  }
  
  setDF(temp);
  max_diff = (temp[grid[1], grid[2]] - temp[1,1]) %>% unlist;
  if(print) {
    cat("Top-Bottom: ", col_1, " and Left-Right: ", col_2, " -- max. diff: ", max_diff, "\n");
  }
  
  if(small_calc == T) {
    rownames(temp) = row_names;
    
    if(dim_names == F) {
      colnames(temp) = rownames(temp) = NULL;
    }
    
    return(temp);
  }
  
  temp[, paste0(grid[2], "_minus_", 1)]  = temp[, grid[2]] - temp[, 1];
  temp[, paste0(grid[2], "_by_", 1)]     = temp[, grid[2]] / temp[, 1];
  temp[grid[1]+1, ]                      = temp[grid[1], ] - temp[1, ]; # [paste0(grid[1], "_minus_", 1), 1]
  temp[grid[1]+2, ]                      = temp[grid[1], ] / temp[1, ];
  setDT(temp);
  
  temp = cbind(c(row_names, paste0(grid[1], c("_minus_", "_by_"), 1)), temp);
  colnames(temp)[1] = "";
  
  if(dim_names == F) {
    colnames(temp) = rownames(temp) = NULL;
  }
  
  return(temp);
}



# to extract proj. R2 from lfe::felm object
myfun = function(f) { summary(f)["P.adj.r.squared"] }




# latex lag vars
latex_lag_vars = function(cols, include = NULL, exclude = NULL) {
  
  if(!is.null(include) & !is.null(exclude)) {
    stop("Only include or exclude should be included!");
  }
  
  ret_val = c();
  
  if(!is.null(include) & length(include) > 0) {
    
    for(i in 1:len(cols)) {
      if(i %in% include) {
        ret_val = c(ret_val, paste0("$", cols[i], "_{t-1}$"));
      } else {
        ret_val = c(ret_val, paste0("$", cols[i], "$"));
      }
    }
    
  } else if(!is.null(exclude) & length(exclude) > 0) {
    
    for(i in 1:len(cols)) {
      if(i %in% exclude) {
        ret_val = c(ret_val, paste0("$", cols[i], "$"));
      } else {
        ret_val = c(ret_val, paste0("$", cols[i], "_{t-1}$"));
      }
    }
    
  } else {
    
    ret_val = paste0("$", cols, "_{t-1}$");
    
  }
  
  return(ret_val);
}

  





# kable and kableExtra fromatting needed for regression tables in Rmarkdown
kable_format = function(f, exp_vars = NULL, is_summary = T, excl_prcnt_exp = F, excl_proj_r2 = F, excl_adj_r2 = F, excl_r2 = T, excl_obs = F, extra_text = NULL, compact = F, col_names_offset = 0, excl_exp_vars = NULL, sig_var = "std_err", force_math_mode = T, force_tnote = F, force_tabularx = F, force_tabularx_r2 = T, exp_var_perm = NULL) {
  
  # f provides the list of regression models or their summaries generated by lfe_summ() function. Provide summary to get the most out of kable_format() function. I have developed and improved it using the summary models
  # exp_vars gives the names of explanatory variables. If missing R variable names from the regression model are used
  # is_summary is set to TRUE if summary is provided in f. If regression model is specified then is_summary must be set to FALSE. is_summary should be set to TRUE unless you know what you are doing!
  # excl_prcnt_exp if TRUE will exclude precentage R^2 explained
  # excl_proj_r2 if TRUE will exclude projected R^2 from the bottom of the table
  # excl_adj_r2 if TRUE will exclude adjusted R^2 from the bottom of the table
  # excl_r2 if TRUE will exclude R^2 from the bottom of the table
  # excl_obs if TRUE will exclude number of obs from the bottom of the table
  # extra_text is used to add extra lines in the bottom of the table (after coefficients and before R^2). It must be a list of vectors with each element of the list for each additional row to be added. The first element of vector would be the identifying text (such as 'Fixed Effects?'). The rest of the vector would be values to be displayed (like "Yes", "No", "N.A." etc).
  # compact if FALSE would add two blank lines before proj.R2 begins. Don't set to TRYE unless required!
  # col_names_offset is the amount by which regression numbers are offset. The numbers start from `1 + col_names_offset`
  # Provide the explanatory variable names that should be omitted from the table. The names here must match the R variable names as specified in the regression model R code
  # sig_var tells which significance variable (out of pval, tstat and std_err) to use.
  # force_math_mode is TRUE sets the coefficient and significances to latex math mode. If FALSE, these values are printed in text mode.
  # exp_var_perm changes the order of appearance of explanatory variables. If NULL, default order is used.
  
  if(is_summary == F & excl_prcnt_exp == F) {
    stop("prcnt_exp can be included (i.e. not excluded) only when providing summaries");
  }
  
  n = sum(!sapply(f, is.null));
  vars = c();
  
  # Just store the summary of f in s
  if(is_summary) {
    s = f;
  } else {
    s = vector("list", n);
    for(i in 1:n) {
      s[[i]] = summary(f[[i]]);
    }
  }
    
  for(i in 1:n) {
    # remove unwanted explanatory variables (if any)
    s[[i]]$coefficients = as.data.frame(s[[i]]$coefficients);
    if(!is.null(excl_exp_vars)) {
      rem_idx = which(rownames(s[[i]]$coefficients) %in% excl_exp_vars);
      if(len(rem_idx) > 0) {
        s[[i]]$coefficients = s[[i]]$coefficients[-rem_idx,];
      }
    }
    
    vars = union(vars, rownames(s[[i]]$coefficients));
  }
  
  if(is.null(exp_vars)) {
    exp_vars = vars;
  }
  
  # permute explanatory variables
  if(!is.null(exp_var_perm)) {
    
    exp_var_perm = unique(exp_var_perm);
    
    if(len(exp_var_perm) != len(vars)) {
      # Stop Invalid permutation
      stop("Provide a valid permutation of these ", len(vars), " vars: ", paste0(vars, collapse = ", "));
    }
    
    if(class(exp_var_perm) == "character" & !setequal(exp_var_perm, vars)) {
      # entries must match all entries in vars
      stop("Provide a valid permutation of these ", len(vars), " vars: ", paste0(vars, collapse = ", "));
    }
    
    # set permutation
    if(class(exp_var_perm) == "character") {
      exp_var_perm = match(exp_var_perm, vars);
    }
    vars = vars[exp_var_perm];
    exp_vars = exp_vars[exp_var_perm];
    
  }
  
  df = data.frame(vars = c(rbind(vars, ""),
                           c(if(!excl_r2) "r2",
                             if(!excl_proj_r2) "proj_r2",
                             if(!excl_adj_r2) "adj_r2",
                             if(!excl_prcnt_exp) "prcnt_explained",
                             if(!excl_obs) "n_obs")
                           ), stringsAsFactors = F);

  df[, paste0("reg_", 1:n)] = NA;
  
  for(i in 1:n) {
    # idx = match(df$vars, rownames(s$coefficients)) %>% {which(!is.na(.))};
    idx = match(rownames(s[[i]]$coefficients), df$vars); # these indices are where df should be filled with reg coef.
    
    p = s[[i]]$coefficients[, 4];
    
    if(force_tnote == T) {
      significance = stars(p, p.val.cutoff = c(0.1, 0.05, 0.01), symbol = "*", pad = "", tnote = force_tnote);
    } else {
      significance = stars(p, p.val.cutoff = c(0.1, 0.05, 0.01), symbol = "*", pad = " ");
    }
    
    if(force_math_mode == T) {
      sprintf_str = "$%0.3f$";
      phantom_str = "$%+0.3f$";
    } else {
      sprintf_str = "%0.3f";
      phantom_str = "%+0.3f";
    }
    
    if(force_tabularx == T) {
      phantom_str = sprintf_str;
    }
    
    coef = paste0(sprintf(phantom_str, s[[i]]$coefficients[, 1]), significance);
    std_err = sprintf(sprintf_str, s[[i]]$coefficients[, 2]);
    tstat = sprintf(phantom_str, s[[i]]$coefficients[, 3]);
    pval = sprintf(sprintf_str, s[[i]]$coefficients[, 4]);
    
    if(force_tabularx == F) {
      # put a phantom for + signs
      coef = gsub("+", "\\phantom{-}", coef, fixed = T);
      tstat = gsub("+", "\\phantom{-}", tstat, fixed = T);
      # put a phantom for positive quantitites
      std_err = paste0("\\phantom{-}", std_err);
      pval = paste0("\\phantom{-}", pval);
    }
    
    # enclose values inside brackets
    tstat = paste0("(", tstat, ")");
    std_err = paste0("(", std_err, ")");
    pval = paste0("(", pval, ")");
    
    df[idx, paste0("reg_", i)] = coef;
    if(sig_var == "pval") {
      df[idx+1, paste0("reg_", i)] = pval;
    } else if (sig_var == "tstat") {
      df[idx+1, paste0("reg_", i)] = tstat;
    } else {
      # default
      df[idx+1, paste0("reg_", i)] = std_err;
    }
    
    if(!excl_prcnt_exp) {
      df[df$vars == "prcnt_explained", i+1] = ifelse(is.na(s[[i]]$prcnt_explained), "", sprintf(phantom_str, s[[i]]$prcnt_explained));
      if(force_tabularx == F) {
        df[df$vars %in% "prcnt_explained", i+1] = 
          gsub("+", "\\phantom{-}", df[df$vars %in% "prcnt_explained", i+1], fixed = T);
      }
      # add surrounding {} for tabularx
      if(force_tabularx_r2 == T) {
        df[df$vars %in% "prcnt_explained", i+1] = paste0("{", df[df$vars %in% "prcnt_explained", i+1], "}");
      }
    }
    
    if(!excl_r2) {
      df[df$vars == "r2", i+1] = sprintf(phantom_str, s[[i]]$r.squared);
      if(force_tabularx == F) {
        df[df$vars %in% "r2", i+1] = 
          gsub("+", "\\phantom{-}", df[df$vars %in% "r2", i+1], fixed = T);
      }
      # add surrounding {} for tabularx
      if(force_tabularx_r2 == T) {
        df[df$vars %in% "r2", i+1] = paste0("{", df[df$vars %in% "r2", i+1], "}");
      }
    }
    
    if(!excl_proj_r2) {
      df[df$vars == "proj_r2", i+1] = sprintf(phantom_str, s[[i]]$P.adj.r.squared);
      if(force_tabularx == F) {
        df[df$vars %in% "proj_r2", i+1] = 
          gsub("+", "\\phantom{-}", df[df$vars %in% "proj_r2", i+1], fixed = T);
      }
      # add surrounding {} for tabularx
      if(force_tabularx_r2 == T) {
        df[df$vars %in% "proj_r2", i+1] = paste0("{", df[df$vars %in% "proj_r2", i+1], "}");
      }
    }
    
    if(!excl_adj_r2) {
      df[df$vars == "adj_r2", i+1] = sprintf(phantom_str, s[[i]]$adj.r.squared);
      if(force_tabularx == F) {
        df[df$vars %in% "adj_r2", i+1] = 
          gsub("+", "\\phantom{-}", df[df$vars %in% "adj_r2", i+1], fixed = T);
      }
      # add surrounding {} for tabularx
      if(force_tabularx_r2 == T) {
        df[df$vars %in% "adj_r2", i+1] = paste0("{", df[df$vars %in% "adj_r2", i+1], "}");
      }
    }
    
    if(!excl_obs) {
      df[df$vars == "n_obs", i+1] = prettyNum(s[[i]]$N, big.mark = ",", scientific = F);
      if(force_tabularx == F) {
        df[df$vars == "n_obs", i+1] = paste0("\\phantom{-}", df[df$vars == "n_obs", i+1]);
      }
      # add surrounding {} for tabularx
      if(force_tabularx_r2 == T) {
        df[df$vars %in% "n_obs", i+1] = paste0("{", df[df$vars %in% "n_obs", i+1], "}");
      }
    }
    
  }
  
  # extra_text is a list of vectors. The first element of vector would be the identifying text (such as 'Fixed Effects?'). The rest of the vector would be values to be displayed (like "Yes", "No", "N.A." etc).
  if(!is.null(extra_text)) {
    for(j in 1:len(extra_text)) {
      df = rbind(df, extra_text[[j]])
    }
  }
  
  
  # The below is just for formatting purposes
  
  # Add two blank lines before proj. R2 begins. Don't do this if compact == T
  df = rbind(df[1:(2*length(exp_vars)),],
             if(!compact) NA,
             if(!compact) NA,
             df[(2*length(exp_vars)+1):nrow(df),]);
  rownames(df) = 1:nrow(df);
  
  # Basically we need line separation for several (exp. var, std err) pairs. We don't need this for last exp. var.
  last = 5 - excl_obs - excl_r2 - excl_adj_r2 - excl_proj_r2 - excl_prcnt_exp; # 5 (one each) for r2, proj_r2, adj_r2, prcnt_exp and observation.
  extra_lines = ifelse(is.null(extra_text), 0, length(extra_text));
  line_sep = c(rep(c("", "\\addlinespace"), length(exp_vars) - 1), rep("", ifelse(compact, 2, 2 + 2) + last + extra_lines)); # 2 for last exp. var., 2 for the lines separating exp. vars and R2 vars, `last` for R2, num_obs etc and lastly `extra_lines` for additional text lines. Skip the 2 separating lines if compact == T
  
  row_names = c( rbind(exp_vars, ""),
                 if(!compact) c("", ""),
                 c(if(!excl_r2) "$R^2$",
                   if(!excl_proj_r2) "$\\textrm{Within} \\: R^2$",
                   if(!excl_adj_r2) "$\\textrm{Adj.} \\: R^2$",
                   if(!excl_prcnt_exp)"$\\% \\: R^2 \\: \\textrm{Explained}$",
                   if(!excl_obs) "$\\textrm{Observations}$")
                );
  
  if(!is.null(extra_text)) {
    row_names = c(row_names, sapply(extra_text, dplyr::first));
  }
  
  # shift oolnames by a phantom character as well
  if(force_tabularx == F) {
    col_names = c("", paste0("\\phantom{-}(", 1:n + col_names_offset, ")"));
  } else {
    col_names = c("", paste0("{(", 1:n + col_names_offset, ")}"));
  }
  
  old_exp_vars = df$vars %>% GRP("^$", invert = T);
  
  df$vars = row_names;
  colnames(df) = col_names;
  
  df[is.na(df)] = "";
  
  n_hline = nrow(df) - (6 - excl_obs - excl_r2 - excl_adj_r2 - excl_proj_r2 - excl_prcnt_exp - compact + extra_lines); # lines up from below where a horizontal line is to be added
  
  
  return(list(df = df, n = n, old_exp_vars = old_exp_vars, col_names = col_names, line_sep = line_sep, n_hline = n_hline));
  
}





kable_CAT = function(obj, file, col_sep_pts = NULL, latex_line_before = NULL, latex_line_after = NULL) {
  
  obj = as.character(obj);
  
  if(!is.null(latex_line_before)) {
    obj = paste0("\n\n", latex_line_before, "\n\n", obj);
  }
  
  # default spacing is 6pt on either side
  if(!is.null(col_sep_pts)) {
    obj = paste0("\n\n", "\\setlength{\\tabcolsep}{", col_sep_pts, "pt}", "\n\n", obj);
  }
  
  if(!is.null(latex_line_after)) {
    obj = paste0(obj, "\n\n", latex_line_after, "\n\n");
  }
  
  cat(obj, file = file);
  saveRDS(obj, paste0(file, ".RData"));
}





# performs univariate sorts on variables in cols. It can be decile or tercile. name of decile vars can be overwritten by names. Summary is provided for summ_var. All numeric outcomes are rounded. Column names can be modified to be more accurate, however it will become difficult to use them in expressions (thus use it only for viewing/printing).
# Output is a host of matrices of size: length(cols) times 10 (or 3). The one to be used directly in documents is print and print_sc which gives avg change in means across deciles along with statistical significance.
univariate_sorts = function(dt, cols, type = "decile", names = NULL, summ_var = "log_turn", round_digits = 3, modify_colnames = T, factor = 1) {
  
  # _sc stands for successive decile changes
  
  float_print = paste0("%0.", round_digits, "f");
  
  if(length(factor) == 1) {
    factor = rep(factor, length(cols));
  }
  
  if(type == "decile") {
    n = 10;
    idx_sc = c(1, 1:9);
    name__ = "D";
  } else {
    n = 3;
    idx_sc = c(1, 1:2);
    name__ = "T";
  }
  
  short_name = substr(cols, 7, nchar(cols) - (1 + nchar(type))); # extra 1 is for underscore in _decile / _tercile
  
  if(is.null(names)) {
    names = short_name;
  }
  
  stat_cols = c("mean_diff", "sum_se", "t.stat", "dof", "pval", "stars", "final", "print");
  stat_cols = c(c("mean", "sd", "n"), stat_cols, paste0(stat_cols, "_sc"));
  
  port_sort = array(NA_real_, dim = c(length(cols), n, length(stat_cols)), dimnames = list(names, paste0(name__, 1:n), stat_cols));
  
  dim = dim(port_sort);
  
  for(i in 1:length(cols)) {
    port_sort[i,,"mean"] = factor[i] * ( decile_view(dt, "ones", cols[i], NULL, T, summ_var, mean, print = F, na.rm = T) %>% unlist );
    port_sort[i,,"sd"]   = factor[i] * ( decile_view(dt, "ones", cols[i], NULL, T, summ_var, sd, print = F, na.rm = T) %>% unlist );
    port_sort[i,,"n"]    = decile_view(dt, "ones", cols[i], NULL, T, summ_var, length, print = F) %>% unlist;
  }
  
  port_sort = lapply(1:dim[3], function(i) port_sort[,,i]);
  names(port_sort) = stat_cols;
  
  # See http://www.stat.yale.edu/Courses/1997-98/101/meancomp.htm for help
  
  # Difference in means, sum of std errors, t.stats, dof and pval
  port_sort$mean_diff = port_sort$mean - port_sort$mean[,1]; # decile_d minus decile_1
  port_sort$sum_se = sqrt( (port_sort$sd^2 / port_sort$n) + (port_sort$sd[,1]^2 / port_sort$n[,1]) ); # decile_d plus decile_1
  port_sort$t.stat = port_sort$mean_diff / port_sort$sum_se;
  port_sort$dof = pmin(port_sort$n - 1, port_sort$n[,1] - 1);
  port_sort$pval = 2*pt(-abs( port_sort$t.stat ), port_sort$dof);
  
  # now do the same for successive decile changes
  port_sort$mean_diff_sc = port_sort$mean - port_sort$mean[,idx_sc]; # decile_d minus decile_1
  port_sort$sum_se_sc = sqrt( (port_sort$sd^2 / port_sort$n) + (port_sort$sd[,idx_sc]^2 / port_sort$n[,idx_sc]) ); # decile_d plus decile_1
  port_sort$t.stat_sc = port_sort$mean_diff_sc / port_sort$sum_se_sc;
  port_sort$dof_sc = pmin(port_sort$n - 1, port_sort$n[,idx_sc] - 1);
  port_sort$pval_sc = 2*pt(-abs( port_sort$t.stat_sc ), port_sort$dof_sc);
  
  # significance, i.e. no. of stars
  port_sort$stars = stars(port_sort$pval, p.val.cutoff = c(0.05, 0.01));
  port_sort$stars_sc = stars(port_sort$pval_sc, p.val.cutoff = c(0.05, 0.01));
  
  # final
  port_sort$final = sprintf(float_print, port_sort$mean_diff) %>% paste0(port_sort$stars) %>% matrix(dim[1], dim[2]);
  port_sort$final_sc = sprintf(float_print, port_sort$mean_diff_sc) %>% paste0(port_sort$stars_sc) %>% matrix(dim[1], dim[2]);
  rownames(port_sort$final) = rownames(port_sort$final_sc) = names;
  
  # For printing in table
  temp = paste0(name__, 1);
  port_sort$print = cbind(temp = sprintf(float_print, port_sort$mean[,1]), port_sort$final[,2:n]);
  port_sort$print_sc = cbind(temp = sprintf(float_print, port_sort$mean[,1]), port_sort$final_sc[,2:n]);
  rownames(port_sort$print) = rownames(port_sort$print_sc) = names;
  
  if(modify_colnames) {
    temp = paste0(name__, 1);
    colnames(port_sort$mean_diff) = c(temp, paste0(name__, 2:n, "_min_", temp));
    colnames(port_sort$sum_se) = paste0(name__, 1:n, "_plus_", temp);
    colnames(port_sort$print) = colnames(port_sort$final) = colnames(port_sort$stars) = colnames(port_sort$pval) = colnames(port_sort$dof) = colnames(port_sort$t.stat) =
      colnames(port_sort$mean_diff);
    
    colnames(port_sort$mean_diff_sc) = c(temp, paste0(name__, 2:n, "_min_", name__, 1:(n-1)));
    colnames(port_sort$sum_se) = paste0(name__, 1:n, "_plus_", name__, idx_sc);
    colnames(port_sort$print_sc) = colnames(port_sort$final_sc) = colnames(port_sort$stars_sc) = colnames(port_sort$pval_sc) = colnames(port_sort$dof_sc) = colnames(port_sort$t.stat_sc) =
      colnames(port_sort$mean_diff_sc);
  }
  
  if(!is.null(round_digits)) {
    round_cols = c("mean", "sd", "n", "mean_diff", "sum_se", "t.stat", "dof", "pval", "mean_diff_sc", "sum_se_sc", "t.stat_sc", "dof_sc", "pval_sc");
    for(c in round_cols) {
      port_sort[[c]] = round(port_sort[[c]], round_digits);
    }
  }
  
  return(port_sort);
}












# Performs bivariate sorts on variable in cols (as of now fixed to 11 controls) and std_dev_funda. The portfolios have already been made. This function only calculates various summary statistics. It can be done both for deciles and terciles. The output is a list of statistics of 3-dimensional arrays. The 3-dimensions are 1) control variable as first sort var, 2) control var. decile and 3) std_dev_funda var decile.
# ``$print gives T1, T2-T1, T3-T1 etc while ``$print_sc gives T1, T2-T1, T3-T2 etc. Notice that the former is difference from first tercile while the latter is successive tercile changes.
bivariate_std_dev_sorts = function(dt, cols, type = "tercile", names = NULL, summ_var = "log_turn", round_digits = 3, modify_colnames = T, second_cols = NULL, factor = 1, short_name = NULL) {
  
  float_print = paste0("%0.", round_digits, "f");
  
  if(length(factor) == 1) {
    factor = rep(factor, length(cols));
  }
  
  first_cols = cols;
  
  if(is.null(short_name)) {
    short_name = substr(cols, 7, nchar(cols) - (1 + nchar(type))); # extra 1 is for underscore in _decile / _tercile
  }
  
  if(is.null(second_cols)) {
    if(type == "decile") {
      second_cols = paste0("LAG_1_10_by_10_DEP_SORT_", short_name, "_AND_std_dev_funda");
    } else {
      second_cols = paste0("LAG_1_DEP_SORT_", short_name, "_AND_std_dev_funda");
    }
  }
  
  if(type == "decile") {
    n = 10;
    idx_sc = c(1, 1:9);
    name__ = "D";
  } else {
    n = 3;
    idx_sc = c(1, 1:2);
    name__ = "T";
  }
  
  if(is.null(names)) {
    names = short_name;
  }
  
  stat_cols = c("mean_diff", "sum_se", "t.stat", "dof", "pval", "stars", "final", "print");
  stat_cols = c(c("mean", "sd", "n"), stat_cols, paste0(stat_cols, "_sc"));

  biport_sort = array(NA_real_, dim = c(length(cols), n, n, length(stat_cols)), dimnames = list(names, paste0(name__, 1:n), paste0(name__, 1:n), stat_cols));
  
  dim = dim(biport_sort);
  
  
  for(i in 1:length(cols)) {
    
    biport_sort[i,,,"mean"] = factor[i] * ( decile_view(dt, first_cols[i], second_cols[i], NULL, T, summ_var, mean, na.rm = T, print = F) %>% as.matrix );
    biport_sort[i,,,"sd"]   = factor[i] * ( decile_view(dt, first_cols[i], second_cols[i], NULL, T, summ_var, sd, na.rm = T, print = F) %>% as.matrix );
    biport_sort[i,,,"n"]    = decile_view(dt, first_cols[i], second_cols[i], NULL, T, summ_var, length, print = F) %>% as.matrix;
    
  }
  
  
  biport_sort = lapply(1:dim[4], function(i) biport_sort[,,,i]);
  names(biport_sort) = stat_cols;
  
  
  temp = paste0(name__, 1);
  
  for(i in 1:length(cols)) {
    
    # Difference in means, sum of std errors, t.stats, dof and pval
    
    biport_sort$mean_diff[names[i],,] = biport_sort$mean[names[i],,] - biport_sort$mean[names[i],,1];

    biport_sort$sum_se[names[i],,] = sqrt( (biport_sort$sd[names[i],,]^2 / biport_sort$n[names[i],,]) + (biport_sort$sd[names[i],,1]^2 / biport_sort$n[names[i],,1]) );
    
    biport_sort$t.stat[names[i],,] = biport_sort$mean_diff[names[i],,] / biport_sort$sum_se[names[i],,];
    
    biport_sort$dof[names[i],,] = pmin(biport_sort$n[names[i],,] - 1, biport_sort$n[names[i],,1] - 1);
    
    biport_sort$pval[names[i],,] = 2*pt(-abs( biport_sort$t.stat[names[i],,] ), biport_sort$dof[names[i],,]);
    
    # now do the same for successive decile changes
    
    biport_sort$mean_diff_sc[names[i],,] = biport_sort$mean[names[i],,] - biport_sort$mean[names[i],,idx_sc];
    
    biport_sort$sum_se_sc[names[i],,] = sqrt( (biport_sort$sd[names[i],,]^2 / biport_sort$n[names[i],,]) + (biport_sort$sd[names[i],,idx_sc]^2 / biport_sort$n[names[i],,idx_sc]) );
    
    biport_sort$t.stat_sc[names[i],,] = biport_sort$mean_diff_sc[names[i],,] / biport_sort$sum_se_sc[names[i],,];
    
    biport_sort$dof_sc[names[i],,] = pmin(biport_sort$n[names[i],,] - 1, biport_sort$n[names[i],,idx_sc] - 1);
    
    biport_sort$pval_sc[names[i],,] = 2*pt(-abs( biport_sort$t.stat_sc[names[i],,] ), biport_sort$dof_sc[names[i],,]);
    
    # significance, i.e. no. of stars
    
    biport_sort$stars[names[i],,] = stars(biport_sort$pval[names[i],,],
                                          p.val.cutoff = c(0.05, 0.01));
    biport_sort$stars_sc[names[i],,] = stars(biport_sort$pval_sc[names[i],,],
                                             p.val.cutoff = c(0.05, 0.01));
    
    # final
    
    biport_sort$final[names[i],,] = sprintf(float_print, biport_sort$mean_diff[names[i],,]) %>% paste0(biport_sort$stars[names[i],,]) %>% matrix(dim[2], dim[3]);
    
    biport_sort$final_sc[names[i],,] = sprintf(float_print, biport_sort$mean_diff_sc[names[i],,]) %>% paste0(biport_sort$stars_sc[names[i],,]) %>% matrix(dim[2], dim[3]);
    
    # For printing in table
    
    biport_sort$print[names[i],,] = cbind(temp = sprintf(float_print, biport_sort$mean[names[i],,1]), biport_sort$final[names[i],,2:n]);
    
    biport_sort$print_sc[names[i],,] = cbind(temp = sprintf(float_print, biport_sort$mean[names[i],,1]), biport_sort$final_sc[names[i],,2:n]);
    
  }
  

  if(modify_colnames) {
    temp = paste0(name__, 1);
    dimnames(biport_sort$mean_diff)[[3]] = c(temp, paste0(name__, 2:n, "_min_", temp));
    dimnames(biport_sort$sum_se)[[3]] = paste0(name__, 1:n, "_plus_", temp);
    dimnames(biport_sort$print)[[3]] = dimnames(biport_sort$final)[[3]] = dimnames(biport_sort$stars)[[3]] = dimnames(biport_sort$pval)[[3]] = dimnames(biport_sort$dof)[[3]] = 
      dimnames(biport_sort$t.stat)[[3]] = dimnames(biport_sort$mean_diff)[[3]];
    
    dimnames(biport_sort$mean_diff_sc)[[3]] = c(temp, paste0(name__, 2:n, "_min_", name__, 1:(n-1)));
    dimnames(biport_sort$sum_se)[[3]] = paste0(name__, 1:n, "_plus_", name__, idx_sc);
    dimnames(biport_sort$print_sc)[[3]] = dimnames(biport_sort$final_sc)[[3]] = dimnames(biport_sort$stars_sc)[[3]] = dimnames(biport_sort$pval_sc)[[3]] = 
      dimnames(biport_sort$dof_sc)[[3]] = dimnames(biport_sort$t.stat_sc)[[3]] = dimnames(biport_sort$mean_diff_sc)[[3]];
  }

  # The below works with multi-dimensional array as well
  if(!is.null(round_digits)) {
    round_cols = c("mean", "sd", "n", "mean_diff", "sum_se", "t.stat", "dof", "pval", "mean_diff_sc", "sum_se_sc", "t.stat_sc", "dof_sc", "pval_sc");
    for(c in round_cols) {
      biport_sort[[c]] = round(biport_sort[[c]], round_digits);
    }
  }
  
  return(biport_sort);
  
}








decimal_align = function(dt, sig_digits = 3) {
  
  fmt = paste0("%0.", sig_digits, "f");
  
  max_char = dt %>% sprintf(fmt, .) %>% nchar %>% max;
  
  pad = rep(" ", max_char) %>% paste0(collapse = "");
  
  dt %>% sprintf(paste0(pad, fmt), .) %>% substr(., nchar(.) - max_char + 1, nchar(.)) %>% gsub("NA", "  ", .);
  
}








# This function finds residuals from a regression specified in formula, f. The data is available in dt. residuals are assigned to new_var. If by_var is specified then regressions are run separately for each by_var. Additionally it is required that each regression is ran where valid_N non-NA data points are available
find_residuals = function(dt, f, new_var, by_var = "cusip", valid_N = 20) {
  
  if(is.null(valid_N)) {
    valid_N = 0;
  }
  
  dep_var = as.character(f)[2];
  
  ind_vars = as.character(f)[3];
  ind_vars = strsplit(ind_vars, "\\+")[[1]] %>% gsub(" ", "", .);
  # interaction terms?
  unique_ind_vars = ind_vars;
  interaction_idx = which(str_detect(ind_vars, ":"));
  if(len(interaction_idx) > 0) {
    for(i in interaction_idx) {
      unique_ind_vars = setdiff(unique_ind_vars, ind_vars[i]);
      unique_ind_vars = c(unique_ind_vars, strsplit(ind_vars[i], "\\:")[[1]] %>% gsub(" ", "", .)); 
    }
  }
  unique_ind_vars = unique(unique_ind_vars);
  
  zero_idx = which(ind_vars == "0");
  if(length(zero_idx) > 0) {
    ind_vars = ind_vars[-zero_idx];
  }
  
  reg_vars = paste0("slp__", 1:length(ind_vars));
  if(length(zero_idx) == 0) {
    reg_vars = c("int__", reg_vars);
  }
  
  expr = paste0("!is.na(", dep_var, ") & !is.na(", paste0(unique_ind_vars, collapse = ") & !is.na("), ")") %>% parse(text = .);
  
  # it may very well happen that there are no valid data to run lm. To counter that use TRY_CATCH
  expr_reg = paste0("(reg_vars) := TRY_CATCH(lm(", deparse(f), ")$coef, print.attempts = F, max.attempts = 1, ret_val = rep(NA_real_, len(reg_vars))) %>% as.list") %>% parse(text = .);
  # expr_reg = paste0("(reg_vars) := lm(", deparse(f), ")$coef %>% as.list") %>% parse(text = .);
  
  if(is.null(by_var)) {
    
    dt[eval(expr), eval(expr_reg)];
    
  } else {
    
    valid_by_vars = dt[eval(expr), .N, by = by_var][N >= valid_N][, by_var, with = F];
    
    setkeyv(dt, by_var);
    setkeyv(valid_by_vars, by_var);
    
    dt[valid_by_vars, eval(expr_reg), by = by_var];
    
  }
  
  # take care of interactions
  expr_ind_vars = str_replace(ind_vars, ":", "*");
  
  if(length(zero_idx) == 0) {
    expr2 = paste0(new_var, " := ", dep_var, " - (", paste0(reg_vars, "*", c(1, expr_ind_vars), collapse = " + "), ")") %>% parse(text = .);
  } else {
    expr2 = paste0(new_var, " := ", dep_var, " - (", paste0(reg_vars, "*", c(expr_ind_vars), collapse = " + "), ")") %>% parse(text = .);
  }
  
  dt[, eval(expr2)];
  dt[, (reg_vars) := NULL];
  
}








# the first gsub chnages existing "\\\\" to "\\". The second gsub changes all "\\" (incl. the existing and new ones) to "\\\\"
text_spec2 = function(x, format = "latex", bold = FALSE, italic = FALSE, monospace = FALSE, 
                      underline = FALSE, strikeout = FALSE, color = NULL, background = NULL, 
                      align = NULL, font_size = NULL, angle = NULL, tooltip = NULL, 
                      popover = NULL, link = NULL, extra_css = NULL, escape = F, 
                      background_as_tile = TRUE, latex_background_in_cell = FALSE) {
  
  kableExtra::text_spec(x = x, format = format, bold = bold, italic = italic, monospace = monospace, 
            underline = underline, strikeout = strikeout, color = color, background = background, 
            align = align, font_size = font_size, angle = angle, tooltip = tooltip, 
            popover = popover, link = link, extra_css = extra_css, escape = escape, 
            background_as_tile = background_as_tile, latex_background_in_cell = latex_background_in_cell) %>%
    gsub("\\\\\\\\", "\\\\", .) %>%
    gsub("\\\\", "\\\\\\\\", .);
  
}









# fill missing dates by cusip
fill_missing_dates = function(dt, date_var = "Date", firm_var = "cusip", vars = NULL) {
  
  # change the column name to something unique
  setnames(dt, firm_var, "cusip__")
  setnames(dt, date_var, "date__")
  
  setorder(dt, cusip__, date__);
  
  
  # keep only those entries which have a valid cusip and date__
  dt = dt[!is.na(cusip__) & !is.na(date__)];
  
  date.bookends = merge(dt[, date__[1], by = cusip__], # first date of each cusip__
                        dt[, date__[.N], by = cusip__], # last date of each cusip__
                        by = "cusip__");
  
  colnames(date.bookends) = c("cusip__", "first", "last");
  setorder(date.bookends, cusip__);
  
  
  unq_trading_dates = sort( dt[, unique(date__)] );
  unq_trading_dates = data.table(date__ = unq_trading_dates, start_date_num = 1:length(unq_trading_dates));
  unq_trading_dates[, end_date_num := start_date_num];
  
  date.bookends = merge(date.bookends, unq_trading_dates[, .(date__, start_date_num)], by.x = "first", by.y = "date__", all.x = T);
  date.bookends = merge(date.bookends, unq_trading_dates[, .(date__, end_date_num)], by.x = "last", by.y = "date__", all.x = T);
  
  setorder(date.bookends, cusip__);
  
  # takes 30 sec max. final is ~ 1.3 GB
  final = lapply(1:nrow(date.bookends), function(i) ( data.table(cusip__ = date.bookends$cusip__[i], date__ = unq_trading_dates$date__[ date.bookends$start_date_num[i] : date.bookends$end_date_num[i] ]) ) ) %>% rbindlist; # don't call rbindlist in a loop. It slows it down considerably (O(n^2))
  
  remove(date.bookends, unq_trading_dates);
  
  # Now merge dt on final on (cusip__, date__) pair. Many of the entries in final data.table will be NAs. Takes ~ 1 min.
  if(is.null(vars)) {
    dt = merge(final, dt, by = c("cusip__", "date__"), all.x = T);
  } else {
    vars = union(c("cusip__", "date__"), vars);
    dt = merge(final, dt[, vars, with = F], by = c("cusip__", "date__"), all.x = T);
  }
  remove(final);
  
  # change the column name to back to original
  setnames(dt, "date__", date_var)
  setnames(dt, "cusip__", firm_var)
  
  return(dt);
}











# fill missing generally
# year quarter usage: fill_missing_general(dt, first_var = "fyearq", second_var = "fqtr", FACTOR = 4)
# year month usage: fill_missing_general(dt, first_var = "year", second_var = "month", FACTOR = 12)
fill_missing_general = function(dt, first_var, second_var, FACTOR = NULL, firm_var = "cusip", vars = NULL) {
  
  dt[, new_date_tmp__ := FACTOR*get(first_var) + get(second_var) - 1];
  
  # now call fill_missing_dates()
  dt = fill_missing_dates(dt = dt, date_var = "new_date_tmp__", firm_var = firm_var, vars = vars);
  
  # convert new_date_tmp__ to first_var and second_var
  dt[, eval(first_var)  := as.integer( new_date_tmp__ / FACTOR )];
  dt[, eval(second_var) := new_date_tmp__ - FACTOR*get(first_var) + 1];
  
  # remove new_date_tmp__
  dt[, new_date_tmp__ := NULL];  
  
  return(dt);

}















# find summary of a lfe::felm type object. It's useful for storing!
lfe_summ = function(f, prcnt_explained = NULL) {
  
  n = sum(!sapply(f, is.null));
  
  if(!is.null(prcnt_explained)) {
    if(len(prcnt_explained) != n) {
      stop("length of prcnt_explained doesn't match with the number of regressions supplied in f.")
    }
  }

  if(sapply(f[1:n], class) %>% unique != "felm") {
    stop("lfe_summ can only accept objects of type: lfe::felm");
  }
  
  summ = vector("list", n);
  
  for(i in 1:n) {
    s = summary(f[[i]]);
    summ[[i]] = list(call = s$call,
                     coefficients = as.data.frame(s$coefficients),
                     r.squared = s$r.squared,
                     adj.r.squared = s$adj.r.squared,
                     P.r.squared = s$P.r.squared,
                     P.adj.r.squared = s$P.adj.r.squared,
                     N = s$N,
                     prcnt_explained = prcnt_explained[i]);
  }
  
  return(summ);
}








# Fama-MacBeth LM. Report time-series average of annual/monthly regressions
FM_lm = function(formula, data, fm_var = NULL, summarize = F) {
  if(is.null(fm_var)) {
    return(lfe::felm(formula = formula, data = data));
  } else {
    dt = data[, .N, by = fm_var][order(get(fm_var))];
    unq_fm_vars = dt[, fm_var, with = F] %>% unlist;
    data_lens = dt[, N] %>% unlist;
    dt = data.frame();
    for(u in unq_fm_vars) {
      f = lfe::felm(formula = formula, data = data[get(fm_var) == u]);
      dt = rbind(dt, c(f$coef, summary(f)$r2) %>% t);
    }
    dt = cbind(unq_fm_vars, dt, data_lens);
    rownames(dt) = 1:nrow(dt);
    colnames(dt) = c(fm_var, rownames(f$coef), "r2", "N");
    dt = as.data.table(dt);
    vars = rownames(f$coef);
    if(summarize == F) {
      return(dt);
    } else {
      means = dt[, lapply(.SD, mean, na.rm = T), .SDcols = setdiff(names(dt), fm_var)];
      se = dt[, lapply(.SD, sd, na.rm = T), .SDcols = vars] / dt[, .N];
      names(se) = paste0("std_err_", names(se));
      ret_dt = cbind(means, count = dt[, .N], se);
      for(v in vars) {
        ret_dt[, eval(paste0("tstat_", v)) := get(v) / get(paste0("std_err_", v))];
      }
      return(ret_dt);
    }
  }
}








# Guidelines for drafting a multi-panel tabel using kable:
# 1) Put the caption and label only for the first table. Remove it for rest.
# 2) After kable styling do: remove_table_tag() to remove table begin/end markers. Do this for all tables
# 3) After step 2 input panel name via: panel_text("aaa"). Do this for all tables.
# 4) For all tables except last, use use blank footnotes. Something like footnote(general_title = "", general = "", escape = F, threeparttable = T). For the last table use footnotes as usual.
# 5) Except for the first table, run fix_begin_tabular() for all other tables to fix the missing space after panel name.
# 6) If two or more panels appear on the same page then a vertical space can be added after all but last panels. Use add_vertical_space(1) to do that. You can choose any positive number of spaces.
# 7) Now enclose all the tables inside \begin{table} and \end{table}. These two lines are outside chunk.
# 8) If landscape mode is required then do not use landscape() after styling. Rather put \begin{landscape} before \begin{table} and \end{landscape} after \end{table}. These also go outside chunks. The ordering is important.
# 9) To reduce margins in a landscape table, use \newgeometry{margin=1in} before \begin{landscape} and \restoregeometry after \end{landscape}


# add a panel name in latex tables generated by kable()
# Panel name should come after \begin{threeparttable}
# Note: the below will only work of the latex output has \begin{threeparttable} in it.
panel_text = function(kable_input, t) {
  gsub("\\\\begin\\{threeparttable\\}", paste0("\\\\begin\\{threeparttable\\}", "\n", "\\\\phantom{ }\\\\\\\\", "\n", t), kable_input);
}

# remove \begin{table} and \end{table}
remove_table_tag = function(kable_input) {
  gsub("\\\\(begin|end)\\{table\\}\\[?[A-Z]?\\]?", "\n", kable_input);
}

# remove \centering
remove_centering = function(kable_input) {
  gsub("\\\\centering", "\n", kable_input);
}

# remove \begingroup and \endgroup{}
remove_grouping = function(kable_input) {
  gsub("\\\\(begingroup|endgroup\\{\\})", "", kable_input);
}

# fixes the missing gap after the second panel name in a two panel table. \begin{tabular} skips the gap but \begin{tabular}[t] adds the gap
fix_begin_tabular = function(kable_input) {
  gsub("\\\\begin\\{tabular\\}\\[t\\]", "", kable_input) %>% gsub("\\\\begin\\{tabular\\}", "\\\\begin\\{tabular\\}\\[t\\]", .);
}

add_vertical_space = function(kable_input, space = 1, extra_text = "") {
  gsub("\\\\end\\{threeparttable\\}", paste0("\\\\end\\{threeparttable\\}\n\n", "\\\\vspace\\{", space, "em\\}\n\n", extra_text), kable_input);
}











# summary stat function
# function_list = c("mean", "sd", "min", "quantile", "median", "quantile", "max", "moments::skewness", "moments::kurtosis", "length");
# function_names = c("Mean", "S.D.", "Min", "p25", "Median", "p75", "Max", "Skewness", "Kurtosis", "N");
function_list = c("mean", "sd", "min", "quantile", "median", "quantile", "max", "moments::skewness", "moments::kurtosis");
function_names = c("Mean", "S.D.", "Min", "p25", "Median", "p75", "Max", "Skewness", "Kurtosis");
function_args = rep("", length(function_list));
function_args[grep("quantile", function_list)] = c("p = 0.25", "p = 0.75");
function_na_rm = rep(T, length(function_list));
function_na_rm[grep("length", function_list)] = F;

summ_fun = function(x) {
  # list(mean(x, na.rm = T),
  #      median(x, na.rm = T),
  #      sd(x, na.rm = T));
  if(!EXISTS(function_names) || length(function_names) == 1 && function_names == "") {
    function_names = sapply(strsplit(function_list, "::"), function(i) ifelse(length(i) == 1, i, tail(i, n = 1)));
  }
  if(!EXISTS(function_args) || length(function_args) == 1 && function_args == "") {
    function_args = rep("", length(function_list));
  } else {
    function_args = ifelse(function_args == "", "", paste0(", ", function_args));
  }
  if(!EXISTS(function_na_rm) || length(function_na_rm) == 1 && function_na_rm == "") {
    function_na_rm = rep(T, length(function_list));
  }
  function_na_rm = ifelse(function_na_rm, ", na.rm = T", "");
  cmd = paste0("list(", paste0(function_names, " = ", function_list, "(x", function_na_rm, function_args, ")", collapse = ", "), ");")
  eval(parse(text = cmd));
}












# my version of column spec
my_column_spec = function(kable_input, width_vector, bold = FALSE, italic = FALSE, 
                          monospace = FALSE, underline = FALSE, strikeout = FALSE, 
                          color = NULL, background = NULL, border_left = FALSE, border_right = FALSE, 
                          width_min = NULL, width_max = NULL, extra_css = NULL, include_thead = FALSE) {
  
  N = len(width_vector);

  ret = kable_input;
  
  for(i in 1:N) {
    
    ret = column_spec(kable_input = ret, column = i, width = width_vector[i], bold = bold, italic = italic, 
                      monospace = monospace, underline = underline, strikeout = strikeout, 
                      color = color, background = background, border_left = border_left, border_right = border_right, 
                      width_min = width_min, width_max = width_max, extra_css = extra_css, include_thead = include_thead);
    
  }
  
  return(ret);
  
}









# force variable to lie within a range
force_range = function(x, x_min = 0, x_max = 1, ret_expr = F, round_digits = 3) {
  
  # new range
  A_1 = x_min;
  A_n = x_max;
  R = A_n - A_1;
  
  if(!(A_n > A_1)) {
    stop("Invalid new range provided!")
  }
  
  a_1 = min(x, na.rm = T);
  a_n = max(x, na.rm = T);
  
  if(is.infinite(a_1) | is.infinite(a_n)) {
    stop("The input has no valid entries!");
  }
  
  if(!(a_n > a_1)) {
    stop("The input vector has no variation!")
  }
  
  r = a_n - a_1; # r > 0 as well as R > 0
  
  # set range to [0,1]
  y = (x - a_1) / r;
  
  # now force range to [A_1, A_n]
  z = y*R + A_1;
  
  # z = y*R + A_1 
  #   = R*(x - a_1)/r + A_1
  #   = (R/r)*x - (R/r)*a_1 + A_1
  #   = (A_1 - (R/r)*a_1) + (R/r)*x
  if(ret_expr == T) {
    int = (A_1 - (R/r)*a_1);
    slp = (R/r);
    sig_str = paste0("%0.", round_digits, "f + %0.", round_digits, "f * x");
    return(sprintf(sig_str, int, slp));
  }
  
  return(z);
  
}








# stars in stastical significance. E.g.
# ifelse(p.val < .0001, "****", ifelse(p.val < .001, "*** ", ifelse(p.val < .01, "**  ", ifelse(p.val < .05, "*   ", "    "))));
# p.val.cutoff must be a vector such that the first entry is the cut-off between 0-star and 1-star; 2nd value is the cut-off btwn 1-star and 2-star; ... ; n-th value is the cut-off between (n-1)-star and n-star
stars = function(x, p.val.cutoff = c(0.05, 0.01, 0.001, 0.0001), symbol = "*", pad = "", tnote = F, tnote_beg_text = "\\tnote{", tnote_end_text = "}") {
  n = len(p.val.cutoff);
  
  if(n < 1) {
    stop("You must provide at-least one valid entry in p.val.cutoff");
  }
  
  p.val.cutoff = c(1, p.val.cutoff, 0);
  
  stars = rep("", len(x));
  
  for(i in 0:n) {
    if(p.val.cutoff[i+1] <= p.val.cutoff[i+2]) {
      stop("p.val.cutoff must be strictly decreasing and lie between 0 and 1 (exclusive).")
    }
    idx = which(x %between% c(p.val.cutoff[i+2], p.val.cutoff[i+1]));
    if(len(idx) > 0) {
      stars[idx] = c(rep(symbol, i), rep(pad, n-i)) %>% paste0(collapse = "");
      if(tnote == T) {
        stars[idx] = paste0(tnote_beg_text, stars[idx], tnote_end_text);
      }
    }
  }
  
  if(is.matrix(x)) {
    stars = matrix(stars, nrow = nrow(x), ncol = ncol(x));
  }
  
  return(stars);
}


















# search string in R files
search_R_files = function(search_str, match_all = T, num_lines = 2, search_dir = "C:/Users/nikhi/Dropbox/R_files_backup", search_file = NULL) {
  if(is.null(search_file)) {
    # find all R files
    all_R_files = list.files(path = search_dir, pattern = "*", recursive = T, full.names = T);
  } else {
    all_R_files = list.files(pattern = search_file);
    search_dir = ".^"; # this will never match anything
  }
  # loop through all files to look for the pattern
  f_cnt = 0
  for (f in all_R_files) {
    tmp = readLines(f, warn = F);
    search_idx = grep(search_str, tmp);
    if(len(search_idx) > 0) {
      # a valid search found in file f. Loop through all searches and print the code
      f_cnt = f_cnt + 1;
      cat(crayon::bgRed("File[#: ", f_cnt, "]: ", stringr::str_replace(f, search_dir, ""), "\n", sep = ""));
      cnt = 0;
      for (s in search_idx) {
        if(is.na(num_lines) | num_lines < 1 | is.infinite(num_lines)) {
          line_idx = 1:len(tmp);
        } else {
          line_idx = seq(s-num_lines+1, s+num_lines-1, 1);
          line_idx = line_idx[line_idx > 0 & line_idx <= len(tmp)];
        }
        cnt = cnt + 1;
        cat(crayon::yellow("Match[#: ", cnt, " @ line: ", s, "]", "\n", sep = ""));
        for(l in line_idx) {
          splits = stringr::str_split(tmp[l], search_str, simplify = T)[1,];
          final_text = c(rbind(splits, rep(crayon::bgGreen(search_str), each = len(splits))));
          final_text = final_text[1:(len(final_text)-1)];
          cat(crayon::bgCyan(l, ":", sep = ""), " ", final_text, "\n", sep = "");
        }
        cat("\n");
        if(match_all == F) {
          return();
        }
      }
      cat("\n\n", sep = "");
    }
  }
}









# The below function performs a overlapping merge between two tables: dt (main table with one date variable) AND names (secondary table with start and end dates). The required variables (oher than id and date vars) from names table should be specified in extra_names_var.
# nomatch argument controls whether non-matching names rows are returned as NA (nomatch = NA) or the corresponding dt rows are skipped altogether (nomatch = 0). Note that nomatch = 0 will change the composition of dt.
overlapping_merge = function(dt, names,
                             names_id = "cusip", names_date_vars = c("start_date", "end_date"), 
                             dt_id = "cusip", dt_date_var = "Date",
                             extra_names_var = c(), nomatch = NA) {
  
  other_dt_vars = setdiff(names(dt), c(dt_id, dt_date_var));
    
  setkeyv(names, c(names_id, names_date_vars));
  # copy date var
  dt_date_var2 = "date2__";
  dt[, eval(dt_date_var2) := get(dt_date_var)];
  # overlapping merge
  dt = foverlaps(dt, names[, c(names_id, names_date_vars, extra_names_var), with = F],
                 by.x = c(dt_id, dt_date_var, dt_date_var2), nomatch = nomatch);
  # remove id and date vars of names
  dt[, c(names_date_vars, dt_date_var2) := NULL];
  # re-order column names
  cols = c(dt_id, dt_date_var, other_dt_vars, extra_names_var);
  dt = dt[, cols, with = F];
  
  return(dt);
  
}







# memory usage of R objects
mem_use = function(n = 10, obj_list = ls(envir = parent.env(environment())), excl_fn = T) {
  ans = rep(NA, len(obj_list));
  for(i in 1:len(obj_list)) {
    if(excl_fn & is.function(get(obj_list[i]))) {
      next;
    }
    ans[i] = round(object.size(get(obj_list[i]))/2^20, 2)
  }
  names(ans) = obj_list;
  ans = ans[!is.na(ans)];
  ans = sort(ans, decreasing = T);
  if(len(obj_list) < n) {
    return(ans);
  } else {
    return(ans[1:n]);
  }
}












# Add spacing in and around \begin and \end table. Works with longtable as well.
set_custom_space = function(obj, space = 1.0) {
  
  # tabularx is not needed since we enclose that within table BUT xltabular is needed as it is not enclosed within longtable
  table_types = c("table", "longtable", "xltabular");
  
  for(t in table_types) {
    # begin part
    obj = stringr::str_replace_all(obj,
                                   paste0("\\\\begin\\{", t, "\\}"),
                                   paste0("\\\\begin\\{spacing\\}",
                                          "\\{", sprintf("%0.1f", space), "\\}",
                                          "\n",
                                          "\\\\begin\\{", t, "\\}"));
    # end part
    obj = stringr::str_replace_all(obj,
                                   paste0("\\\\end\\{", t, "\\}"),
                                   paste0("\\\\end\\{", t, "\\}",
                                          "\n",
                                          "\\\\end\\{spacing\\}"));
  }
  
  return(obj);
  
  
  # Previous Method
  
  # \\1 captures what comes before and \\3 captures what comes after. \\2 takes care of both table and longtable. (.*) matches any string of len >= 0. (.*?) matches shortest string of len >= 0.
  
  # obj %>% gsub(paste0("(.*)", "\\\\begin\\{", "(.*?)", "table\\}", "(.*)"),
  #              paste0("\\1",
  #                     "\\\\begin{spacing}{", sprintf("%0.1f", space), "}", "\n",
  #                     "\\\\begin{", "\\2", "table}",
  #                     "\\3"),
  #              .) %>%
  #         gsub(paste0("(.*)", "\\\\end\\{", "(.*?)", "table\\}", "(.*)"),
  #              paste0("\\1",
  #                     "\\\\end{", "\\2", "table}", "\n",
  #                     "\\\\end{spacing}",
  #                     "\\3"),
  #        .);
}



















