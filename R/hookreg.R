#' A function to calculate the slope and intercept of an amplification curve data
#' from a quantitative PCR experiment at the end of the data stream.
#'
#' \code{hookreg} is a function to calculate the slope and intercept of an
#' amplification curve data from a quantitative PCR experiment. The idea is that
#' a strong negative slope at the end of an amplification curve is indicative for
#' a hook effect (see Barratt and Mackay 2002).
#' @return gives a \code{numeric} (S3 class, type of \code{double}) as output 
#' for the detection of a hook
#'
#' @param x is the cycle numbers (x-axis).
#' @param y is the cycle dependent fluorescence amplitude (y-axis).
#' @param normalize is a logical parameter indicating if the data should be
#' normalized to the 0.999 quantile
#' @param sig.level defines the significance level to test for a significant
#' regression
#' @param CI.level confidence level required for the slope
#' @param robust is a logical parameter indicating if the data should be
#' analyzed be a robust linear regression (\code{lmrob}).
#' @author Stefan Roediger, Michal Burdukiewcz
#' @references K. Barratt, J.F. Mackay, \emph{Improving Real-Time PCR Genotyping
#' Assays by Asymmetric Amplification}, J. Clin. Microbiol. 40 (2002) 1571--1572.
#' doi:10.1128/JCM.40.4.1571-1572.2002.
#' @keywords slope intercept hook
#' @examples
#' default.par <- par(no.readonly = TRUE)
#' # Calculate slope and intercept on noise (negative) amplification curve data
#' # for the last eight cycles.
#'
#' library(qpcR)
#'
#' res_hook <- data.frame(sample=colnames(boggy)[-1], 
#'                        t(sapply(2:ncol(boggy), function(i) {
#'                        hookreg(x=boggy[, 1], y=boggy[, i])})))
#' res_hook
#'
#' data_colors <- rainbow(ncol(boggy[, -1]), alpha=0.5)
#' cl <- kmeans(na.omit(res_hook[, 2:3]), 2)$cluster
#'
#' par(mfrow=c(1,2))
#' matplot(x=boggy[, 1], y=boggy[, -1], xlab="Cycle", ylab="RFU",
#'  main="boggy Data Set", type="l", lty=1, lwd=2, col=data_colors)
#'  legend("topleft", as.character(res_hook$sample), pch=19,
#'          col=data_colors, bty="n")
#'
#' plot(res_hook$intercept, res_hook$slope, pch=19, cex=2, col=data_colors,
#'  xlab="intercept", ylab="Slope",
#'  main="Clusters of Amplification Curves with an Hook Effect-like Curvature\nboggy Data Set")
#'  points(res_hook$intercept, res_hook$slope, col=cl, pch=cl, cex=cl)
#'  legend("topright", c("Strong Hook effect", " Weak Hook effect"), pch=c(1,2), col=c(1,2), bty="n")
#'  text(res_hook$intercept, res_hook$slope, res_hook$sample)
#'
#' par(default.par)
#' @export hookreg

hookreg_strict <- function(x, y, normalize = TRUE, sig.level = 0.0005, CI.level = 0.999, 
                           min_hook_delta = 8, robust = FALSE) {
    # Remove NA values
    data <- na.omit(data.frame(x = x, y = y))
    x <- data[, "x"]
    y <- data[, "y"]
    
    # Normalize the y values if required
    if (normalize) {
        y <- y / quantile(y, 0.999)
    }
    
    # Define the range starting from the point where y is in the top 25% of its values
    hook_quantile_range <- which(y >= quantile(y, c(0.75)))[1]:length(y)
    hook_max_range <- which(y == max(y[hook_quantile_range]))[1]
    range <- hook_max_range:length(x)
    hook_delta <- length(range)
    
    # Check if the range is long enough for meaningful analysis
    if (hook_max_range < length(x) && hook_delta >= min_hook_delta) {
        # Perform robust or regular linear regression on the selected range
        if (robust) {
            res_lm_fit <- try(lmrob(y[range] ~ x[range]), silent = TRUE)
        } else {
            res_lm_fit <- try(lm(y[range] ~ x[range]), silent = TRUE)
        }
        
        # Extract regression results if the fit was successful
        if (inherits(res_lm_fit, "try-error")) {
            return(rep(NA, 10))
        }
        
        res_lm_fit_summary <- summary(res_lm_fit)$coefficients[2, 4]  # p-value for slope
        res_lm_fit_coefficients <- coefficients(res_lm_fit)
        res_lm_fit_confint <- confint(res_lm_fit, level = CI.level)  # confidence intervals
        
        # Apply stricter check: Ensure both confidence intervals are below 0 and even below -0.85
        res_lm_fit_confint_decision <- ifelse(
          res_lm_fit_confint[2, 1] < -0.85 && res_lm_fit_confint[2, 2] < -0.85, 
          TRUE, FALSE
        )
        
        # Significance test on p-value
        res_hook_significance <- ifelse(res_lm_fit_summary < sig.level, TRUE, FALSE)
        
        # Final decision based on stricter criteria
        dec_hook <- ifelse(res_hook_significance == TRUE || res_lm_fit_confint_decision == TRUE, TRUE, FALSE)
        
        # Return the results
        res_hookreg <- c(
            res_lm_fit_coefficients[[1]], res_lm_fit_coefficients[[2]], 
            hook_max_range, hook_delta, res_lm_fit_summary, 
            res_lm_fit_confint[2, 1], res_lm_fit_confint[2, 2], 
            res_hook_significance, res_lm_fit_confint_decision, dec_hook
        )
    } else {
        # Return default values if no hook is detected
        res_hookreg <- c(0, 0, 0, 0, NA, NA, NA, FALSE, FALSE, FALSE)
    }
    
    # Name the result fields
    names(res_hookreg) <- c("intercept", "slope", "hook.start", "hook.delta", "p.value", 
                            "CI.low", "CI.up", "hook.fit", "hook.CI", "hook")
    return(res_hookreg)
}
