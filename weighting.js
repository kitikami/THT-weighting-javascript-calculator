

$(document).ready(function() {
    $("#parameters").submit(calculate_weight);
    $("#parameters_regr").submit(calculate_regression);
});



var v_t = ""
    , v_x = ""
    , r = ""
    , d = ""
    , mu = ""
    , size = ""
    , w = ""
    ;

// run when submitting weighting form
function calculate_weight(event) {
    event.preventDefault();
    get_values();
    var good = check_values(v_t,v_x,r,.5);  // return true if inputs are valid, false otherwise
    if(good) {
        w = optimize(reliability,lower=0,upper=1,tol=.00000001,maximum=true)[0];    // find w that maximizes correlation
        console.log(w,reliability(w));  // view unrounded figures in console.log
        $("#weight").text(Math.round(w*100000)/100000);
        $("#weight_formula").text("weight for n days ago = " + Math.round(w*100000)/100000 + "^n");
        $("#regression_submit").attr("disabled", false);    // regression form is disabled until you calculate a decay factor
    } else if(Number(r) == 0 ) {    // set weight = 0 if there is no d-t-d correlation in talent
        $("#weight").text("0");
        $("#weight_formula").text("");
    } else if(Number(r) == 1 ) {    // set weight = 1 if there are no changes in talent
        $("#weight").text("1");
        $("#weight_formula").text("");
    } else {
        $("#weight").text("");
        $("#weight_formula").text("");
    }
}

// run when submitting regression form
function calculate_regression(event) {
    event.preventDefault();
    get_values();
    var good = check_values(v_t,v_x,r,mu);
    if(good) {
        // scale variances from variance of rates to variance of successes
        var v_t_scaled = v_t * Math.pow(size,2);
        var v_x_scaled = v_x * Math.pow(size,2);
        var regr = 0, regr_i = 0;   // regr = regression constant, regr_i = regression constant without changes in talent
        [regr,regr_i] = regression_constant(v_t_scaled,v_x_scaled,mu,w,r,d);
        console.log(regr,regr_i);   // view unrounded figures in console.log
//        console.log(reliability(w));
        $("#regression_constant").text(Math.round(regr)+" ("+Math.round(regr_i)+" without changes in talent)");
    } else {
        $("#regression_constant").text("");
    }
}


// retrieve parameter values from user inputs
function get_values() {
    v_t = Number(document.getElementById('variance_true_input').value);
    v_x = Number(document.getElementById('variance_x_input').value);
    r = Number(document.getElementById('cor_input').value);
    d = document.getElementById('days_input').value;

    size = document.getElementById('size_input').value;
    mu = document.getElementById('u_input').value;
}

// verify that user inputs are valid
function check_values(v_t,v_x,cor,u) {
    var good = 0;   // counter that keeps track of how many inputs have been validated
    var hide_variance_true_error = true;
    var hide_variance_x_error = true;
    var hide_cor_error = true;
    if(Number(v_t) > 0) {   // variance must be positive
        hide_variance_true_error = true;
        good++;
    } else if(v_t=""){
        hide_variance_true_error = true;
    } else {
        hide_variance_true_error = false;
    }

    if(Number(v_x) > (v_t)) {   // overall variance must be greater than talent variance
        hide_variance_x_error = true;
        good++;
    } else if(v_x=""){
        hide_variance_x_error = true;
    } else {
        hide_variance_x_error = false;
    }

    if(Number(cor) > 0 && Number(cor) < 1) {    // d-t-d correlation in talent must be between 0 and 1 (negative values make no sense in this context)
        hide_cor_error = true;
        good++;
    } else if(cor="" || Number(cor) == 0 || Number(cor) == 1){  // handle r=0 and r=1 as special cases
        hide_variance_true_error = true;
    } else {
        hide_cor_error = false;
    }

    if(Number(u) > 0 && Number(u) < 1) {    // population mean for binomial rate must be between 0 and 1
        hide_u_error = true;
        good++;
    } else if(cor=""){
        hide_u_error = true;
    } else {
        hide_u_error = false;
    }

    $("#variance_true_error").attr("hidden", hide_variance_true_error);
    $("#variance_x_error").attr("hidden", hide_variance_x_error);
    $("#cor_error").attr("hidden", hide_cor_error);
    $("#u_error").attr("hidden", hide_u_error);


    if(good == 4){  // approve if all four inputs are valid
        return true;
    } else {
        return false;}
}

/**
 * MATH STUFF
 */

 // calculate standard deviation of weighted results
function calculate_sd_x(w,v_t,v_x,r,d) {
    sd = Math.sqrt(
        v_t * 2 *
        (
            (r * w - Math.pow((r * w),d)) / ((1 - r * w) * (1 - Math.pow(w,2))) +
            (Math.pow(w,(2 * d)) * (Math.pow((r / w),d) - r / w)) / ((1 - r / w) * (1 - Math.pow(w,2)))
        ) +
        v_x * (1 - Math.pow(w,(2 * d))) / (1 - Math.pow(w,2))
    );
    return sd;
}

// calculate expected correlation between weighted results and talent
function reliability(w) {
    get_values();
    sd_x = calculate_sd_x(w,v_t,v_x,r,d);
    sd_y = Math.sqrt(v_t);
    cov = v_t * (1-Math.pow((r*w),d))/(1-r*w);
//    return -1 * cov/(sd_x*sd_y);  // use if optimize function can only find minimum and not maximum
    return cov/(sd_x*sd_y);
}


// calculate regression constant
function regression_constant(v_t,v_x,mu,w,r,d) {
    // calculate regression constant based on reliability of weighted sample to itself
    // projects weighted average of player's talent over sample
    var sd_t = Math.sqrt(v_t)/size;                             // standard deviation in current talent
    var sd_x = calculate_sd_x(w,v_t,v_x,r,d);                   // standard deviation in weighted results
    var trials = size * (1-Math.pow(w,d))/(1-w);                // effective sample size after weighting
    var bias = trials/(trials-1);                               // adjustment necessary to make var_true + var_rand = var_observed
    var var_rnd = size * mu*(1-mu)*(1-Math.pow(w,(2*d)))/(1-Math.pow(w,2));     // random binomial variance for weighted smample
    var var_true = bias * (Math.pow(sd_x,2) - var_rnd);         // variance in true talent of full results
    var var_explained = bias * (1 - var_rnd/Math.pow(sd_x,2));  // proportion of variance explained by spread in talent
    var regr = trials * (1-var_explained) / var_explained;      // preliminary regression constant (will have to be adjusted)
    var regr_i =  mu*(1-mu)/Math.pow(sd_t,2);                   // regression constant ignoring changes in talent

    // test values to troubleshoot if something goes wrong
    // console.log(v_t,v_x,bias,var_rnd,var_true,var_explained,regr,regr_i);

    // adjust regression constant to project current talent rather than talent over weighted sample
    var A = (1 - w)*(1 - Math.pow((r*w),d))/((1 - r*w)* (1 - Math.pow(w,d)));
    var B = 2 *
			(
				(r*w - Math.pow((r*w),d)) / ( (1-r*w) * (1-Math.pow(w,2)) ) +
				(Math.pow(w,(2*d))*(Math.pow((r/w),d) - r/w)) / ( (1-r/w) * (1-Math.pow(w,2)) )
			) +
		(1 - Math.pow(w,(2*d)))/(1-Math.pow(w,2));
    var C = 2 *
			(
				(w - Math.pow(w,d)) / ( (1-w) * (1-Math.pow(w,2)) ) +
				(Math.pow(w,(2*d))*(Math.pow((1/w),d) - 1/w)) / ( (1-1/w) * (1-Math.pow(w,2)) )
			) +
		(1 - Math.pow(w,(2*d)))/(1-Math.pow(w,2));

    // regression for talent from weighted sample to single day
    var sample_size_correction = A/(B/C);
    var_explained = var_true * Math.pow(sample_size_correction,4)/(var_true * Math.pow(sample_size_correction,4) + var_rnd)
    regr = trials * (1-var_explained) / var_explained

    // test values to troubleshoot if something goes wrong
    // console.log(A,B,C,sample_size_correction);

    return [regr,regr_i];
}



