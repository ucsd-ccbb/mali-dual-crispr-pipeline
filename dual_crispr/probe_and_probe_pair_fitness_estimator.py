# third-party libraries
import numpy
import pandas

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


# FYI: Throughout, "tuple" is used as shorthand to refer to a tuple of (gene_name, probe_name) that uniquely
# identifies a gene/probe pair.  Thus, a data object that is "per_tuple" is per probe (given the assumption that
# each probe can only be for one gene), while a data object that is "by_tuple_by_tuple" is matrix-like and symmetric,
# with each axis indexed with entries for all gene/probe pairs in the experiment, and cell values being values of
# whatever the measured data is for the construct made up of the intersecting gene/probe-gene/probe pair.  Note that not
# all possible gene/probe-gene/probe pairs exist in actual constructs, so some cells will have NaN or default (usually
# zero) values.
#
# Alternately, when I designate a data object "per_probe" or "by_probe_by_probe", I mean that it is NOT indexed--these
# objects are generally either one-d or two-d numpy arrays, respectively. The dimensions of their axes represent the
# exact same thing as the "per_tuple" or "by_tuple_by_tuple" data objects (and in the same order), but since they are
# numpy arrays the slots along those dimensions are not explicitly labeled with what gene/probe pair they pertain to.
# Thus, these should be used with care.

# Function generalizes functionality from Roman's MethodII.R lines 161-255, shown below with additional comments by ab:
#
# ab: fc is a symmetric probe-by-probe matrix where the value in each cell is the fc value for the construct including
# ab: those two probes, or zero if there is no construct for that probe-by-probe combination.
# ab: w0 is a symmetric probe-by-probe matrix of initial weights, where the value in each cell is zero for "bad"
# ab: and/or non-existent constructs and one for "good" constructs.
# ab: The manuscript gives eqn 2 as:
# ab: fc = fp + fp' + pi sub pp'
# ab: and states "The probe fitnesses are found by robust fitting of Eq. (2).
# ab: The probe-level pi-scores are the residuals of the robust fit."
# #iterative robust least squares
# irls<-function(fc,w0,probes,ag=2,tol=1e-3,maxiter=30) {
# # w0 is the physical goodness of constructs. It is not subject to change.
# # It is used to modify w, to silence bad constructs
#
#    ab:upper.tri seems to be getting the upper triangle of the symmetric probe-by-probe fc matrix
#    expressed_utri<-upper.tri(fc) & w0>0
#    ab: lessee--that would be number of probes; note that n here is *different* than n outside the scope of this
#    ab: function, where it is the number of genes :(
#    n<-dim(fc)[1]
#    ab: probe-by-probe matrix with default values of 1 everywhere except on diagonal, where they are zero
#    w<-matrix(1,nrow=n,ncol=n) #initial weights
#    diag(w)<-0
#    fij<-matrix(0,nrow=n,ncol=n) #initial weights
#    eij<-matrix(0,nrow=n,ncol=n) #initial weights
#    b<-rep(0,n) #rhs
#
#    #iteration step 0
#    w<-w*w0  # I think here we are "silenc[ing] bad constructs", which have w0 = 0
#    A<-w
#    ab: for each probe
#    for (i in 1:n) {
#       ab: sum of fc*weight for all probes paired with this probe.  Why set i as the second index into in the
#       ab: probe-by-probe matrix rather than the first?  Isn't the matrix symmetric?
#       b[i]<-sum(fc[,i]*w[,i])
#       ab: reminder: small<-1e-6 ; it is a parameter hardcoded at the top of the whole notebook.
#       ab: so, A = w and w is 1 everywhere except on the diagonal and at probe-by-probe pairs where the analogous
#       ab: construct is either bad or non-existent.  This seems to be setting the values along the diagonal to the
#       ab: sum of the weights for that probe (across all other probes it is paired with) plus a tiny value, presumably
#       ab: to make sure the diagonal is now *never* zero.
#       A[i,i]<-sum(w[,i])
#    }
#
#    ab: apparently, "solve() function solves equation a %*% x = b for x, where b is a vector or matrix."
#    ab: (from http://www.endmemo.com/program/R/solve.php). Note that %*% is apparently R notation for
#    ab: "multiply matrices".  Ok, so figure out what matrix you need to multiply A by to get b, and call that matrix y
#    ab: (again, note y here is *different* from y outside this function, where it is #log2 frequencies :( ) Here, y is
#    ab: a matrix with 1 column and as many rows as there are probes, since b is an array (essentially, a matrix with
#    ab: 1 row and as many columns as there are probes)
#    y<-solve(A,b)
#    names(y)<-probes
#    ab: ok, here's that same logic for getting all pairs of probes
#    for (i in 1:(n-1)) {
#       for (j in (i+1):n) {
#          fij[i,j]<-y[i]+y[j]
#       }
#    }
#    fij<-fij+t(fij)  # ab: make fij symmetric
#    eij<-fc-fij #residuals
#    ab: end iteration step 0
#
#    l<-1 #counter
#    rel<-1  # ab: apparently, to judge by Roman's comment at the end of this while loop, "rel" = "relative error"
#    while (rel > tol & l < maxiter) {#iterate until tol is reached or maxit
#       ab: iterate until tol is reached or maxit
#       ab: calc stddev of residuals for all expressed probes
#       s<-sd(eij[expressed_utri]) #calculate sd only from expressed constructs
#       yold<-y  # ab: archive y, the matrix w/ one column and rows for each probe
#
#       w<-(1-eij^2/(ag*s)^2)^2 #something like Tukey's biweight
#       w[abs(eij)>(ag*s)]<-0
#       diag(w)<-0
#
#       w<-w*w0
#
#       A<-w
#       for (i in 1:n) {
#          b[i]<-sum(fc[,i]*w[,i])
#          A[i,i]<-sum(w[,i])
#       }
#       y<-solve(A,b)
#       names(y)<-probes
#       fij<-matrix(0,nrow=n,ncol=n) #initial weights
#       for (i in 1:(n-1)) {
#          for (j in (i+1):n) {
#             fij[i,j]<-y[i]+y[j]
#          }
#       }
#       fij<-fij+t(fij)
#       eij<-fc-fij #residuals
#
#       ab: ok, in iteration 0, rel = one and l (the counter) = one
#       ab: here, the counter = the counter plus one
#       ab: rel: we're calculating the difference of the y value for this time through the loop
#       ab: from the y value the last time through the loop and then squaring it (presumably to get rid of
#       ab: any negatives); we do this across all probes and sum the squared values.  Then we divide that sum
#       ab: of squares by either one or the sum across all probes of the square of the y values from the last time
#       ab: through the loop, whichever is greater.  Then we take the sqrt of that ratio and say it is the
#       ab: "relative error".
#
#       rel<-sqrt(sum((y-yold)^2)/max(1,sum(yold^2))) #relative error
#       ab: reminder: sqrtsum<-function(y) sqrt(sum(y^2))
#       cat(l,sqrtsum(yold),sqrtsum(y-yold),"\n")  # ab: print out some values to the screen
#       l<-l+1
#
#       ab: we break out of this loop when either:
#       ab: a) the relative error has fallen below the "tol" value (what does "tol" stand for?) OR
#       ab: b) we have performed the maximum number of trips through the loop.
#       ab: Both tol and maxiter are parameters passed into this function.
#     }
#    ab: Roman's comment when y is unpacked from this returned list is "#these are probe fitnesses fp"
#    ab: Roman's comment when eij is unpacked from this returned list is "#raw pi-scores per construct"
#    vl<-list(y, fij, eij)
#    return(vl)
# }
#
# TODO: Find out from Roman what in the heck the "ag" variable *is* so I can give it a better name
def _calc_unnormed_probe_fitnesses_and_probe_pair_pi_scores(unnormed_construct_fitnesses_by_tuple_by_tuple_df,
                                                            initial_weights_by_tuple_by_tuple_df, ag=2,
                                                            max_relative_error=1e-3, stop_at_num_iterations=30):
    # i.e., inputs are fc, w0, tol, and maxiter

    # initialize the per-probe fitness fits and the per-probe residuals (pi scores)
    unnormed_fitness_per_tuple_series = _calc_unnormed_fitness_per_tuple_series(
        unnormed_construct_fitnesses_by_tuple_by_tuple_df, initial_weights_by_tuple_by_tuple_df, None)

    pi_scores_by_tuple_by_tuple_df = _calc_pi_scores_by_tuple_by_tuple_df(
        unnormed_construct_fitnesses_by_tuple_by_tuple_df, unnormed_fitness_per_tuple_series)
    num_iterations = 1
    relative_error = 1  # initialize to large number
    # i.e., expressed_utri
    by_probe_by_probe_upper_triangle_w_nonzero_weights_mask = numpy.logical_and(
        numpy.triu(unnormed_construct_fitnesses_by_tuple_by_tuple_df), initial_weights_by_tuple_by_tuple_df > 0)

    # Do iterative robust least squares fitting:
    while relative_error > max_relative_error and num_iterations < stop_at_num_iterations:
        # stddev_of_residuals_for_expressed_probe_pairs (i.e., s) and ag_times_stddev are scalars.
        # i.e., w <- (1 - eij ^ 2 / (ag * s) ^ 2) ^ 2
        stddev_of_residuals_for_expressed_probe_pairs = numpy.std(
            pi_scores_by_tuple_by_tuple_df[by_probe_by_probe_upper_triangle_w_nonzero_weights_mask])
        ag_times_stddev = ag * stddev_of_residuals_for_expressed_probe_pairs
        current_weights_by_probe_by_probe_2darray = numpy.square(
            (1 - numpy.square(pi_scores_by_tuple_by_tuple_df.values)) / numpy.square(ag_times_stddev))

        # i.e., w[abs(eij) > (ag * s)] <- 0
        large_residuals_mask_by_probe_by_probe_array = numpy.absolute(
            pi_scores_by_tuple_by_tuple_df.values) > ag_times_stddev
        current_weights_by_probe_by_probe_2darray[large_residuals_mask_by_probe_by_probe_array] = 0
        # i.e., diag(w) <- 0
        diag_indices_array = numpy.diag_indices_from(current_weights_by_probe_by_probe_2darray)
        current_weights_by_probe_by_probe_2darray[diag_indices_array] = 0

        # recalculate the per-probe fitness fits and the per-probe residuals (pi scores)
        previous_unnormed_fitness_per_tuple_series = unnormed_fitness_per_tuple_series
        unnormed_fitness_per_tuple_series = _calc_unnormed_fitness_per_tuple_series(
            unnormed_construct_fitnesses_by_tuple_by_tuple_df, initial_weights_by_tuple_by_tuple_df,
            current_weights_by_probe_by_probe_2darray)
        pi_scores_by_tuple_by_tuple_df = _calc_pi_scores_by_tuple_by_tuple_df(
            unnormed_construct_fitnesses_by_tuple_by_tuple_df, unnormed_fitness_per_tuple_series)

        # i.e. rel <- sqrt(sum((y - yold) ^ 2) / max(1, sum(yold ^ 2)))
        sum_of_diffs_in_new_vs_old_unnormed_fitness_per_probe_array = numpy.sum(
            unnormed_fitness_per_tuple_series.values - previous_unnormed_fitness_per_tuple_series.values, axis=0)
        sum_squared_previous_unnormed_fitnesses = numpy.sum(numpy.square(
            previous_unnormed_fitness_per_tuple_series.values))  # sum_squared_previous_unnormed_fitnesses is a scalar
        denominator = max([1, sum_squared_previous_unnormed_fitnesses])  # denominator is a scalar
        relative_error = numpy.sqrt(numpy.square(sum_of_diffs_in_new_vs_old_unnormed_fitness_per_probe_array) /
                                    denominator)  # relative error is a scalar

        # TODO: add logging of cat statement info?
        # cat(l, sqrtsum(yold), sqrtsum(y - yold), "\n")
        num_iterations += 1
    # end while loop

    # Roman's original code also sends back fitted_fc_by_probe_by_probe_2darray (fij), but it is never actually used
    return unnormed_fitness_per_tuple_series, pi_scores_by_tuple_by_tuple_df  # i.e., y, eij


# Function generalizes functionality from Roman's MethodII.R lines 168-169 and 175-182, shown below with additional
# comments by ab:
#
# ab: w0 is a symmetric probe-by-probe matrix of initial weights, where the value in each cell is zero for "bad"
# ab: and/or non-existent constructs and one for "good" constructs.
#
#    ab: probe-by-probe matrix with default values of 1 everywhere except on diagonal, where they are zero
#    w<-matrix(1,nrow=n,ncol=n) #initial weights
#    diag(w)<-0
#
#    w<-w*w0  # I think here we are "silenc[ing] bad constructs", which have w0 = 0
#    A<-w
#    ab: for each probe
#    for (i in 1:n) {
#       ab: sum of fc*weight for all probes paired with this probe.  Why set i as the second index into in the
#       ab: probe-by-probe matrix rather than the first?  Isn't the matrix symmetric?
#       b[i]<-sum(fc[,i]*w[,i])
#       ab: reminder: small<-1e-6 ; it is a parameter hardcoded at the top of the whole notebook.
#       ab: so, A = w and w is 1 everywhere except on the diagonal and at probe-by-probe pairs where the analogous
#       ab: construct is either bad or non-existent.  This seems to be setting the values along the diagonal to the
#       ab: sum of the weights for that probe (across all other probes it is paired with) plus a tiny value, presumably
#       ab: to make sure the diagonal is now *never* zero.
#       A[i,i]<-sum(w[,i])
#    }
#
#    ab: apparently, "solve() function solves equation a %*% x = b for x, where b is a vector or matrix."
#    ab: (from http://www.endmemo.com/program/R/solve.php). Note that %*% is apparently R notation for
#    ab: "multiply matrices".  Ok, so figure out what matrix you need to multiply A by to get b, and call that matrix y
#    ab: (again, note y here is *different* from y outside this function, where it is #log2 frequencies :( ) Here, y is
#    ab: a matrix with 1 column and as many rows as there are probes, since b is an array (essentially, a matrix with
#    ab: 1 row and as many columns as there are probes)
#    y<-solve(A,b)
#    names(y)<-probes
#
# Note that the R code makes its calculations one probe at a time, while this calculates for all at once.
#
# TODO: unclear if this functionality is even relevant to single-probe construct case, because in that case,
# the unnormed_construct_fitnesses ARE the unnormed_probe_fitnesses, so seems like all this reweighting and refitting
# to minimize residuals may be unnecessary?? Check w/Roman.
def _calc_unnormed_fitness_per_tuple_series(unnormed_construct_fitnesses_by_tuple_by_tuple_df,  # i.e., fc
                                            initial_weights_by_tuple_by_tuple_df,               # i.e., w0
                                            current_weights_by_probe_by_probe_2darray):         # i.e., w

    num_probes = _get_num_probes(unnormed_construct_fitnesses_by_tuple_by_tuple_df.values)

    # TODO: Think this would have to be modified to make code work for single-probe constructs--current_weights
    # would be just a 1d array by probe, populated with ones
    if current_weights_by_probe_by_probe_2darray is None:
        current_weights_by_probe_by_probe_2darray = numpy.ones((num_probes, num_probes))
        current_weights_by_probe_by_probe_2darray[numpy.diag_indices(num_probes)] = 0

    # Think this code would work with single-probe constructs, assuming change above to current weights was made and
    # that initial_weights_by_tuple_by_tuple_df was really initial_weights_by_tuple_df (a series represented as a df)
    working_weights_by_probe_by_probe_2darray = current_weights_by_probe_by_probe_2darray * \
                                                initial_weights_by_tuple_by_tuple_df.values

    # I do not think the below call needs to be made in the single-probe construct case; since there is only one
    # probe in each row of the working weights "matrix", the summing of probe weights that is done in this call is
    # irrelevant
    revised_weights_matrix_2darray = _reset_diagonal_of_weights_by_probe_by_probe(
        working_weights_by_probe_by_probe_2darray)

    # Think this code would work with single-probe constructs, assuming change above to weights was made and
    # that unnormed_construct_fitnesses_by_tuple_by_tuple_df was really unnormed_construct_fitnesses_by_tuple_df
    # (a series represented as a df)
    sum_fc_times_weight_across_probe_per_probe_array = numpy.sum(
        unnormed_construct_fitnesses_by_tuple_by_tuple_df.values * working_weights_by_probe_by_probe_2darray, axis=1)

    # TODO: Think this would have to be modified to make code work for single-probe constructs--R solve function is
    # based on matrix assumptions, so probably need a fundamentally different approach. Maybe just a simple linear
    # regression is called for?
    unnormed_fitness_per_probe_array = _calc_per_probe_fits_in_r(revised_weights_matrix_2darray,
                                                                 sum_fc_times_weight_across_probe_per_probe_array)

    # Rest of this should work fine with single-probe constructs
    unnormed_fitness_per_tuple_series = pandas.Series(unnormed_fitness_per_probe_array,
                                                      index=unnormed_construct_fitnesses_by_tuple_by_tuple_df.index)
    return unnormed_fitness_per_tuple_series  # i.e., y


# Function generalizes functionality from Roman's MethodII.R lines 175-180, shown below with additional comments by ab:
#
# ab: w = current weights times initial weights
# ab: n = number of probes, i = index of current probe
# ab: fc is a symmetric probe-by-probe matrix where the value in each cell is the fc value for the construct including
# ab: those two probes, or zero if there is no construct for that probe-by-probe combination.
#
#    A<-w
#    ab: for each probe
#    for (i in 1:n) {
#       ...
#       ab: reminder: small<-1e-6 ; it is a parameter hardcoded at the top of the whole notebook.
#       ab: so, A = w and w is 1 everywhere except on the diagonal and at probe-by-probe pairs where the analogous
#       ab: construct is either bad or non-existent.  This seems to be setting the values along the diagonal to the
#       ab: sum of the weights for that probe (across all other probes it is paired with) plus a tiny value, presumably
#       ab: to make sure the diagonal is now *never* zero.
#       A[i,i]<-sum(w[,i])
#    }
#
# TODO: why does scoring_main.R say A[i, i] <- sum(w[, i]) + small but MethodII.R says just A[i,i]<-sum(w[,i]) ?
# Ask Roman which it should be!  This code implements the version that adds small.
#
# Note that the R code makes its calculations one probe at a time, while this calculates for all at once.
#
# I do not think the below code needs to be called in the single-probe construct case; since there is only one
# probe in each row of the working weights "matrix", the summing of probe weights that is done in this function is
# irrelevant (and the setting of "diagonal" values would be nonsensical and probably error out).
def _reset_diagonal_of_weights_by_probe_by_probe(working_weights_by_probe_by_probe_2darray,  # i.e., w
                                                 small_positive_value=1e-6):                 # i.e., small
    # TODO: is goal here to seed revised_weights_matrix from working_weights_matrix but not have latter change as we
    # make modifications to the former?  If so, do I need to do something more complex to ensure this in Python? copy?
    revised_weights_matrix_2darray = working_weights_by_probe_by_probe_2darray
    # note that this nonzero sum is NOT what is used to get sum_fc_times_weight_across_probe_per_probe ...

    sum_weight_across_probe_made_nonzero_per_probe = numpy.sum(working_weights_by_probe_by_probe_2darray, axis=1) + \
        small_positive_value
    num_probes = _get_num_probes(working_weights_by_probe_by_probe_2darray)
    revised_weights_matrix_2darray[numpy.diag_indices(num_probes)] = sum_weight_across_probe_made_nonzero_per_probe
    return revised_weights_matrix_2darray   # i.e., A


# I think the code in this function will work with single-probe constructs if the by_probe_by_probe array
# are actually just by_probe--so a 1d array.
def _get_num_probes(any_by_probe_by_probe_2darray):
    # TODO: Add code to ensure both dimensions are the same if there is >1 dimension
    return any_by_probe_by_probe_2darray.shape[0]


# Function generalizes functionality from Roman's MethodII.R lines 181-182, shown below with additional comments by ab:
#
#    ab: A = current weights * initial weights
#    ab: b = sum of construct fitness*weight for all probes paired with this probe.
#    ab: apparently, "solve() function solves equation a %*% x = b for x, where b is a vector or matrix."
#    ab: (from http://www.endmemo.com/program/R/solve.php). Note that %*% is apparently R notation for
#    ab: "multiply matrices".  Ok, so figure out what matrix you need to multiply A by to get b, and call that matrix y
#    ab: (again, note y here is *different* from y outside this function, where it is #log2 frequencies :( ) Here, y is
#    ab: a matrix with 1 column and as many rows as there are probes, since b is an array (essentially, a matrix with
#    ab: 1 row and as many columns as there are probes)
#    y<-solve(A,b)
#    names(y)<-probes
def _calc_per_probe_fits_in_r(revised_weights_matrix_2darray, sum_fc_times_weight_across_probe_per_probe_array):

    # TODO: R code seems to expect sum_fc_times_weight_across_probe_per_probe_array to be  1 row long and num_probes
    # columns wide, but is currently num_probes rows long and 1 column wide ... need to transpose before calling solve?

    # TODO: fill in code; determine if need to do in R natively or can do in rpy2
    # y <- solve(A, b)
    # A = revised_weights_matrix_2darray
    # b = sum_fc_times_weight_across_probe_per_probe_array
    raise NotImplementedError
    return unnormed_fitness_per_probe_array  # i.e., y


# Function generalizes functionality from Roman's MethodII.R lines 183-189, shown below with additional comments by ab:
#
# ab: n = number of probes
# ab: y = a matrix with 1 column and as many rows as there are probes, holding the fitted fitness for each probe
# ab: fc is a symmetric probe-by-probe matrix where the value in each cell is the fc value for the construct including
# ab: those two probes, or zero if there is no construct for that probe-by-probe combination.
#    for (i in 1:(n-1)) {
#       for (j in (i+1):n) {
#          fij[i,j]<-y[i]+y[j]
#       }
#    }
#    fij<-fij+t(fij)  # ab: make fij symmetric
#    eij<-fc-fij #residuals
#
# I think this code (with the exception of the call to _extract_fitted_fc_by_probe_by_probe_2darray) will work with
# single-probe constructs if all the by_tuple_by_tuple data structures are actually just by_tuple, with only one
# column (so, in fact, equivalent to series, but represented as dataframes).
def _calc_pi_scores_by_tuple_by_tuple_df(unnormed_construct_fitnesses_by_tuple_by_tuple_df,  # i.e., fc
                                         unnormed_fitness_per_tuple_series):                 # i.e., y

    """

    Returns:
        pandas.DataFrame:
    """
    fitted_fc_by_probe_by_probe_2darray = _extract_fitted_fc_by_probe_by_probe_2darray(
        unnormed_fitness_per_tuple_series)
    current_residuals_by_probe_by_probe_2darray = unnormed_construct_fitnesses_by_tuple_by_tuple_df.values - \
                                                  fitted_fc_by_probe_by_probe_2darray
    pi_scores_by_tuple_by_tuple_df = pandas.DataFrame(
        current_residuals_by_probe_by_probe_2darray,
        columns=unnormed_construct_fitnesses_by_tuple_by_tuple_df.columns.values,
        index=unnormed_construct_fitnesses_by_tuple_by_tuple_df.index)
    return pi_scores_by_tuple_by_tuple_df  # i.e., eij


# Function generalizes functionality from Roman's MethodII.R lines 183-186, shown below with additional comments by ab:
#
# ab: n = number of probes
# ab: y = a matrix with 1 column and as many rows as there are probes, holding the fitted fitness for each probe
# ab: fij = fitted_fitness_by_probe_by_probe_2darray
#    for (i in 1:(n-1)) {
#       for (j in (i+1):n) {
#          fij[i,j]<-y[i]+y[j]
#       }
#    }
#    fij<-fij+t(fij)  # ab: make fij symmetric
#
# TODO: This function WILL NOT WORK if data is single-probe, since it explicitly expects two probes per construct--
# will need to be refactored to support single-probe case
# TODO: Should I integrate 2d case with scoring_data.make_by_tuple_by_tuple_annotated_symmetric_df?  Seems repetitive
def _extract_fitted_fc_by_probe_by_probe_2darray(unnormed_fitness_per_tuple_series):  # i.e., y
    # set fitted_fc for constructs with *non-identical* probe pairs to the sum of the fits for the two different
    # probes in the construct.  Do nothing for constructs with identical probes--leave values as zero (default).
    # Note that this only fills one half (one triangle) of a square matrix; is made symmetric afterwards.
    unnormed_fitness_per_probe_array = unnormed_fitness_per_tuple_series.values
    num_probes = unnormed_fitness_per_probe_array.shape[0]
    fitted_fc_by_probe_by_probe_2darray = numpy.zeros((num_probes, num_probes))
    for first_probe_index in range(0, num_probes-1):
        for second_probe_index in range(1, num_probes):
            fitted_fc_by_probe_by_probe_2darray[first_probe_index, second_probe_index] = \
                unnormed_fitness_per_probe_array[first_probe_index] + \
                unnormed_fitness_per_probe_array[second_probe_index]

    # make this 2d array symmetric by adding it to its own transpose
    fitted_fc_by_probe_by_probe_2darray = fitted_fc_by_probe_by_probe_2darray + \
                                          fitted_fc_by_probe_by_probe_2darray.transpose()
    return fitted_fc_by_probe_by_probe_2darray  # i.e., fij
