
# multilevelMatching 0.1.0.9002

- Added `multiMatch()` function to carry out all types of matching. This effectively combines `multilevelGPSMatch()` and `multilevelMatchX()` into one function. Features include:
  - Better output: tidier estimates - see `estimateTau()` - plus more information from under-the-hood - see the `impute_mat` object for matrix off all imputed potential outcomes.
  - Now allows for one-to-many matches for the main matching procedure (and imputing potential outomces). The user may specify `M_matches >=1`.
  - Now allows for one-to-many matches for the variance component when `J_var_matches >1`
  - Allowed for more user-specified arguments (for fitting PS models) i.e. `model_options`
- Divergence: Using `multiMatch()` with `match_on='existing'` does not always  return the same results as using `multilevelGPSMatch()` for matching on the existing (user-specified) generalized propensity scores. 
- Added S3 methods for `print` and `summary` for the `multiMatch` class
- Added Brian Barkley (@BarkleyBG) as co-author and maintainer
 

# multilevelMatching 0.1.0

- Original package version, released to GitHub by Shu Yang (@shuyang1987)
