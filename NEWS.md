
# multilevelMatching 1.0.0

- Original CRAN release version by Shu Yang ([@shuyang1987](https://github.com/shuyang1987/)) and Brian G Barkley ([@BarkleyBG](https://github.com/BarkleyBG/)).
- Added `multiMatch()` function to carry out all types of matching. This effectively combines `multilevelGPSMatch()` and `multilevelMatchX()` into one function. Features include:
  - Better output: 
     - tidier estimates: see `estimateTau()`
     - plus more information from under-the-hood: see the `impute_mat` object for matrix off all imputed potential outcomes.
  - Now allows for one-to-many matches for the main matching procedure (and imputing potential outomces). The user may specify `M_matches >=1`.
  - Now allows for one-to-many matches to estimate the variance component. The user may specify `J_var_matches >=1`.
  - Allowed for more user-specified arguments (for fitting PS models) i.e. `model_options`
- Divergence: Using `multiMatch()` with `match_on='existing'` does not always  return the same results as using `multilevelGPSMatch()` for matching on the existing (user-specified) generalized propensity scores. 
- Added S3 methods for `print` and `summary` for the `multiMatch` class
- Added S3 methods for `print` and `summary` for the `multiMatch` class
- Users can apply the `estimateTrtModel()` function before using `multiMatch()` to verify that the model fitted in `multiMatch()` is the same as the user desires

# multilevelMatching 0.1.0

- Original package version, released to GitHub by Shu Yang (@shuyang1987). See [Release v0.1.0](https://github.com/shuyang1987/multilevelMatching/releases/tag/v0.1.0).
