# multilevelMatching 0.1.0.9000



- Added `multiMatch()` function to carry out all types of matching. This effectively combines `multilevelGPSMatch()` and `multilevelMatchX()` into one function.
- `multiMatch()` has a number of small improvements built in.
- New development on matching for the foreseeable future will be on the `multiMatch()` function.
- `multiMatch()` does not implement stratification, which is still carried out through `multilevelGPSStratification()` 
    - FWIW, `multilevelGPSStratification()` needs work and has not improved much since v0.1.
- Some subfunctions to carry out procedures have changed in `multiMatch()` from the `multilevelGPSMatch()` and `multilevelMatchX()`. 
    - For example, `prepareData()` is used in `multiMatch()`
    - whereas `prepareData_legacy()` and `estimateTau_legacy()` are used in `multilevelGPSMatch()` and `multilevelMatchX()`. 
- Breaking: Output has become tidier. See `estimateTau()` for details
- Now allows for one-to-many matches for the main matching procedure (and imputing potential outomces). The user may specify `M_matches >=1`
- Now allows for one-to-many matches for the variance component when `J_var_matches >1`
- Added `calcSigSqAI2006()` (and unit tests and defensive programming) to implement the one-to-many matching for $\hat{\sigma}^2(X_i, W_i)$ estimator as introduced in Abadie and Imbens 2006 Econometrica
- Added + exported `calcKMVarFactor()` function for a variance component. It passes unit tests.
- Allowed for more user-specified arguments (for fitting PS models) i.e. `model_options`
- Names/rownames from the `X`, `Y`, or `W` args should be handled and treated as identifying information for the study units, and passed on to some of the output information.
- Added stable references to functions from different namespaces/packages with `::`
- Cleaned up package dependencies
- `multilevelGPSStratification()` now sorts by treatment level `W`
- Vignette now illustrates how user can supply propsensity scores via `GPSM="existing"`
- Add matrix of imputed potential outcomes to the output
- Implemented some unit tests, checks, error-avoidances, defensive programming, etc. 


#### A note on matching on existing GPS:

- Removes warning for using `multiMatch()` with existing GPS 
    - (Fixes BarkleyBG/multilevelMatching/#4)
- Using `multilevelGPSMatch()` on existing GPS produces same output as original version 0.1.0
- Using `multiMatch()` with `match_on='existing'` does NOT ALWAYS return the same output as original version 0.1.0
   - It will not necessarily return the same output when two units have the same GPS for a treatment level, and ties need to be broken. This is because the matching procedure takes place in a different order than in `multilevelGPSMatch()`, and so random number generation may/will be different.
   - This is likely a problem for all different methods of matching (when ties need to be broken)
   - However, without ties, the two functions may return identical estimates
   - As such, I'm removing the warnings and proceeding as is.
- Adds some notes in the README and vignette about this issue.

 

# multilevelMatching 0.1.0

- Link: the [original version by Shu Yang](https://github.com/shuyang1987/multilevelMatching)
- Instructions for downloading original version available [here](README.md) and [here](https://github.com/shuyang1987/multilevelMatching/blob/master/README.md)
