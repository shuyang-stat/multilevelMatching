# multilevelMatching 0.2.2

- Now allows for one-to-many matches for the variance component when `J_var_matches >1`
- Added `calcSigSqAI2006()` (and unit tests and defensive programming) to implement the one-to-many matching for $\hat{\sigma}^2(X_i, W_i)$ estimator as introduced in Abadie and Imbens 2006 Econometrica

## Package information

- Author and maintainer of the package is Shu Yang
- Improvements from version 0.1 to version 0.2.2 were contributed by Brian Barkley
- The original version (v0.1) can be accessed from:
    - forked repo: https://github.com/BarkleyBG/multilevelMatching
    - original repo: https://github.com/shuyang1987/multilevelMatching
    
    
## Planned improvements

See [GH issues](https://github.com/BarkleyBG/multilevelMatching/issues)

# multilevelMatching 0.2.1

- Added + exported `calcKMFactor()` function for a variance component. It passes unit tests.

    
# multilevelMatching 0.2

## New development for v0.2: the `multiMatch()` function

- Added `multiMatch()` function to carry out all types of matching. This effectively combines `multilevelGPSMatch()` and `multilevelMatchX()` into one function.
- `multiMatch()` has a number of small improvements built in.
- New development on matching for the foreseeable future will be on the `multiMatch()` function.
- `multiMatch()` does not implement stratification, which is still carried out through `multilevelGPSStratification()` 
    - FWIW, `multilevelGPSStratification()` needs work and has not improved much since v0.1.
- Some subfunctions to carry out procedures have changed in `multiMatch()` from the `multilevelGPSMatch()` and `multilevelMatchX()`. 
    - For example, `prepareData()` is used in `multiMatch()`
    - whereas `prepareData_legacy()` and `estimateTau_legacy()` are used in `multilevelGPSMatch()` and `multilevelMatchX()`. 

## Other improvements

- Cleaned code w/ DRY principle
- Allowed for more user-specified arguments (for fitting PS models) i.e. `model_options`
- Names/rownames from the `X`, `Y`, or `W` args should be handled and treated as identifying information for the study units, and passed on to some of the output information.


# multilevelMatching 0.1.5

- Added stable references to functions with `::`

# multilevelMatching 0.1.4

- Cleaned up package dependencies
- `multilevelGPSStratification()` now sorts by treatment level `W`
- Vignette now illustrates how user can supply propsensity scores via `GPSM="existing"`

# multilevelMatching 0.1.3

- Better sorting of dataset
- Sorting the matrix of imputed potential outcomes

# multilevelMatching 0.1.2
 
- Add matrix of imputed potential outcomes to output

# multilevelMatching 0.1.1

## Breaking changes

- Output has become tidier. See `estimateTau()` for details

## Small improvements

- Implemented some unit tests
- Some checks and bug-avoidances are introduced
- Added `estimateTau()` to unify the main 3 functions (Dont Repeat Yourself principle)


# multilevelMatching 0.1

- Link: the [original version by Shu Yang](https://github.com/shuyang1987/multilevelMatching)
- Instructions for downloading original version available [here](README.md) and [here](https://github.com/shuyang1987/multilevelMatching/blob/master/README.md)
