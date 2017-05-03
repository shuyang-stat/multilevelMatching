# multilevelMatching 0.1.1 (under development)

## This is a forked repo 

- [Original version by Shu Yang](https://github.com/shuyang1987/multilevelMatching)

## Breaking changes

- Output has become tidier. See `estimateTau()` for details

## Small improvements

- Implemented some unit tests
- Some checks and bug-avoidances are introduced
- Cleaned up package dependencies
- Added `estimateTau()` to unify the main 3 functions (Dont Repeat Yourself principle)
- `multilevelGPSStratification()` now sorts by treatment level `W`
- Vignette now illustrates how user can supply propsensity scores via `GPSM="existing"`
- Matrix of all imputed potential outcomes included in output
- Allowing for user to specify reference category of multinomial logistic model via `model_options`

## Planned improvements

- cleanup code w/ DRY principle
- perhaps merge main 3 user-facing functions into 1
- Allow for more user-specified arguments (for fitting PS models) i.e. `model_options`

# multilevelMatching 0.1

- Link to the [original version by Shu Yang](https://github.com/shuyang1987/multilevelMatching)
- Instructions for downloading original version available [here](README.md) and [here](https://github.com/shuyang1987/multilevelMatching/blob/master/README.md)
