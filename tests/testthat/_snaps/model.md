# invalid input gives an error

    Code
      run_basin_model(bad)
    Condition
      Error:
      ! phi0 must lie between 0 and 1

---

    Code
      run_basin_model()
    Condition
      Error:
      ! no stratigraphic units supplied (use 'site' or 'units')

