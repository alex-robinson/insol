# insol

Currently there is one public function available in the insol module: `calc_insol_day`. This function returns the average insolation for a given day of year, latitude, and time between 5 Ma years before present and 1 Ma into the future. The function can calculate the insolation for a point, a vector of latitudes, or a 2D array of latitudes, returning a double with the same size as the input. 

```fortran
double precision function calc_insol_day(day,lats,time_bp,[S0],[day_year],[fldr])
```

```
day        Day of the year
lats       point, vector or 2D array of latitude values of interest 
           (-90:90 degrees)
time_bp    "time before present", ie, year relative to 1950
S0         Solar constant (optional, default 1365.0 W/m2)
day_year   Total number of days in the year, currently only
           handles day_year=360 (optional)
fldr       Folder containing table of precomputed orbital 
           parameter values (optional, default=input/)
```

Additional interfaces to be developed may include:

```fortran
double precision function calc_insol_daily(days,lat,time_bp,[S0],[day_year],[fldr])
```
... to calculate daily insolation over a year for a given latitude.

```fortran
double precision function calc_insol_hour(hour,day,lats,time_bp,[S0],[day_year],[fldr])
```
... to calculate the hourly insolation for a given day and latitude(s).

More?