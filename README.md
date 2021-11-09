# Interference Pattern Fitting for Slit Diffraction
Python class to analyze diffraction slit experimental data.

TODO:
- add peak/min detection for calculating observed wavelengths automatically
- add possible cutoff parameter: this is the far field diffraction limit after all
- add more customizable plotting options, title, colors etc

# Output
After applying to the example file in this repository, we have the following outputs.
Output plot:
![slitgraph](https://user-images.githubusercontent.com/42879357/140615957-8e085456-5bc2-494d-887e-b81359aacf21.PNG)
Output data:
```
=========
Observed peak intensity: 0.169
Fitted peak intensity: 0.1759691987726578 +/- 6.584519466019747e-05
=========
True slit spacing: 0.00025
Fitted slit spacing: 0.00024726060262317 +/- 4.746176627367676e-08
=========
True slit width: 4e-05
Fitted slit width: 3.43003371250152e-05 +/- 5.96292694103109e-08
=========
```
