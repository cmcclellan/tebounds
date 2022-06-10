## Treatment Effect Bounds

This is a R port of the Stata commmand 'tebounds'.  It calculates the average treatment effect bounds for a binary outcome and binary treatment under endogeniety and misreporting.  It implements three sets of assumptions to tighten worst case bounds: montone treatment selection, monotone treatment response, and monotone instrumental variables.  Options and functionality adhere as closely as possible to Stata's 'tebounds' and further information can be found in that supporting documentation.<sup>1</sup>

### Installation
To install
```markdown
library(devtools)
install_github("cmcclellan/tebounds")
```
### References
<sup>1</sup>McCarthy, Ian, Daniel L. Millimet, and Manan Roy. "Bounding treatment effects: A command for the partial identification of the average treatment effect with endogenous and misreported treatment assignment." The Stata Journal 15.2 (2015): 411-4
