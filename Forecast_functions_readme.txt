To make the algorithm work, there are five variables and two functions that you might want to change in the file "House_thermal_model_GA.py":

-------------------------------
Variable		| Line
-------------------------------
Start day		| 597
End day			| 598
Prediction horizon	| 599
Forecast method		| 605
Date/time list format	| 606
-------------------------------

For the forecast method and date/time list format, there is a comment with the options available. 
The functions I made for the prediction are:

- date_time_list (line 303): creates a list with the time signals. There are two formats, a single value, referenced to the 1st of January, and a list with the format [day (referenced to the 1st of January), hour, minute].
- TS_forecast (line 321): creates a forecast for a time series given. This might need some update later to determine the parameters for each parameter (temperature, PV, load...).

The optimization is executed in line 715, in case you want to review the main logic.