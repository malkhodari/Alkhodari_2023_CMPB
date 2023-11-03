# Circadian assessment of heart failure using explainable deep learning and novel multi-parameter polar images

This repository includes the main codes needed to run the work presented at "Alkhodari et al. (2022) Circadian assessment of heart failure using explainable deep learning and novel multi-parameter polar images. Computer Methods and Programs in Biomedicine. https://doi.org/XXX".

The "methodology" includes the following,
1) Load ECG and clinical information in .m and .csv, respectively.
2) Extract HRV using the Pan-Tompkins algorithm, filter HRV, extract hourly HRV, and adjust starting hour to 00:00.
3) Extract HRV features in time, frequency, non-linear, and fragmentation.
4) Generate HRV images with filled-in color-coded clinical information.
5) Predict heart failure stage using the pre-trained deep learning model (model structure is also available).
6) Interpret the predictions using GRAD-CAM, most important time hours, and most important clinical information.

For any questions, please do not hesitate to contact the corresponding author (Mohanad Alkhodari).
