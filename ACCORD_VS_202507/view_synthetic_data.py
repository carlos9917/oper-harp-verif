import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

# Load analysis and forecast fields
analysis = xr.open_dataset('synthetic_tp_analysis.nc')['tp_analysis'].values
fc_5 = xr.open_dataset('synthetic_tp_forecast_5.nc')['tp_forecast'].values
fc_10 = xr.open_dataset('synthetic_tp_forecast_10.nc')['tp_forecast'].values
fc_20 = xr.open_dataset('synthetic_tp_forecast_20.nc')['tp_forecast'].values
fc_30 = xr.open_dataset('synthetic_tp_forecast_30.nc')['tp_forecast'].values

fields = [analysis, fc_5, fc_10, fc_20, fc_30]
titles = ['Analysis', 'Forecast 5', 'Forecast 10', 'Forecast 20', 'Forecast 30']

plt.figure(figsize=(15, 3))
for i, (field, title) in enumerate(zip(fields, titles)):
    print(field.shape)
    plt.subplot(1, 5, i+1)
    plt.imshow(field.squeeze(), cmap='viridis')
    plt.title(title)
    plt.axis('off')
plt.colorbar()
plt.tight_layout()
plt.show()
