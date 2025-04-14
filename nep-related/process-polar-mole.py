import numpy as np
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt

def calculate_metrics_and_save_plots_molecular(predicted_values, actual_values, suffix):
    """
    Calculate RMSE and R2 metrics for predicted vs. actual molecular polarizability values
    and save the diagonal plots as images.
    
    Args:
    - predicted_values: Numpy array of predicted polarizability values.
    - actual_values: Numpy array of actual (DFT) polarizability values.
    - suffix: Suffix for the output plot filenames to differentiate them.
    """
    # Split the data for axx ayy azz and axy axz ayz
    axx_ayy_azz_pred = predicted_values[:, :3].flatten()
    axx_ayy_azz_actual = actual_values[:, :3].flatten()
    axy_axz_ayz_pred = predicted_values[:, 3:].flatten()
    axy_axz_ayz_actual = actual_values[:, 3:].flatten()

    # Metrics for axx ayy azz
    rmse_axx_ayy_azz = np.sqrt(mean_squared_error(axx_ayy_azz_actual, axx_ayy_azz_pred))
    r2_axx_ayy_azz = r2_score(axx_ayy_azz_actual, axx_ayy_azz_pred)

    # Metrics for axy axz ayz
    rmse_axy_axz_ayz = np.sqrt(mean_squared_error(axy_axz_ayz_actual, axy_axz_ayz_pred))
    r2_axy_axz_ayz = r2_score(axy_axz_ayz_actual, axy_axz_ayz_pred)

    # Plot for axx ayy azz
    plt.figure(figsize=(6, 6))
    plt.scatter(axx_ayy_azz_actual, axx_ayy_azz_pred, alpha=0.5)
    min_val = min(axx_ayy_azz_actual.min(), axx_ayy_azz_pred.min()) -0.2
    max_val = max(axx_ayy_azz_actual.max(), axx_ayy_azz_pred.max()) +0.2
    plt.xlim(min_val, max_val)
    plt.ylim(min_val, max_val)
    plt.plot([min_val, max_val], [min_val, max_val], 'k--', lw=2) 
    #plt.plot([axx_ayy_azz_actual.min(), axx_ayy_azz_actual.max()], 
             #[axx_ayy_azz_actual.min(), axx_ayy_azz_actual.max()], 'k--', lw=2)
    plt.xlabel('Actual Values (axx ayy azz)')
    plt.ylabel('Predicted Values (axx ayy azz)')
    plt.title('Diagonal Plot for axx ayy azz') 
    axx_ayy_azz_plot_path = f'plot_mole_diag_{suffix}.jpg'
    plt.savefig(axx_ayy_azz_plot_path, format='jpg')
    plt.close()

    # Plot for axy axz ayz
    plt.figure(figsize=(6, 6))
    plt.scatter(axy_axz_ayz_actual, axy_axz_ayz_pred, alpha=0.5)
    min_val = min(axy_axz_ayz_actual.min(), axy_axz_ayz_pred.min()) -0.2
    max_val = max(axy_axz_ayz_actual.max(), axy_axz_ayz_pred.max()) +0.2
    plt.xlim(min_val, max_val)
    plt.ylim(min_val, max_val)
    plt.plot([min_val, max_val], [min_val, max_val], 'k--', lw=2) 
    #plt.plot([axy_axz_ayz_actual.min(), axy_axz_ayz_actual.max()], 
             #[axy_axz_ayz_actual.min(), axy_axz_ayz_actual.max()], 'k--', lw=2)
    plt.xlabel('Actual Values (axy axz ayz)')
    plt.ylabel('Predicted Values (axy axz ayz)')
    plt.title('Diagonal Plot for axy axz ayz')
    plt.axis('equal') 
    axy_axz_ayz_plot_path = f'plot_mole_offdiag_{suffix}.jpg'
    plt.savefig(axy_axz_ayz_plot_path, format='jpg')
    plt.close()

    return {
        "RMSE_axx_ayy_azz": rmse_axx_ayy_azz,
        "R2_axx_ayy_azz": r2_axx_ayy_azz,
        "RMSE_axy_axz_ayz": rmse_axy_axz_ayz,
        "R2_axy_axz_ayz": r2_axy_axz_ayz,

    }

# Load the data
file_path = 'polarizability_train.out'
data = np.loadtxt(file_path)

# Split the data into predicted and actual values
predicted_values = data[:, :6]  # First half of each row
actual_values = data[:, 6:]  # Second half of each row

# Calculate metrics and save plots
metrics_and_plots = calculate_metrics_and_save_plots_molecular(predicted_values, actual_values, "molecular")

# Print metrics and plot paths for download
print(metrics_and_plots)
