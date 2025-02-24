import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import pandas as pd
import numpy as np
from scipy.stats import t
from sklearn.preprocessing import MinMaxScaler

#downloading and sclaing data - min max - for sigmoid function in nn
 

#scale data
scaler = MinMaxScaler()
data_scaled = scaler.fit_transform(data)
data_scaled
data_scaled = pd.DataFrame(data_scaled)
data_scaled

# Step 1: Create the mask for NaN values
nan_mask = ~np.isnan(data_scaled.values)  # True for valid values, False for NaNs

# Replace NaNs with 0 only temporarily for input; they will be ignored in loss
data_filled = np.nan_to_num(data_scaled.values, nan=0)  

# Convert data and mask to PyTorch tensors
data_tensor = torch.tensor(data_filled, dtype=torch.float32)
mask_tensor = torch.tensor(nan_mask, dtype=torch.float32)

# Step 2: Update the DataLoader to include the mask
dataset = TensorDataset(data_tensor, mask_tensor)
dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)



# Define the Autoencoder
class Autoencoder(nn.Module):
    def __init__(self, input_dim, latent_dim):
        super(Autoencoder, self).__init__()
        # Encoder
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 64),
            nn.ReLU(),
            nn.Linear(64, latent_dim),
            nn.ReLU()
        )
        # Decoder
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, 64),
            nn.ReLU(),
            nn.Linear(64, input_dim),
            nn.Sigmoid()  # Output values in [0, 1]
        )
    
    def forward(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return decoded

# Model Parameters
input_dim = data.shape[1]  # Number of proteins
latent_dim = 10  # Dimensionality of the latent space
model = Autoencoder(input_dim, latent_dim)

# Define Loss and Optimizer
criterion = nn.MSELoss(reduction='none')  # Reconstruction loss
optimizer = optim.Adam(model.parameters(), lr=0.001)

# Train the Autoencoder
epochs = 50
for epoch in range(epochs):
    total_loss = 0
    for batch in dataloader:
        inputs, mask = batch
        optimizer.zero_grad()

        # Forward pass through the model
        outputs = model(inputs)

        # Calculate the MSE loss for each element
        loss = criterion(outputs, inputs)

        # Mask the loss: Only keep the loss for valid (non-NaN) values
        masked_loss = loss * mask  # Multiply the loss by the mask (zeros out loss for NaNs)

        # Average the loss over all valid values
        loss = masked_loss.sum() / mask.sum()  # Sum of valid loss values, then normalize by count of valid elements

        loss.backward()
        optimizer.step()

        total_loss += loss.item()
    
    print(f"Epoch {epoch + 1}/{epochs}, Loss: {total_loss / len(dataloader):.4f}")

# Save the trained model
torch.save(model.state_dict(), "autoencoder.pth")

protein_to_ko = 'PBRM1'

protein_index = data.columns.get_loc(protein_to_ko)
protein_index

# Knock out the protein by setting its values to NaN across all samples
knockout_data = data_tensor.clone()
knockout_data[:, protein_index] = float('nan')  # Set the selected protein to NaN across all samples

# Fill the NaN values in the knockout data (we'll fill them with 0 for the model)
knockout_data_filled = torch.nan_to_num(knockout_data, nan=0.0)

# Use the trained autoencoder to reconstruct the knockout data
model.eval()  # Set the model to evaluation mode
with torch.no_grad():
    reconstructed = model(knockout_data_filled)

# Convert original and reconstructed data to DataFrames for easier viewing
original_df = pd.DataFrame(knockout_data_filled.numpy(), columns=data.columns)
reconstructed_df = pd.DataFrame(reconstructed.numpy(), columns=data.columns)

# Show the first 10 rows for both
print("Original data (after knockout):")
print(original_df.head(10))  # Original data with the protein knocked out

print("\nReconstructed data (after knockout):")
print(reconstructed_df.head(10))  # Reconstructed data after knockout