import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset

# Define the Conditional VAE
class ConditionalVAE(nn.Module):
    def __init__(self, input_dim, cond_dim, hidden_dim, latent_dim):
        """
        input_dim: Number of protein features.
        cond_dim: Dimensionality of the condition vector (e.g., same as input_dim if you use a binary mask).
        hidden_dim: Number of hidden units.
        latent_dim: Size of the latent space.
        """
        super(ConditionalVAE, self).__init__()
        self.input_dim = input_dim
        self.cond_dim = cond_dim
        
        # Encoder: input is concatenation of data and condition.
        self.fc1 = nn.Linear(input_dim + cond_dim, hidden_dim)
        self.fc_mu = nn.Linear(hidden_dim, latent_dim)
        self.fc_logvar = nn.Linear(hidden_dim, latent_dim)
        
        # Decoder: latent vector concatenated with condition.
        self.fc3 = nn.Linear(latent_dim + cond_dim, hidden_dim)
        self.fc4 = nn.Linear(hidden_dim, input_dim)
    
    def encode(self, x, c):
        # Concatenate the protein data and the condition vector.
        x_cond = torch.cat([x, c], dim=1)
        h = F.relu(self.fc1(x_cond))
        mu = self.fc_mu(h)
        logvar = self.fc_logvar(h)
        return mu, logvar
    
    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std
    
    def decode(self, z, c):
        # Concatenate the latent code and the condition vector.
        z_cond = torch.cat([z, c], dim=1)
        h = F.relu(self.fc3(z_cond))
        x_recon = self.fc4(h)
        return x_recon
    
    def forward(self, x, c):
        mu, logvar = self.encode(x, c)
        z = self.reparameterize(mu, logvar)
        x_recon = self.decode(z, c)
        return x_recon, mu, logvar

# Loss function: Reconstruction loss (MSE for continuous data) + KL divergence
def CVAE_loss_function(x_recon, x, mu, logvar):
    # Reconstruction loss: Mean Squared Error summed over features.
    recon_loss = F.mse_loss(x_recon, x, reduction='sum')
    # KL divergence loss
    kl_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
    return recon_loss + kl_loss

